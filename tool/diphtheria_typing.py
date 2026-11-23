#!/usr/bin/env python3
"""
diphtheria_typing.py - tox detection, tox allele assignment and 7-gene MLST
for Corynebacterium diphtheriae genome assemblies.

Database layout (see database/README in the repository):

    <db>/VERSION.txt                 single line, e.g. Pasteur_2026-10
    <db>/tox/tox_alleles.fasta       >tox_<allele_id>
    <db>/mlst/profiles.tsv           ST <locus1> ... <locusN> [extra cols]
    <db>/mlst/alleles/<locus>.fasta  ><locus>_<allele_id>

Method: every allele (MLST loci + tox) is searched against each assembly
with blastn. Exact = 100% identity over 100% of the allele with no gaps.

Only dependency outside the standard library is BLAST+ (blastn, makeblastdb).
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

__version__ = "1.2.0"

OUTPUT_COLUMNS = [
    "Sample", "tox_gene", "tox_allele", "ST", "mlst_profile",
    "Identity", "Coverage", "db_version", "notes",
]


def output_columns(loci, layout="combined"):
    """Column order for the chosen allele layout.

    combined: one mlst_profile column (atpA(1);dnaE(1);...), the default.
    separate: one column per locus, ALLELE1..ALLELEn (n = number of loci,
    7 for this scheme), each holding one call such as atpA(1).
    """
    if layout != "separate":
        return list(OUTPUT_COLUMNS)
    i = OUTPUT_COLUMNS.index("mlst_profile")
    alleles = [f"ALLELE{k + 1}" for k in range(len(loci))]
    return OUTPUT_COLUMNS[:i] + alleles + OUTPUT_COLUMNS[i + 1:]

BLAST_FIELDS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                "qlen", "slen", "qseq", "sseq"]

# Stop codons, for judging whether a tox gene is intact.
STOP_CODONS = ("TAA", "TAG", "TGA")


# --------------------------------------------------------------------------
# Database
# --------------------------------------------------------------------------

class TypingDB:
    def __init__(self, path):
        self.path = path
        self.version = self._read_version()
        self.tox_fasta = os.path.join(path, "tox", "tox_alleles.fasta")
        self.profiles_tsv = os.path.join(path, "mlst", "profiles.tsv")
        self.alleles_dir = os.path.join(path, "mlst", "alleles")
        for p in (self.tox_fasta, self.profiles_tsv, self.alleles_dir):
            if not os.path.exists(p):
                sys.exit(f"ERROR: database is incomplete, missing: {p}")
        self.loci, self.profiles = self._read_profiles()
        for locus in self.loci:
            if not os.path.exists(self.locus_fasta(locus)):
                sys.exit(f"ERROR: no allele file for MLST locus {locus}")
        self.disrupted_tox = self._read_disrupted_tox()

    def _read_version(self):
        vf = os.path.join(self.path, "VERSION.txt")
        if not os.path.exists(vf):
            return "unknown"
        with open(vf) as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    return line
        return "unknown"

    def locus_fasta(self, locus):
        return os.path.join(self.alleles_dir, f"{locus}.fasta")

    def _read_disrupted_tox(self):
        """{allele: first stop codon} for bundled tox alleles that are
        themselves disrupted - a stop codon before the end of the gene.

        Read from the database rather than hardcoded, so it stays correct when
        the Pasteur snapshot is updated. In Pasteur_2026-10b, 14 of the 76 tox
        alleles are frameshifted: real tox-bearing but non-toxin-producing
        variants. A near match to one of those must not be reported as a gene
        this isolate has newly broken.
        """
        disrupted = {}
        for seqid, seq in read_fasta(self.tox_fasta):
            _, allele = split_allele_id(seqid)
            status, stop, _ = cds_status(seq)
            if status == "disrupted":
                disrupted[allele] = stop
        return disrupted

    def _read_profiles(self):
        with open(self.profiles_tsv, newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            if header[0].upper() != "ST":
                sys.exit("ERROR: profiles.tsv must start with an 'ST' column")
            # Locus columns are those with a matching allele file; anything
            # else (clonal_complex, LINcode, ...) is ignored.
            loci = [h for h in header[1:]
                    if os.path.exists(os.path.join(self.alleles_dir, f"{h}.fasta"))]
            idx = [header.index(l) for l in loci]
            profiles = {}
            for row in reader:
                if not row or not row[0].strip():
                    continue
                key = tuple(row[i].strip() for i in idx)
                profiles[key] = row[0].strip()
        if not loci:
            sys.exit("ERROR: no MLST loci found in profiles.tsv header")
        return loci, profiles


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0] if len(line) > 1 else "", []
            elif line:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def validate_fasta(path):
    """Return (n_contigs, total_length) or raise ValueError."""
    n, total = 0, 0
    valid = set("ACGTNRYKMSWBDHVacgtnrykmswbdhv-")
    for name, seq in read_fasta(path):
        n += 1
        total += len(seq)
        bad = set(seq) - valid
        if bad:
            raise ValueError(f"non-nucleotide characters in {name}: "
                             f"{''.join(sorted(bad))[:10]}")
    if n == 0 or total == 0:
        raise ValueError("no sequences found (is this a FASTA file?)")
    return n, total


def split_allele_id(seqid):
    """'atpA_12' -> ('atpA', '12'); 'tox_3' -> ('tox', '3')."""
    if "_" not in seqid:
        return seqid, seqid
    locus, allele = seqid.rsplit("_", 1)
    return locus, allele


def first_stop(seq):
    """1-based codon number of the first stop codon, or None if there is none.

    Reads seq in frame from its first base, so the caller must pass a
    sequence that starts at the gene's start codon.
    """
    seq = seq.upper()
    for i in range(0, len(seq) - 2, 3):
        if seq[i:i + 3] in STOP_CODONS:
            return i // 3 + 1
    return None


def cds_status(seq):
    """Judge a whole gene sequence read from its own start codon.

    Returns (status, stop_codon, total_codons). status is "intact" when the
    only stop codon is the last one, else "disrupted" - a stop codon before
    the end (whether from a nonsense change or from a frameshift), or no stop
    codon at all. stop_codon is None in that last case.

    This looks only at the sequence given, never at another allele, so it is
    unaffected by whether the closest known allele happens to be in frame.
    """
    total = len(seq) // 3
    stop = first_stop(seq)
    if stop is not None and stop == total and len(seq) % 3 == 0:
        return "intact", stop, total
    return "disrupted", stop, total


FASTA_EXTS = (".fasta", ".fa", ".fna", ".fas", ".fsa", ".contigs", ".scaffolds")


def clean_name(name):
    """'CD001.fasta' -> 'CD001'; also drops .gz and repeated extensions."""
    base = name.strip()
    changed = True
    while changed:
        changed = False
        for ext in (".gz",) + FASTA_EXTS:
            if base.lower().endswith(ext) and len(base) > len(ext):
                base = base[:-len(ext)]
                changed = True
    return base or name


def run(cmd):
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True)
    if res.returncode != 0:
        sys.exit(f"ERROR running {' '.join(cmd)}:\n{res.stderr}")
    return res.stdout


class Hit:
    __slots__ = BLAST_FIELDS + ["locus", "allele"]

    def __init__(self, fields):
        for k, v in zip(BLAST_FIELDS, fields):
            setattr(self, k, v)
        for k in ("pident", "evalue", "bitscore"):
            setattr(self, k, float(getattr(self, k)))
        for k in ("length", "mismatch", "gapopen", "qstart", "qend",
                  "sstart", "send", "qlen", "slen"):
            setattr(self, k, int(getattr(self, k)))
        self.locus, self.allele = split_allele_id(self.qseqid)

    @property
    def coverage(self):
        """% of the allele covered by the alignment."""
        return 100.0 * (self.qend - self.qstart + 1) / self.qlen

    @property
    def identity(self):
        return self.pident

    @property
    def exact(self):
        return (self.pident == 100.0 and self.gapopen == 0
                and self.qstart == 1 and self.qend == self.qlen
                and self.length == self.qlen)

    @property
    def at_contig_edge(self):
        lo, hi = sorted((self.sstart, self.send))
        return lo <= 1 or hi >= self.slen


def diff_counts(hit):
    """(substitutions, bases inserted, bases deleted) over the aligned region.

    The allele is the query and the assembly the subject, so a gap in qseq is
    a base the assembly has gained and a gap in sseq one it has lost.
    """
    snps = inserted = deleted = 0
    for q, s in zip(hit.qseq.upper(), hit.sseq.upper()):
        if q == "-":
            inserted += 1
        elif s == "-":
            deleted += 1
        elif q != s:
            snps += 1
    return snps, inserted, deleted


def diff_summary(hit):
    """Differences from the allele in words: '1 SNP', '2 SNPs, 3 bp deletion'."""
    snps, inserted, deleted = diff_counts(hit)
    parts = []
    if snps:
        parts.append(f"{snps} SNP{'s' if snps != 1 else ''}")
    if deleted:
        parts.append(f"{deleted} bp deletion")
    if inserted:
        parts.append(f"{inserted} bp insertion")
    return ", ".join(parts) or "no differences"


def full_length(hit):
    """True when the alignment covers the whole allele, end to end."""
    return hit.qstart == 1 and hit.qend == hit.qlen


def blast_assembly(assembly, query_fasta, workdir, threads):
    dbprefix = os.path.join(workdir, "asm")
    run(["makeblastdb", "-in", assembly, "-dbtype", "nucl",
         "-out", dbprefix])
    out = run(["blastn", "-task", "blastn", "-query", query_fasta,
               "-db", dbprefix, "-outfmt", "6 " + " ".join(BLAST_FIELDS),
               "-evalue", "1e-20", "-perc_identity", "80",
               "-dust", "no", "-soft_masking", "false",
               "-max_target_seqs", "100000", "-max_hsps", "1",
               "-num_threads", str(threads)])
    hits = defaultdict(list)
    for line in out.splitlines():
        if line.strip():
            h = Hit(line.split("\t"))
            hits[h.locus].append(h)
    return hits


def best(hits):
    return max(hits, key=lambda h: (h.bitscore, h.identity, h.qlen))


# --------------------------------------------------------------------------
# Calling
# --------------------------------------------------------------------------

def call_mlst_locus(hits, min_id, min_cov):
    """Return (allele_call, note). Conventions follow tseemann/mlst:
       N   exact match
       ~N  novel allele, closest is N
       ?N  partial match to N (e.g. contig break)
       -   locus not found
       N,M multiple distinct exact copies"""
    if not hits:
        return "-", None
    exact = [h for h in hits if h.exact]
    if exact:
        # Distinct contig locations carrying exact alleles
        locs = {}
        for h in sorted(exact, key=lambda h: -h.qlen):
            key = (h.sseqid, min(h.sstart, h.send) // 50)
            locs.setdefault(key, h)
        alleles = sorted({h.allele for h in locs.values()},
                         key=lambda a: (len(a), a))
        if len(alleles) > 1:
            return ",".join(alleles), "multiple exact alleles"
        return alleles[0], None
    b = best(hits)
    if b.coverage >= min_cov and b.identity >= min_id:
        return f"~{b.allele}", f"{diff_summary(b)} from allele {b.allele}"
    if b.identity >= min_id:
        return f"?{b.allele}", (f"partial match to allele {b.allele} "
                                f"({b.coverage:.1f}% aligned)")
    return "-", None


def call_st(db, calls):
    if any(not c.isdigit() for c in calls):
        return "NA"
    return db.profiles.get(tuple(calls), "novel")


def tox_novel_note(hit, disrupted_tox):
    """The note for a novel tox allele: how it differs from the closest known
    allele, and whether the gene the assembly carries looks intact.

    Gene integrity is judged on the found sequence alone (cds_status), not by
    comparing indel lengths with the closest allele: 14 of the bundled alleles
    are themselves frameshifted, so such a comparison would both blame a
    reported SNP for a disruption the allele already had and, for a change that
    restores the frame, call an intact gene broken.
    """
    head = f"closest tox allele {hit.allele} ({hit.identity:.1f}%"
    if hit.coverage < 100:
        head += f", {hit.coverage:.1f}% aligned"
    note = f"{head}): {diff_summary(hit)}"

    if full_length(hit):
        found = hit.sseq.replace("-", "")
        status, stop, total = cds_status(found)
        if status == "disrupted":
            if stop is None:
                where = "no stop codon found"
            elif len(found) % 3:
                # The codon count of a frameshifted sequence means nothing,
                # so name the cause instead of quoting a total.
                where = f"frameshift, premature stop at codon {stop}"
            else:
                where = f"premature stop at codon {stop} of {total}"
            note += f" - {where}, possibly disrupted gene"
    else:
        note += "; gene not full length, integrity not checked"

    if hit.allele in disrupted_tox:
        note += (f"; allele {hit.allele} is itself a disrupted tox allele "
                 f"(stop at codon {disrupted_tox[hit.allele]})")
    return note


def call_tox(hits, min_id, min_cov, disrupted_tox=None):
    """Return dict with tox_gene, tox_allele, Identity, Coverage, note."""
    disrupted_tox = disrupted_tox or {}
    out = {"tox_gene": "tox not present", "tox_allele": "NA",
           "Identity": 0.0, "Coverage": 0.0, "note": None}
    if not hits:
        return out
    exact = [h for h in hits if h.exact]
    if exact:
        alleles = sorted({h.allele for h in exact}, key=lambda a: (len(a), a))
        out.update(tox_gene="tox", tox_allele=",".join(alleles),
                   Identity=100.0, Coverage=100.0)
        if len(alleles) > 1:
            out["note"] = "multiple exact tox alleles"
        return out
    b = best(hits)
    out["Identity"] = round(b.identity, 2)
    out["Coverage"] = round(b.coverage, 2)
    if b.identity >= min_id and b.coverage >= min_cov:
        out.update(tox_gene="novel tox allele", tox_allele=f"~{b.allele}")
        out["note"] = tox_novel_note(b, disrupted_tox)
        return out
    # Something tox-like is there but fails thresholds: never silently
    # report this as absent.
    if b.coverage >= 30:
        out.update(tox_gene="review", tox_allele=f"?{b.allele}")
        reason = []
        if b.coverage < min_cov:
            reason.append("partial coverage")
            if b.at_contig_edge:
                reason.append("hit at contig edge (possible assembly break)")
        if b.identity < min_id:
            reason.append("identity below threshold")
        out["note"] = "tox-like hit: " + "; ".join(reason)
        return out
    out["note"] = (f"weak tox-like fragment only "
                   f"({b.coverage:.1f}% cov, {b.identity:.1f}% id)")
    return out


def type_sample(name, assembly, db, query_fasta, args):
    row = {c: "" for c in OUTPUT_COLUMNS}
    row.update(Sample=name, db_version=db.version)
    notes = []
    try:
        n_contigs, total = validate_fasta(assembly)
    except ValueError as e:
        row.update(tox_gene="ERROR", tox_allele="NA", ST="NA",
                   Identity="NA", Coverage="NA", notes=f"invalid FASTA: {e}")
        return row
    if total < 1_000_000:
        notes.append(f"assembly only {total:,} bp - check completeness")

    work = tempfile.mkdtemp(prefix="diphtyp_", dir=args.tmpdir)
    try:
        hits = blast_assembly(assembly, query_fasta, work, args.threads)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    calls = []
    for locus in db.loci:
        c, note = call_mlst_locus(hits.get(locus, []),
                                  args.mlst_min_identity,
                                  args.mlst_min_coverage)
        calls.append(c)
        if note:
            notes.append(f"{locus}: {note}")
    row["ST"] = call_st(db, calls)
    row["mlst_profile"] = ";".join(f"{l}({c})" for l, c in zip(db.loci, calls))
    for k, (l, c) in enumerate(zip(db.loci, calls)):
        row[f"ALLELE{k + 1}"] = f"{l}({c})"
    if row["ST"] == "novel":
        notes.append("exact alleles at all loci but profile not in database")
    elif row["ST"] == "NA":
        missing = [l for l, c in zip(db.loci, calls) if c == "-"]
        if len(missing) == len(db.loci):
            notes.append("no MLST loci found - check species/assembly")

    tox = call_tox(hits.get("tox", []), args.tox_min_identity,
                   args.tox_min_coverage, db.disrupted_tox)
    row.update(tox_gene=tox["tox_gene"], tox_allele=tox["tox_allele"],
               Identity=f"{tox['Identity']:.1f}",
               Coverage=f"{tox['Coverage']:.1f}")
    if tox["note"]:
        notes.append(tox["note"])
    row["notes"] = "; ".join(notes) if notes else ""
    return row


def build_query(db, workdir):
    """Concatenate tox + MLST alleles into one query FASTA."""
    path = os.path.join(workdir, "all_alleles.fasta")
    with open(path, "w") as out:
        sources = [db.tox_fasta] + [db.locus_fasta(l) for l in db.loci]
        for src in sources:
            for name, seq in read_fasta(src):
                out.write(f">{name}\n{seq}\n")
    return path


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="tox detection, tox allele and 7-gene MLST typing for "
                    "Corynebacterium diphtheriae assemblies.")
    p.add_argument("-i", "--input", action="append", required=True,
                   help="Assembly FASTA (repeatable)")
    p.add_argument("-n", "--name", action="append",
                   help="Sample name for the matching --input (repeatable; "
                        "default: file name without extension)")
    p.add_argument("-d", "--db", required=True, help="Database directory")
    p.add_argument("-o", "--output", default="-", help="Output TSV (default stdout)")
    p.add_argument("--tox-min-identity", type=float, default=99.0)
    p.add_argument("--tox-min-coverage", type=float, default=95.0)
    p.add_argument("--mlst-min-identity", type=float, default=95.0)
    p.add_argument("--mlst-min-coverage", type=float, default=95.0)
    p.add_argument("--allele-layout", choices=["combined", "separate"],
                   default="combined",
                   help="combined (default): all MLST calls in one "
                        "mlst_profile column; separate: one column per "
                        "locus, ALLELE1..ALLELE7")
    p.add_argument("-t", "--threads", type=int, default=1)
    p.add_argument("--tmpdir", default=None)
    p.add_argument("--version", action="version", version=__version__)
    args = p.parse_args(argv)
    if args.name and len(args.name) != len(args.input):
        p.error("--name must be given once per --input")
    return args


def main(argv=None):
    args = parse_args(argv)
    for exe in ("blastn", "makeblastdb"):
        if shutil.which(exe) is None:
            sys.exit(f"ERROR: {exe} not found on PATH (install BLAST+)")
    db = TypingDB(args.db)
    names = [clean_name(n) for n in
             (args.name or [os.path.basename(f) for f in args.input])]

    qdir = tempfile.mkdtemp(prefix="diphtyp_q_", dir=args.tmpdir)
    try:
        query = build_query(db, qdir)
        rows = [type_sample(n, f, db, query, args)
                for n, f in zip(names, args.input)]
    finally:
        shutil.rmtree(qdir, ignore_errors=True)

    out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
    try:
        w = csv.DictWriter(out,
                           fieldnames=output_columns(db.loci, args.allele_layout),
                           delimiter="\t", lineterminator="\n",
                           extrasaction="ignore", restval="")
        w.writeheader()
        w.writerows(rows)
    finally:
        if out is not sys.stdout:
            out.close()


if __name__ == "__main__":
    main()
