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

__version__ = "1.0.0"

OUTPUT_COLUMNS = [
    "Sample", "tox_gene", "tox_allele", "ST", "mlst_profile",
    "Identity", "Coverage", "db_version", "notes",
]

BLAST_FIELDS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                "qlen", "slen"]


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
        return f"~{b.allele}", None
    if b.identity >= min_id:
        return f"?{b.allele}", None
    return "-", None


def call_st(db, calls):
    if any(not c.isdigit() for c in calls):
        return "NA"
    return db.profiles.get(tuple(calls), "novel")


def call_tox(hits, min_id, min_cov):
    """Return dict with tox_gene, tox_allele, Identity, Coverage, note."""
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
        out["note"] = f"closest tox allele {b.allele}; confirm novel allele"
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
    if row["ST"] == "novel":
        notes.append("exact alleles at all loci but profile not in database")
    elif row["ST"] == "NA":
        missing = [l for l, c in zip(db.loci, calls) if c == "-"]
        if len(missing) == len(db.loci):
            notes.append("no MLST loci found - check species/assembly")

    tox = call_tox(hits.get("tox", []), args.tox_min_identity,
                   args.tox_min_coverage)
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
        w = csv.DictWriter(out, fieldnames=OUTPUT_COLUMNS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    finally:
        if out is not sys.stdout:
            out.close()


if __name__ == "__main__":
    main()
