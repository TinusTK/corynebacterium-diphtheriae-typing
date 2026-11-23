"""
Unit tests for diphtheria_typing.py using a SYNTHETIC database (random
sequences). These test the calling logic, not biology. Real-genome
validation belongs in docs/VALIDATION.md.

Run: python -m unittest discover -s tests -v   (needs BLAST+ on PATH)
"""
import csv
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "tool", "diphtheria_typing.py")
sys.path.insert(0, os.path.join(ROOT, "tool"))
import diphtheria_typing as dt            # noqa: E402  (needs sys.path above)
LOCI = ["atpA", "dnaE", "dnaK", "fusA", "leuA", "odhA", "rpoB"]
rng = random.Random(42)


def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, positions):
    s = list(seq)
    for p in positions:
        s[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[p]]
    return "".join(s)


STOPS = ("TAA", "TAG", "TGA")


def rand_orf(codons):
    """A plausible coding sequence: GTG, random sense codons, then TGA.

    tox is a real gene, and the tool now judges whether the gene an assembly
    carries is intact, so the synthetic tox allele has to be a genuine open
    reading frame rather than random bases.
    """
    body = []
    while len(body) < codons - 2:
        c = "".join(rng.choice("ACGT") for _ in range(3))
        if c not in STOPS:
            body.append(c)
    return "GTG" + "".join(body) + "TGA"


def mutate_in_frame(seq, positions):
    """Substitute at each position without creating a premature stop codon."""
    s = list(seq)
    for p in positions:
        start = (p // 3) * 3
        for alt in "ACGT":
            if alt == s[p]:
                continue
            trial = s[:]
            trial[p] = alt
            if "".join(trial[start:start + 3]) not in STOPS:
                s = trial
                break
        else:
            raise AssertionError(f"no silent-of-stop substitution at {p}")
    return "".join(s)


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


class TypingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if shutil.which("blastn") is None:
            raise unittest.SkipTest("BLAST+ not installed")
        cls.tmp = tempfile.mkdtemp()
        db = cls.db = os.path.join(cls.tmp, "db")
        os.makedirs(os.path.join(db, "tox"))
        os.makedirs(os.path.join(db, "mlst", "alleles"))
        with open(os.path.join(db, "VERSION.txt"), "w") as fh:
            fh.write("# test\nSYNTHETIC_TEST\n")

        # 3 alleles per locus, each a few SNPs apart
        cls.alleles = {}
        for locus in LOCI:
            base = rand_seq(450)
            cls.alleles[locus] = {"1": base,
                                  "2": mutate(base, [10, 200]),
                                  "3": mutate(base, [50, 300, 400])}
            write_fasta(os.path.join(db, "mlst", "alleles", f"{locus}.fasta"),
                        [(f"{locus}_{k}", v) for k, v in cls.alleles[locus].items()])
        # tox_1 and tox_2 are intact genes; tox_3 mirrors the 14 real Pasteur
        # alleles that are themselves frameshifted (here tox_1 with base 235
        # deleted, exactly as tox_24 relates to tox_1 in the real database).
        tox = rand_orf(561)
        cls.tox = {"1": tox,
                   "2": mutate_in_frame(tox, [100, 900, 1500]),
                   "3": tox[:234] + tox[235:]}
        write_fasta(os.path.join(db, "tox", "tox_alleles.fasta"),
                    [(f"tox_{k}", v) for k, v in cls.tox.items()])
        with open(os.path.join(db, "mlst", "profiles.tsv"), "w") as fh:
            fh.write("ST\t" + "\t".join(LOCI) + "\tclonal_complex\n")
            fh.write("8\t" + "\t".join(["1"] * 7) + "\tCC8\n")
            fh.write("26\t" + "\t".join(["2"] * 7) + "\t\n")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    # helpers ---------------------------------------------------------------
    def genome(self, name, profile, tox=None, split_tox=False, extra=None):
        """Build an assembly: background + MLST alleles + optional tox."""
        contigs, chunks = [], [rand_seq(3000)]
        for locus, a in zip(LOCI, profile):
            seq = a if len(a) > 5 else self.alleles[locus][a]
            chunks += [seq, rand_seq(800)]
        contigs.append(("contig1", "".join(chunks)))
        if tox and not split_tox:
            contigs.append(("contig2", rand_seq(2000) + tox + rand_seq(2000)))
        elif tox and split_tox:
            contigs.append(("contig2", rand_seq(2000) + tox[:900]))
            contigs.append(("contig3", tox[900:] + rand_seq(2000)))
        if extra:
            contigs += extra
        path = os.path.join(self.tmp, f"{name}.fasta")
        write_fasta(path, contigs)
        return path

    def run_tool(self, *files, extra_args=()):
        out = os.path.join(self.tmp, "out.tsv")
        cmd = [sys.executable, SCRIPT, "-d", self.db, "-o", out, *extra_args]
        for f in files:
            cmd += ["-i", f]
        subprocess.run(cmd, check=True)
        with open(out) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    # tests -----------------------------------------------------------------
    def test_columns_match_spec(self):
        rows = self.run_tool(self.genome("cols", ["1"] * 7))
        self.assertEqual(list(rows[0].keys()),
                         ["Sample", "tox_gene", "tox_allele", "ST",
                          "mlst_profile", "Identity", "Coverage",
                          "db_version", "notes"])
        self.assertEqual(rows[0]["db_version"], "SYNTHETIC_TEST")
        self.assertEqual(rows[0]["Sample"], "cols")

    def test_separate_allele_columns(self):
        rows = self.run_tool(self.genome("sep", ["1", "2", "1", "1", "1", "1", "1"]),
                             extra_args=["--allele-layout", "separate"])
        self.assertEqual(list(rows[0].keys()),
                         ["Sample", "tox_gene", "tox_allele", "ST",
                          "ALLELE1", "ALLELE2", "ALLELE3", "ALLELE4",
                          "ALLELE5", "ALLELE6", "ALLELE7",
                          "Identity", "Coverage", "db_version", "notes"])
        self.assertEqual(rows[0]["ALLELE1"], "atpA(1)")
        self.assertEqual(rows[0]["ALLELE2"], "dnaE(2)")
        self.assertEqual(rows[0]["ALLELE7"], "rpoB(1)")

    def test_separate_allele_columns_invalid_fasta(self):
        bad = os.path.join(self.tmp, "bad2.fasta")
        with open(bad, "w") as fh:
            fh.write("not fasta\n")
        r = self.run_tool(bad, extra_args=["--allele-layout", "separate"])[0]
        self.assertEqual(r["tox_gene"], "ERROR")
        self.assertEqual(r["ALLELE1"], "")

    def test_toxigenic_known_st(self):
        r = self.run_tool(self.genome("CD001", ["1"] * 7, tox=self.tox["2"]))[0]
        self.assertEqual((r["tox_gene"], r["tox_allele"], r["ST"]), ("tox", "2", "8"))
        self.assertEqual((r["Identity"], r["Coverage"]), ("100.0", "100.0"))

    def test_tox_absent(self):
        r = self.run_tool(self.genome("CD002", ["2"] * 7))[0]
        self.assertEqual((r["tox_gene"], r["tox_allele"], r["ST"]),
                         ("tox not present", "NA", "26"))
        self.assertEqual((r["Identity"], r["Coverage"]), ("0.0", "0.0"))

    def test_novel_tox_allele(self):
        novel = mutate_in_frame(self.tox["1"], [700])  # 1 SNP, 99.94%
        r = self.run_tool(self.genome("nov", ["1"] * 7, tox=novel))[0]
        self.assertEqual(r["tox_gene"], "novel tox allele")
        self.assertEqual(r["tox_allele"], "~1")
        self.assertGreaterEqual(float(r["Identity"]), 99.0)

    # SNP detail on novel alleles -------------------------------------------
    def test_novel_tox_reports_one_snp_and_no_disruption(self):
        novel = mutate_in_frame(self.tox["1"], [700])
        r = self.run_tool(self.genome("snp1", ["1"] * 7, tox=novel))[0]
        self.assertIn("closest tox allele 1 (99.9%): 1 SNP", r["notes"])
        self.assertNotIn("disrupted", r["notes"])

    def test_novel_tox_premature_stop_flagged(self):
        # Turn codon 100 into a stop codon: an in-frame nonsense change.
        seq = self.tox["1"]
        nonsense = seq[:297] + "TAA" + seq[300:]
        r = self.run_tool(self.genome("stop", ["1"] * 7, tox=nonsense))[0]
        self.assertEqual(r["tox_gene"], "novel tox allele")
        self.assertIn("premature stop at codon 100 of 561", r["notes"])
        self.assertIn("possibly disrupted gene", r["notes"])

    def test_novel_tox_frameshift_flagged(self):
        seq = self.tox["1"]
        shifted = seq[:411] + seq[412:]          # 1 bp deletion
        r = self.run_tool(self.genome("fs", ["1"] * 7, tox=shifted))[0]
        self.assertEqual(r["tox_gene"], "novel tox allele")
        self.assertIn("1 bp deletion", r["notes"])
        self.assertIn("frameshift, premature stop at codon", r["notes"])
        self.assertIn("possibly disrupted gene", r["notes"])

    def test_closest_allele_itself_disrupted_is_named(self):
        """A near match to tox_3, which is frameshifted in the database: the
        row must say the disruption belongs to the allele, not to this SNP."""
        novel = mutate(self.tox["3"], [900])
        r = self.run_tool(self.genome("inh", ["1"] * 7, tox=novel))[0]
        self.assertEqual(r["tox_allele"], "~3")
        self.assertIn("allele 3 is itself a disrupted tox allele", r["notes"])

    def test_novel_tox_notes_give_no_instruction(self):
        novel = mutate_in_frame(self.tox["1"], [700])
        r = self.run_tool(self.genome("noadvice", ["1"] * 7, tox=novel))[0]
        self.assertNotIn("confirm novel allele", r["notes"])
        self.assertNotIn("phenotypic", r["notes"])

    def test_novel_mlst_allele_reports_snp_count(self):
        prof = ["1"] * 7
        prof[3] = mutate(self.alleles["fusA"]["1"], [222])
        r = self.run_tool(self.genome("nalc", prof))[0]
        self.assertIn("fusA: 1 SNP from allele 1", r["notes"])

    def test_split_tox_flagged_for_review_not_absent(self):
        r = self.run_tool(self.genome("split", ["1"] * 7, tox=self.tox["1"],
                                      split_tox=True))[0]
        self.assertEqual(r["tox_gene"], "review")
        self.assertIn("contig edge", r["notes"])

    def test_divergent_tox_review(self):
        div = mutate(self.tox["1"], list(range(5, 1683, 40)))   # ~97.5%
        r = self.run_tool(self.genome("div", ["1"] * 7, tox=div))[0]
        self.assertEqual(r["tox_gene"], "review")
        self.assertIn("identity below threshold", r["notes"])

    def test_novel_st_combination(self):
        r = self.run_tool(self.genome("nst", ["1", "2", "1", "1", "1", "1", "1"]))[0]
        self.assertEqual(r["ST"], "novel")
        self.assertIn("dnaE(2)", r["mlst_profile"])

    def test_novel_mlst_allele(self):
        prof = ["1"] * 7
        prof[3] = mutate(self.alleles["fusA"]["1"], [222])
        r = self.run_tool(self.genome("nal", prof))[0]
        self.assertEqual(r["ST"], "NA")
        self.assertIn("fusA(~1)", r["mlst_profile"])

    def test_no_loci(self):
        r = self.run_tool(os.path.join(ROOT, "tool", "test-data", "no_loci.fasta"))[0]
        self.assertEqual((r["tox_gene"], r["ST"]), ("tox not present", "NA"))
        self.assertIn("no MLST loci found", r["notes"])

    def test_invalid_fasta_does_not_crash_batch(self):
        bad = os.path.join(self.tmp, "bad.fasta")
        with open(bad, "w") as fh:
            fh.write("this is not a fasta file\n")
        good = self.genome("good", ["1"] * 7)
        rows = self.run_tool(bad, good)
        self.assertEqual(rows[0]["tox_gene"], "ERROR")
        self.assertEqual(rows[1]["ST"], "8")

    def test_sample_name_extension_stripped(self):
        a = self.genome("ext", ["1"] * 7)
        rows = self.run_tool(a, extra_args=["-n", "CD009.contigs.fasta"])
        self.assertEqual(rows[0]["Sample"], "CD009")

    def test_multiple_samples_with_names(self):
        a = self.genome("a", ["1"] * 7, tox=self.tox["1"])
        b = self.genome("b", ["2"] * 7)
        rows = self.run_tool(a, b, extra_args=["-n", "S1", "-n", "S2"])
        self.assertEqual([r["Sample"] for r in rows], ["S1", "S2"])


class FakeHit:
    """Just the fields diff_summary / tox_novel_note read."""

    def __init__(self, qseq, sseq, allele="1", identity=99.9, coverage=100.0,
                 qstart=1, qend=None):
        self.qseq, self.sseq, self.allele = qseq, sseq, allele
        self.identity, self.coverage = identity, coverage
        self.qstart = qstart
        self.qend = qend if qend is not None else len(qseq.replace("-", ""))
        self.qlen = len(qseq.replace("-", ""))


class HelperTests(unittest.TestCase):
    """The difference wording and the gene-integrity check, on their own."""

    def test_diff_summary_counts(self):
        self.assertEqual(dt.diff_summary(FakeHit("ACGT", "ACGT")),
                         "no differences")
        self.assertEqual(dt.diff_summary(FakeHit("ACGT", "ACTT")), "1 SNP")
        self.assertEqual(dt.diff_summary(FakeHit("ACGTA", "ACTTG")), "2 SNPs")
        self.assertEqual(dt.diff_summary(FakeHit("ACGT", "AC-T")),
                         "1 bp deletion")
        self.assertEqual(dt.diff_summary(FakeHit("AC-GT", "ACAGT")),
                         "1 bp insertion")
        self.assertEqual(dt.diff_summary(FakeHit("ACGTA", "AC-TG")),
                         "1 SNP, 1 bp deletion")

    def test_first_stop(self):
        self.assertIsNone(dt.first_stop("GTGAAACCC"))
        self.assertEqual(dt.first_stop("GTGTAACCC"), 2)
        self.assertEqual(dt.first_stop("GTGAAATGA"), 3)
        # Out of frame, so not counted as a stop codon.
        self.assertIsNone(dt.first_stop("GTGATAACCC"))

    def test_cds_status(self):
        self.assertEqual(dt.cds_status("GTGAAATGA"), ("intact", 3, 3))
        self.assertEqual(dt.cds_status("GTGTGAAAATGA"), ("disrupted", 2, 4))
        # No stop codon at all.
        self.assertEqual(dt.cds_status("GTGAAACCC"), ("disrupted", None, 3))

    def test_integrity_judged_on_found_gene_not_on_closest_allele(self):
        """A change that restores the reading frame gives an intact gene even
        when the closest allele is itself frameshifted.

        This is why integrity is read off the found sequence: comparing indel
        lengths against the closest allele would call this one disrupted.
        """
        intact = rand_orf(60)                    # 180 bp, clean ORF
        broken = intact[:30] + intact[31:]       # the database allele: -1 bp
        self.assertEqual(dt.cds_status(broken)[0], "disrupted")
        self.assertEqual(dt.cds_status(intact)[0], "intact")

        # The assembly carries `intact`; the closest allele is `broken`, so the
        # alignment has a 1 bp insertion relative to it.
        hit = FakeHit(qseq=broken[:30] + "-" + broken[30:], sseq=intact,
                      allele="3")
        note = dt.tox_novel_note(hit, {"3": 11})
        self.assertIn("1 bp insertion", note)
        self.assertNotIn("possibly disrupted gene", note)
        self.assertIn("allele 3 is itself a disrupted tox allele", note)

    def test_partial_alignment_is_not_judged(self):
        seq = rand_orf(60)
        hit = FakeHit(qseq=seq, sseq=seq, qend=150, coverage=83.3)
        note = dt.tox_novel_note(hit, {})
        self.assertIn("gene not full length, integrity not checked", note)
        self.assertNotIn("possibly disrupted gene", note)


if __name__ == "__main__":
    unittest.main()
