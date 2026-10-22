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
LOCI = ["atpA", "dnaE", "dnaK", "fusA", "leuA", "odhA", "rpoB"]
rng = random.Random(42)


def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, positions):
    s = list(seq)
    for p in positions:
        s[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[p]]
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
        tox = rand_seq(1683)
        cls.tox = {"1": tox, "2": mutate(tox, [100, 900, 1500])}
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
        self.assertEqual(list(rows[0].keys())[:8],
                         ["Sample", "Species", "tox_gene", "tox_allele", "ST",
                          "Identity", "Coverage", "db_version"])
        self.assertEqual(rows[0]["db_version"], "SYNTHETIC_TEST")
        self.assertEqual(rows[0]["Sample"], "cols")

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
        novel = mutate(self.tox["1"], [700])           # 1 SNP, 99.94%
        r = self.run_tool(self.genome("nov", ["1"] * 7, tox=novel))[0]
        self.assertEqual(r["tox_gene"], "novel tox allele")
        self.assertEqual(r["tox_allele"], "~1")
        self.assertGreaterEqual(float(r["Identity"]), 99.0)

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

    def test_multiple_samples_with_names(self):
        a = self.genome("a", ["1"] * 7, tox=self.tox["1"])
        b = self.genome("b", ["2"] * 7)
        rows = self.run_tool(a, b, extra_args=["-n", "S1", "-n", "S2"])
        self.assertEqual([r["Sample"] for r in rows], ["S1", "S2"])


if __name__ == "__main__":
    unittest.main()
