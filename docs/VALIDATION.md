# Validation record

Reference genomes in `validation/expected.tsv` are re-checked automatically
by the *Validate* workflow on every database update.

The unit tests use synthetic sequences and only verify calling logic.
Before clinical or public-health use, validate against characterised genomes
and record results here for each database version.

| Sample / accession | Expected ST | Expected tox | Result ST | Result tox | DB version | Pass |
|---|---|---|---|---|---|---|
| NCTC 13129 (GCF_000195815.1) | 8 | positive | 8 — atpA(3);dnaE(5);dnaK(6);fusA(5);leuA(3);odhA(3);rpoB(6) | tox, allele 1 (100.0% id, 100.0% cov) | Pasteur_2026-10 | ✅ 2026-10-23 |

Suggested set: reference strains, local tox-positive and tox-negative
isolates, an NTTB isolate if available, an assembly with tox split across
contigs, and a non-*diphtheriae* *Corynebacterium* negative control.
