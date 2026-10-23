# Validation record

Reference genomes in `validation/expected.tsv` are re-checked automatically
by the *Validate* workflow on every database update.

The unit tests use synthetic sequences and only verify calling logic.
Before clinical or public-health use, validate against characterised genomes
and record results here for each database version.

Results for database **Pasteur_2026-10** (tool v1.0, 2026-10-23). Full output:
[`validation/results/validation_Pasteur_2026-10.tabular`](../validation/results/validation_Pasteur_2026-10.tabular).

| Genome | Accession | Expected tox | Result tox (allele) | Expected ST | Result ST | Pass |
|---|---|---|---|---|---|---|
| NCTC 13129 | GCF_000195815.1 | positive | tox (1) | 8 | 8 | ✅ |
| C7(β) | CP003210 | positive | tox (2) | 26 | 26 | ✅ |
| PW8 | CP003216 | positive | tox (2) | – | 44 | ✅ tox |
| CDC-E8392 | CP003211 | positive | tox (4) | – | 50 | ✅ tox |
| 31A | CP003206 | positive | tox (3) | – | 584 | ✅ tox |
| 241 | CP003207 | negative | not present | – | 175 | ✅ tox |
| VA01 | CP003217 | negative | not present | – | 80 | ✅ tox |
| HC01 | CP003212 | negative | not present | – | 175 | ✅ tox |
| *C. glutamicum* ATCC 13032 (negative control) | GCF_000011325.1 | negative | not present | none | NA, no MLST loci | ✅ |

All tox hits were exact (100% identity, 100% coverage), with no QC flags.

**Summary:** tox status 9/9 correct (5 positive, 3 negative, 1 negative control);
ST 2/2 correct where an independently published ST exists. "–" means no
independently published ST was found, so the result is recorded but not
counted as verified. Strains 241 and HC01 both typed as ST175 with identical
profiles, consistent with Trost et al. describing them as the most closely
related pair in their study.

**Sources for expected values**
- NCTC 13129 ST8, tox+: Guglielmini *et al.* 2021, *J Clin Microbiol*, doi:10.1128/JCM.01581-21
- tox status of C7(β), PW8, CDC-E8392, 31A, 241, VA01, HC01: Trost *et al.* 2012,
  *J Bacteriol* 194:3199, doi:10.1128/JB.00183-12
- C7(β) ST26: Choudhury *et al.*, *J Clin Diagn Res* (MLST of diphtheria isolates, Assam)

Suggested set: reference strains, local tox-positive and tox-negative
isolates, an NTTB isolate if available, an assembly with tox split across
contigs, and a non-*diphtheriae* *Corynebacterium* negative control.
