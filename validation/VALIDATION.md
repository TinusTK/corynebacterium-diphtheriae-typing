# Validation record

The typing tool was validated against public, fully sequenced genomes whose tox
status, and where available ST, have been published independently [1–5]. The
check is automated: the *Validate* workflow (`.github/workflows/validate.yml`)
downloads every genome in [`expected.tsv`](expected.tsv) from NCBI [8, 9],
types it with the current database, and fails if any result differs from the
expected value. It runs whenever the database or the typing script changes, so
every database update is re-validated before release.

The unit tests (`tests/`) use synthetic sequences and only check the calling
logic. This record is the real-genome evidence.

## Files in this folder

| File | Contents |
|---|---|
| `VALIDATION.md` | This record: method, results, references |
| [`expected.tsv`](expected.tsv) | Genomes, accessions, expected ST and tox status, with reference numbers |
| [`validation_Pasteur_2026-10.tabular`](validation_Pasteur_2026-10.tabular) | Full tool output for all genomes below (Galaxy tabular format) |

## Results: database Pasteur_2026-10, tool v1.0 (2026-10-23)

| Genome | Accession | Expected tox | Result tox (allele) | Expected ST | Result ST | Pass |
|---|---|---|---|---|---|---|
| NCTC 13129 | GCF_000195815.1 [3] | positive [1] | tox (1) | 8 [1] | 8 | ✅ |
| C7(β) | CP003210 [2] | positive [2] | tox (2) | 26 [4] | 26 | ✅ |
| PW8 | CP003216 [2] | positive [2] | tox (2) | – | 44 | ✅ tox |
| CDC-E8392 | CP003211 [2] | positive [2] | tox (4) | – | 50 | ✅ tox |
| 31A | CP003206 [2] | positive [2] | tox (3) | – | 584 | ✅ tox |
| 241 | CP003207 [2] | negative [2] | not present | – | 175 | ✅ tox |
| VA01 | CP003217 [2] | negative [2] | not present | – | 80 | ✅ tox |
| HC01 | CP003212 [2] | negative [2] | not present | – | 175 | ✅ tox |
| *C. glutamicum* ATCC 13032 (negative control) | GCF_000011325.1 [5] | negative | not present | none | NA, no MLST loci | ✅ |

All tox hits were exact (100% identity, 100% coverage), with no QC flags.

### Summary

- **tox status:** 9/9 correct (5 tox-positive, 3 tox-negative, 1 non-*diphtheriae* negative control).
- **ST:** 2/2 correct where an independently published ST exists (NCTC 13129 ST8 [1], C7(β) ST26 [4]).
  "–" means no independently published ST was found; those STs are recorded but not counted as verified.
- **Consistency check:** strains 241 and HC01 both typed as ST175 with identical
  allele profiles, consistent with Trost *et al.* [2] describing these two tox-negative
  strains as the most closely related pair in their study.
- **Specificity:** the *C. glutamicum* negative control returned no tox hit and no MLST
  loci, so the tool does not report false types for a related non-*diphtheriae* species.

### MLST allele profiles (from the output file)

| Genome | atpA | dnaE | dnaK | fusA | leuA | odhA | rpoB | ST |
|---|---|---|---|---|---|---|---|---|
| NCTC 13129 | 3 | 5 | 6 | 5 | 3 | 3 | 6 | 8 |
| C7(β) | 8 | 2 | 16 | 1 | 3 | 3 | 12 | 26 |
| PW8 | 14 | 2 | 23 | 4 | 2 | 14 | 2 | 44 |
| CDC-E8392 | 2 | 2 | 4 | 1 | 3 | 3 | 2 | 50 |
| 31A | 2 | 2 | 124 | 3 | 3 | 2 | 3 | 584 |
| 241 | 2 | 1 | 37 | 19 | 24 | 3 | 4 | 175 |
| VA01 | 2 | 1 | 32 | 19 | 13 | 3 | 6 | 80 |
| HC01 | 2 | 1 | 37 | 19 | 24 | 3 | 4 | 175 |

## Method

1. Genomes were downloaded from NCBI: RefSeq assemblies through the NCBI Datasets
   API [8], GenBank complete chromosomes through NCBI E-utilities [9].
2. Each genome was typed with `diphtheria_typing.py` (v1.0.0) against the bundled
   BIGSdb-Pasteur database [6, 7], version `Pasteur_2026-10`, using BLAST+ [10]
   with the default thresholds (tox: ≥ 99% identity, ≥ 95% coverage; MLST exact match).
3. Results were compared automatically with `expected.tsv`. tox status counts as
   correct if a `tox` or `novel tox allele` call matches an expected positive, and
   `tox not present` matches an expected negative.

## Limitations

- Only 2 of the 8 *C. diphtheriae* genomes have an independently published ST, so
  ST accuracy is verified on a small set. tox detection is verified on all 9.
- All genomes are complete, closed assemblies. Draft assemblies from routine
  sequencing are more fragmented; the tool flags a tox gene split across contigs
  as `review` instead of calling it absent, but this case is covered only by the
  synthetic unit tests, not by a real genome.
- The database was built without BIGSdb-Pasteur API authentication, so it holds
  definitions curated up to the end of 2024. Newer alleles or STs are reported as
  novel or untyped.
- tox gene detection does not prove toxin production: non-toxigenic tox-bearing
  (NTTB) strains exist, so toxigenicity needs phenotypic confirmation (e.g. the
  modified Elek test).

## Adding genomes

Append a row to `expected.tsv`:
`sample  accession  expected_ST (or - if unpublished)  positive|negative  source`.
Use only publicly available genomes, never patient isolates, and give a
reference for every expected value.

## References

1. Guglielmini J, Hennart M, Badell E, *et al.*, Brisse S. Genomic epidemiology and
   strain taxonomy of *Corynebacterium diphtheriae*. *J Clin Microbiol*.
   2021;59(12):e01581-21. doi:10.1128/JCM.01581-21
2. Trost E, Blom J, Soares SC, *et al.* Pangenomic study of *Corynebacterium
   diphtheriae* that provides insights into the genomic diversity of pathogenic
   isolates from cases of classical diphtheria, endocarditis, and pneumonia.
   *J Bacteriol*. 2012;194(12):3199–3215. doi:10.1128/JB.00183-12
3. Cerdeño-Tárraga AM, Efstratiou A, Dover LG, *et al.* The complete genome sequence
   and analysis of *Corynebacterium diphtheriae* NCTC13129. *Nucleic Acids Res*.
   2003;31(22):6516–6523. doi:10.1093/nar/gkg874
4. Choudhury G, Gogoi N, Nath R, Sargiary P, Das PP, Bhuyan B, Borkotoki U.
   Molecular characterisation of *Corynebacterium diphtheriae* isolates of faucial
   diphtheria cases from Assam: a cross-sectional study. *J Clin Diagn Res*.
   2022;16(5):DC24–DC26. doi:10.7860/JCDR/2022/52934.16379
5. Kalinowski J, Bathe B, Bartels D, *et al.* The complete *Corynebacterium
   glutamicum* ATCC 13032 genome sequence and its impact on the production of
   L-aspartate-derived amino acids and vitamins. *J Biotechnol*.
   2003;104(1–3):5–25. doi:10.1016/S0168-1656(03)00154-8
6. BIGSdb-Pasteur, *Corynebacterium diphtheriae* species complex database.
   Institut Pasteur. https://bigsdb.pasteur.fr/diphtheria/
7. Jolley KA, Maiden MCJ. BIGSdb: Scalable analysis of bacterial genome variation
   at the population level. *BMC Bioinformatics*. 2010;11:595.
   doi:10.1186/1471-2105-11-595
8. NCBI Datasets. National Center for Biotechnology Information.
   https://www.ncbi.nlm.nih.gov/datasets/
9. NCBI Entrez Programming Utilities (E-utilities). National Center for
   Biotechnology Information. https://www.ncbi.nlm.nih.gov/books/NBK25501/
10. Camacho C, Coulouris G, Avagyan V, *et al.* BLAST+: architecture and
    applications. *BMC Bioinformatics*. 2009;10:421. doi:10.1186/1471-2105-10-421
