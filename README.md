# *Corynebacterium diphtheriae* Typing Tool

I built this tool for public health use, to provide rapid typing of
*Corynebacterium diphtheriae* genome assemblies:

- **tox** gene detection and **tox allele** assignment
- **7-gene MLST** sequence type (ST)
- match quality (identity / coverage) and the **database version** on every row
- one tabular report per run, ready for Excel and downstream reporting

> **Prefer no command line?**
> [Genomics Workflow Microbials (GWM)](https://github.com/TinusTK/GenomicsWorkflowMicrobials)
> is a free desktop app for Windows, Mac and Linux with a point-and-click
> interface. *C. diphtheriae* typing is built in: choose a folder of reads or
> assemblies, press Run and read the results as tables. Nothing to install
> by hand, and no server needed: it all runs locally on your computer.

Database: [BIGSdb-Pasteur *C. diphtheriae* species complex](https://bigsdb.pasteur.fr/diphtheria/),
refreshed twice yearly through a reviewed pull request.

> **tox gene ≠ toxin production.** Non-toxigenic *tox*-bearing (NTTB) strains
> exist. Confirm toxigenicity phenotypically (e.g. modified Elek test) as per
> laboratory procedure.

---

Available on Galaxy Australia.

## Output

One tab-separated row per assembly. Example run on five assemblies
(values are illustrative, chosen to show each kind of result):

```tsv
Sample	tox_gene	tox_allele	ST	mlst_profile	Identity	Coverage	db_version	notes
CD001	tox	23	534	atpA(1);dnaE(2);dnaK(1);fusA(1);leuA(10);odhA(4);rpoB(3)	100.0	100.0	Pasteur_2026-10	
CD002	tox not present	NA	534	atpA(1);dnaE(2);dnaK(1);fusA(1);leuA(10);odhA(4);rpoB(3)	0.0	0.0	Pasteur_2026-10	
CD003	novel tox allele	~23	8	atpA(1);dnaE(1);dnaK(1);fusA(1);leuA(1);odhA(1);rpoB(1)	99.9	100.0	Pasteur_2026-10	closest tox allele 23; confirm novel allele
CD004	review	?23	novel	atpA(1);dnaE(2);dnaK(1);fusA(7);leuA(10);odhA(4);rpoB(3)	100.0	53.5	Pasteur_2026-10	exact alleles at all loci but profile not in database; tox-like hit: partial coverage; hit at contig edge (possible assembly break)
CD005	tox not present	NA	NA	atpA(1);dnaE(2);dnaK(1);fusA(~4);leuA(10);odhA(4);rpoB(3)	0.0	0.0	Pasteur_2026-10	
```

The same report opened in Excel or Galaxy:

| Sample | tox_gene | tox_allele | ST | mlst_profile | Identity | Coverage | db_version | notes |
|---|---|---|---|---|---:|---:|---|---|
| CD001 | tox | 23 | 534 | atpA(1);dnaE(2);dnaK(1);fusA(1);leuA(10);odhA(4);rpoB(3) | 100.0 | 100.0 | Pasteur_2026-10 | |
| CD002 | tox not present | NA | 534 | atpA(1);dnaE(2);dnaK(1);fusA(1);leuA(10);odhA(4);rpoB(3) | 0.0 | 0.0 | Pasteur_2026-10 | |
| CD003 | novel tox allele | ~23 | 8 | atpA(1);dnaE(1);dnaK(1);fusA(1);leuA(1);odhA(1);rpoB(1) | 99.9 | 100.0 | Pasteur_2026-10 | closest tox allele 23; confirm novel allele |
| CD004 | review | ?23 | novel | atpA(1);dnaE(2);dnaK(1);fusA(7);leuA(10);odhA(4);rpoB(3) | 100.0 | 53.5 | Pasteur_2026-10 | exact alleles at all loci but profile not in database; tox-like hit: partial coverage; hit at contig edge (possible assembly break) |
| CD005 | tox not present | NA | NA | atpA(1);dnaE(2);dnaK(1);fusA(~4);leuA(10);odhA(4);rpoB(3) | 0.0 | 0.0 | Pasteur_2026-10 | |

How to read it:

- **CD001** — tox detected, exact match to allele 23, ST534.
- **CD002** — no tox, ST534.
- **CD003** — tox with a single SNP difference from allele 23: a candidate new allele to confirm.
- **CD004** — tox found on only 53.5% of its length at the end of a contig, so it is flagged for review rather than called absent. All seven MLST alleles are known but the combination is new (`novel` ST).
- **CD005** — no tox; fusA is a new allele (`~4`, closest 4), so no ST can be assigned.

### Columns

| Column | What it is | Meaning |
|---|---|---|
| `Sample` | The assembly's name: the input file name with the FASTA extension stripped (`.fasta`, `.fa`, `.fna`, `.fas`, `.fsa`, `.contigs`, `.scaffolds`, and `.gz`), or in Galaxy the dataset or collection element name. | Identifies which genome the row belongs to. Worth a glance when you run a collection, since two files whose names differ only by extension end up with the same label. |
| `tox_gene` | One of `tox` (exact match to a known allele), `novel tox allele`, `review`, `tox not present`, or `ERROR` if the file could not be read as FASTA. Rules in *tox calls* below. | The headline result: whether the diphtheria toxin gene is present in the assembly. This is a prediction from sequence, not a test of toxin production, so a `tox` or `novel tox allele` call still needs laboratory confirmation before an isolate is called toxigenic. `review` means look at the assembly before accepting any answer. |
| `tox_allele` | The allele number of the matched tox sequence: a plain number for an exact match, `N,M` when two or more different alleles match exactly, `~N` for a near match (closest allele is N), `?N` for a partial match to N, `NA` when nothing tox-like was found. | Ties the hit to a named allele in the Pasteur database, so the result can be compared with other labs and with earlier isolates. `~N` is a candidate new allele worth submitting; `N,M` or `?N` more often points at the assembly — a duplicated or broken gene — than at the biology. |
| `ST` | The 7-gene MLST sequence type, looked up in the database's profile table: a number when all seven alleles are exact and the combination is listed, `novel` when all seven are exact but the combination is not listed, `NA` when any locus is not a single exact allele. | Places the isolate in the global MLST scheme, which is what allows comparison with published strains. It is a table lookup, so it is only as complete as the bundled database — see `db_version`. `novel` is a candidate new ST; `NA` usually means an incomplete assembly or not *C. diphtheriae*. ST is too coarse to settle transmission questions on its own. |
| `mlst_profile` | The call at each of the seven loci in database order, as `locus(call)` joined by `;` — e.g. `atpA(1);dnaE(2);...`. Notation follows `tseemann/mlst`: `N` exact, `~N` novel (closest N), `?N` partial, `-` not found or below the identity threshold, `N,M` two exact copies. See *MLST calls* below. | Shows which locus caused a `novel` or `NA` ST, so you can tell a near miss from a failed assembly. One `~N` among six exact calls is a candidate new allele; several `-` suggests the wrong species or poor coverage. |
| `ALLELE1`–`ALLELE7` | Only with the separate allele layout, in place of `mlst_profile`: the same seven calls, one per column, each `locus(call)` such as `atpA(1)`, same notation. Blank on an `ERROR` row. | The same information as `mlst_profile`, laid out so one locus can be sorted or filtered on its own in Excel or read by a downstream tool. |
| `Identity` | % nucleotide identity of the best tox hit, to one decimal: `100.0` for an exact match, `0.0` when no tox-like hit was found, `NA` if the file could not be read. Refers to tox only, not to the MLST loci. | With `Coverage`, says how good the tox match is and why a call fell short of `tox`. High identity with low coverage points to a gene split by an assembly gap; low identity to a divergent or related sequence. |
| `Coverage` | % of the matched tox allele's length covered by the best hit, to one decimal: `100.0` exact, `0.0` no hit, `NA` unreadable file. Also tox only. | The usual reason a tox-positive genome is flagged `review` rather than called `tox`: a partial gene. Read it alongside `notes`, which says whether the hit sits at a contig edge. |
| `db_version` | The label of the bundled Pasteur database, read from `database/VERSION.txt` — e.g. `Pasteur_2026-10b`; `unknown` if that file is missing. | Every allele number and ST in the row comes from this version, so keep it with the result. A `novel` call may just be an allele or ST added to the database after this snapshot, so re-run on a newer release before reporting it as new. |
| `notes` | QC flags and the reason for a `review` or `novel` call, joined by `; `. Blank when there is nothing to flag. For example `assembly only 812,345 bp - check completeness`, `closest tox allele 23; confirm novel allele`, `tox-like hit: partial coverage; hit at contig edge (possible assembly break)`, `exact alleles at all loci but profile not in database`, `no MLST loci found - check species/assembly`, `atpA: multiple exact alleles`, `weak tox-like fragment only (12.0% cov, 88.0% id)`. | Read this first on any row that is not a clean call: it says what to check. A `weak tox-like fragment` note is the one case where an apparently tox-negative row still deserves a look. |

To get one column per MLST locus instead of `mlst_profile`, choose **MLST
allele layout → One column per allele** in Galaxy, or add
`--allele-layout separate` on the command line. `mlst_profile` is then
replaced by `ALLELE1` to `ALLELE7`, each holding one call, e.g. `atpA(1)`.

### tox calls

| `tox_gene` | Rule | `tox_allele` |
|---|---|---|
| `tox` | 100% identity over 100% of a known allele | allele number |
| `novel tox allele` | identity ≥ 99% and coverage ≥ 95%, not exact | `~N` (closest N) |
| `review` | tox-like hit (≥ 30% coverage) failing thresholds — e.g. gene split across contigs, or divergent | `?N` |
| `tox not present` | no tox-like hit ≥ 30% coverage (a weaker fragment is still recorded in `notes`) | `NA` |
| `ERROR` | the file could not be read as FASTA | `NA` (and `ST`, `Identity`, `Coverage` all `NA`) |

`review` exists so a gene broken by an assembly gap is **never silently
reported as absent**. Thresholds are adjustable in Galaxy.

### MLST calls

Per-locus notation (as in `tseemann/mlst`): `N` exact, `~N` novel allele
(closest N), `?N` partial, `-` missing. `ST` is the sequence type, `novel`
(all alleles known, combination new) or `NA` (any locus not exact).

---

## Software and data used

| Component | Used for | Version |
|---|---|---|
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) (`makeblastdb`, `blastn`) | Builds a BLAST database from each assembly and searches every tox and MLST allele against it; identity and coverage come from these alignments | 2.16.0 in Galaxy (any ≥ 2.9 standalone) |
| Python 3 | Runs `diphtheria_typing.py`: allele calling, ST lookup, tox rules, report. Standard library only, no extra packages | 3.12 in Galaxy (≥ 3.8 standalone) |
| [BIGSdb-Pasteur](https://bigsdb.pasteur.fr/diphtheria/) *C. diphtheriae* species complex database | tox alleles, MLST alleles and ST profiles, bundled with each release | see `database/VERSION.txt` |

**How it works:** each assembly is turned into a BLAST database, and all tox and
MLST alleles are searched against it with `blastn` (e-value ≤ 1e-20, ≥ 80%
identity, no low-complexity filtering). An exact call needs 100% identity across
the full allele length with no gaps. Otherwise the best hit is reported with its
identity and coverage, and the calling rules above apply.

## Validation

Validated against nine public genomes: tox status correct for 9/9 (five
tox-positive, three tox-negative, one non-*diphtheriae* negative control), with
exact allele matches, and ST correct for both strains with an independently
published ST (NCTC 13129 ST8, C7(β) ST26). Validation re-runs automatically on
every database update. Results, method, limitations and references:
[validation/VALIDATION.md](validation/VALIDATION.md).

## Running without Galaxy

The typing script is a normal command-line program; the Galaxy wrapper is just
a thin layer on top. Both use the same code and the same database.

```bash
# 1. Requirements: BLAST+ and Python 3.8+ (no other Python packages)
sudo apt install ncbi-blast+          # or: conda install -c bioconda blast

# 2. Get the code and a database
git clone https://github.com/TinusTK/corynebacterium-diphtheriae-typing
cd corynebacterium-diphtheriae-typing
#    either download a release tarball (database included), or build it:
pip install requests requests-oauthlib
python scripts/update_database.py --db database   # see Authentication below

# 3a. Run directly
python tool/diphtheria_typing.py -d database -i CD001.fasta -i CD002.fasta -o results.tsv

# 3b. ...or install the `cdiph-typing` command
pip install .
cdiph-typing -d database -i assemblies/*.fasta -o results.tsv
```

Batch a folder: `cdiph-typing -d database $(for f in assemblies/*.fasta; do echo -i $f; done) -o results.tsv`.
Options: `cdiph-typing --help` (thresholds, threads, sample names, allele layout).

Tests: `python -m unittest discover -s tests -v`

## Database updates

`.github/workflows/update_database.yml` runs at 00:00 UTC on **1 April and
1 October** (or manually from the Actions tab). It downloads the current
alleles and profiles, compares them with `database/`, runs the tests and opens
a **pull request** with a change table and reviewer checklist. Nothing is
deployed automatically.

**Before each release, bump the tool version** in `tool/tool.xml`. The Tool
Shed rejects a version that's already published:

- database update or wrapper change only → raise the suffix:
  `1.0.0+galaxy0` → `1.0.0+galaxy1`
- change to `diphtheria_typing.py` → raise `@TOOL_VERSION@` (also in
  `pyproject.toml` and `__version__`) and reset the suffix: `1.1.0+galaxy0`

CI shows a warning on every push until the version is bumped. After review and
merge, tag a release named after the tool version with a `v` in front
(`git tag v1.1.0+galaxy0 && git push --tags`, or *Releases → Draft a new
release* on GitHub). The GitHub release is titled with the Tool Shed version
and its notes give the database version. `release.yml` then bundles the
tool with the database, runs the Galaxy tool tests, creates a GitHub release,
and publishes the new revision to the Galaxy Tool Shed.

| Release | Database |
|---|---|
| v1.0.0+galaxy0 | Pasteur_2026-10 |
| v1.0.0+galaxy1 | Pasteur_2026-10 |
| v1.0.0+galaxy2 | Pasteur_2026-10b (first authenticated download) |
| v1.1.0+galaxy0 | Pasteur_2026-10b (adds the one-column-per-allele option) |

### Authentication (important)

Since 1 January 2025 BIGSdb-Pasteur only serves data curated after that date to
authenticated API users. Without credentials the updater still works but
**newer alleles and STs will be missing**, and the PR says so in bold.

## Citation

- Bolt F. *et al.* (2010) Multilocus sequence typing identifies evidence for
  recombination and two distinct lineages of *Corynebacterium diphtheriae*.
  *J Clin Microbiol* 48:4177–4185.
- BIGSdb-Pasteur, Institut Pasteur — https://bigsdb.pasteur.fr/diphtheria/
- Camacho C. *et al.* (2009) BLAST+: architecture and applications.
  *BMC Bioinformatics* 10:421.
- Jolley K.A. & Maiden M.C.J. (2010) BIGSdb: Scalable analysis of bacterial
  genome variation at the population level. *BMC Bioinformatics* 11:595.
- TinusTK. *Corynebacterium diphtheriae* Typing Tool. GitHub —
  https://github.com/TinusTK/corynebacterium-diphtheriae-typing
