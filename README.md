# Galaxy *Corynebacterium diphtheriae* Typing Tool

Rapid typing of *C. diphtheriae* genome assemblies for Galaxy Australia:

- **tox** gene detection and **tox allele** assignment
- **7-gene MLST** sequence type (ST)
- match quality (identity / coverage) and the **database version** on every row
- one tabular report per run, ready for Excel and downstream reporting

Database: [BIGSdb-Pasteur *C. diphtheriae* species complex](https://bigsdb.pasteur.fr/diphtheria/),
refreshed twice yearly through a reviewed pull request.

> **tox gene ≠ toxin production.** Non-toxigenic *tox*-bearing (NTTB) strains
> exist. Confirm toxigenicity phenotypically (e.g. modified Elek test) as per
> laboratory procedure. This tool does not determine species — the Species
> column is taken from the laboratory identification.

---

## Output

```tsv
Sample  Species                      tox_gene         tox_allele  ST   Identity  Coverage  db_version       mlst_profile                                   notes
CD001   Corynebacterium diphtheriae  tox              23          534  100.0     100.0     Pasteur_2026-10  atpA(..);dnaE(..);...
CD002   Corynebacterium diphtheriae  tox not present  NA          534  0.0       0.0       Pasteur_2026-10  atpA(..);dnaE(..);...
```

The first eight columns are the agreed report format; `mlst_profile` and
`notes` are appended for traceability and QC.

### tox calls

| `tox_gene` | Rule | `tox_allele` |
|---|---|---|
| `tox` | 100% identity over 100% of a known allele | allele number |
| `novel tox allele` | identity ≥ 99% and coverage ≥ 95%, not exact | `~N` (closest N) |
| `review` | tox-like hit (≥ 30% coverage) failing thresholds — e.g. gene split across contigs, or divergent | `?N` |
| `tox not present` | no tox-like hit ≥ 30% coverage | `NA` |

`review` exists so a gene broken by an assembly gap is **never silently
reported as absent**. Thresholds are adjustable in Galaxy.

### MLST calls

Per-locus notation (as in `tseemann/mlst`): `N` exact, `~N` novel allele
(closest N), `?N` partial, `-` missing. `ST` is the sequence type, `novel`
(all alleles known, combination new) or `NA` (any locus not exact).

---

## Repository layout

```text
tool/        Galaxy wrapper (tool.xml) + diphtheria_typing.py + test-data
database/    VERSION.txt, tox/, mlst/  (populated by the updater, kept out of code)
scripts/     update_database.py, bundle_tool.sh
tests/       unit tests (synthetic database - logic only)
workflows/   Galaxy workflow (.ga)
docs/        validation record
.github/     tests, twice-yearly database update, release
```

## Getting started

```bash
# 1. BLAST+ and Python 3.8+
sudo apt install ncbi-blast+          # or: conda install -c bioconda blast

# 2. Populate the database (see Authentication below)
pip install requests requests-oauthlib
python scripts/update_database.py --db database

# 3. Run
python tool/diphtheria_typing.py -d database -i CD001.fasta -i CD002.fasta -o results.tsv

# 4. Tests
python -m unittest discover -s tests -v
```

## Database updates

`.github/workflows/update_database.yml` runs at 00:00 UTC on **1 April and
1 October** (or manually from the Actions tab). It downloads the current
alleles and profiles, compares them with `database/`, runs the tests and opens
a **pull request** with a change table and reviewer checklist. Nothing is
deployed automatically.

After review and merge, tag a release (`git tag v1.1 && git push --tags`);
`release.yml` builds a Galaxy-ready tarball with the database bundled.

| Release | Database |
|---|---|
| v1.0 | Pasteur_2026-10 |
| v1.1 | Pasteur_2027-04 |
| v1.2 | Pasteur_2027-10 |

### Authentication (important)

Since 1 January 2025 BIGSdb-Pasteur only serves data curated after that date to
authenticated API users. Without credentials the updater still works but
**newer alleles and STs will be missing**, and the PR says so in bold.

Register at bigsdb.pasteur.fr, obtain OAuth1 API tokens, and add these
repository secrets (Settings → Secrets and variables → Actions):

`BIGSDB_CONSUMER_KEY`, `BIGSDB_CONSUMER_SECRET`, `BIGSDB_TOKEN`, `BIGSDB_TOKEN_SECRET`

Also enable *Settings → Actions → General → Allow GitHub Actions to create and
approve pull requests*.

## Galaxy installation

```bash
bash scripts/bundle_tool.sh            # -> build/cdiphtheriae_typing/
```

Install the bundled directory as a local tool, or use the release tarball.
Inputs accept one or more FASTA datasets or a collection; sample names come
from the dataset / collection element names.

## Roadmap

- **1.1** improved reporting, Galaxy visualisations
- **1.2** cgMLST, novel allele reporting, additional QC metrics
- **2.0** spuA, pbp and further virulence / resistance markers

## Citation

- Bolt F. *et al.* (2010) Multilocus sequence typing identifies evidence for
  recombination and two distinct lineages of *Corynebacterium diphtheriae*.
  *J Clin Microbiol* 48:4177–4185.
- BIGSdb-Pasteur, Institut Pasteur — https://bigsdb.pasteur.fr/diphtheria/
- Camacho C. *et al.* (2009) BLAST+: architecture and applications.
  *BMC Bioinformatics* 10:421.
