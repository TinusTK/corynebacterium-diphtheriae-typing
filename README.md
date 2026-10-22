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
> laboratory procedure.

---

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

| Column | Meaning |
|---|---|
| `Sample` | file name, or dataset / collection element name in Galaxy |
| `tox_gene` | `tox`, `novel tox allele`, `review` or `tox not present` |
| `tox_allele` | allele number; `~N` novel (closest N); `?N` partial; `NA` |
| `ST` | sequence type; `novel` (known alleles, new combination); `NA` |
| `mlst_profile` | per-locus allele calls (see MLST calls below) |
| `Identity` | % identity of the best tox hit |
| `Coverage` | % of the tox allele covered by the best hit |
| `db_version` | database used, e.g. `Pasteur_2026-10` |
| `notes` | QC flags and reasons for review; blank when nothing to flag |

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

## Running without Galaxy

The typing script is a normal command-line program; the Galaxy wrapper is just
a thin layer on top. Both use the same code and the same database.

```bash
# 1. Requirements: BLAST+ and Python 3.8+ (no other Python packages)
sudo apt install ncbi-blast+          # or: conda install -c bioconda blast

# 2. Get the code and a database
git clone https://github.com/TinusTK/galaxy-corynebacterium-diphtheriae-typing
cd galaxy-corynebacterium-diphtheriae-typing
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
Options: `cdiph-typing --help` (thresholds, threads, sample names).

Tests: `python -m unittest discover -s tests -v`

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

#### Pasteur API credentials (one-time setup)

1. **Create an account** at https://bigsdb.pasteur.fr/register/ and register
   for the *C. diphtheriae* species complex **sequence definitions** database
   (`pubmlst_diphtheria_seqdef`).
2. **Request an API key** by emailing bigsdb@pasteur.fr (academic / public
   health) with your username, affiliation, institutional email, sector,
   the database name, and how you'll use it (automated twice-yearly typing
   database updates). See https://bigsdb.pasteur.fr/requesting-api-key.
   They reply with a **consumer key** and **consumer secret**.
3. **Get an access token** on your own computer (once):
   ```bash
   pip install requests requests-oauthlib
   python scripts/bigsdb_auth.py
   ```
   Paste the consumer key/secret, open the printed link, log in, click
   Authorize, paste the code back. It prints the access token and secret.
4. **Add four repository secrets** (Settings → Secrets and variables → Actions
   → New repository secret):

   | Name | Value |
   |---|---|
   | `BIGSDB_CONSUMER_KEY` | consumer key from Pasteur |
   | `BIGSDB_CONSUMER_SECRET` | consumer secret from Pasteur |
   | `BIGSDB_ACCESS_TOKEN` | from `bigsdb_auth.py` |
   | `BIGSDB_ACCESS_SECRET` | from `bigsdb_auth.py` |

The access token doesn't expire unless revoked; the updater swaps it for a
12-hour session token on each run. If an update fails with a session-token
error, re-run step 3 and replace the two access secrets.

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
