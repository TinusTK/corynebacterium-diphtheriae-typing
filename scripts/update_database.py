#!/usr/bin/env python3
"""
update_database.py - refresh database/ from BIGSdb-Pasteur (C. diphtheriae
species complex seqdef database) and write a change summary.

    python scripts/update_database.py --db database --summary CHANGES_DB.md

Downloads:
  - 7-gene MLST scheme loci, allele FASTAs and profiles
  - tox locus alleles

Authentication
--------------
Since 1 January 2025 BIGSdb-Pasteur only returns data curated after that date
to authenticated API users. Without credentials this script still runs but
newer alleles/STs will be MISSING, and the summary says so in bold.

To authenticate (see README "Pasteur API credentials"):
  1. get a consumer key/secret from the BIGSdb-Pasteur team
  2. run scripts/bigsdb_auth.py once to get a long-lived access token
  3. set BIGSDB_CONSUMER_KEY, BIGSDB_CONSUMER_SECRET,
         BIGSDB_ACCESS_TOKEN, BIGSDB_ACCESS_SECRET
This script exchanges the access token for a 12-hour session token itself.
(requires `pip install requests requests-oauthlib`).

Exit codes: 0 = updated (changes written), 3 = no changes, other = error.
"""

import argparse
import csv
import datetime as dt
import io
import json
import os
import shutil
import sys
import tempfile
import time
import urllib.request

API = "https://bigsdb.pasteur.fr/api/db/pubmlst_diphtheria_seqdef"
NO_CHANGES = 3


# --------------------------------------------------------------------------
# HTTP
# --------------------------------------------------------------------------

class Client:
    def __init__(self, base):
        self.base = base.rstrip("/")
        self.session = None
        keys = [os.environ.get(k) for k in (
            "BIGSDB_CONSUMER_KEY", "BIGSDB_CONSUMER_SECRET",
            "BIGSDB_ACCESS_TOKEN", "BIGSDB_ACCESS_SECRET")]
        if any(keys) and not all(keys):
            raise RuntimeError("some BIGSDB_* credentials are set but not all four")
        self.authenticated = all(keys)
        if self.authenticated:
            from requests_oauthlib import OAuth1Session
            # Access tokens are long-lived but only grant session tokens;
            # API resources need a session token (valid ~12 h).
            acc = OAuth1Session(keys[0], client_secret=keys[1],
                                resource_owner_key=keys[2],
                                resource_owner_secret=keys[3])
            r = acc.get(f"{self.base}/oauth/get_session_token", timeout=60)
            if r.status_code != 200:
                raise RuntimeError(
                    f"could not get session token ({r.status_code}): {r.text[:300]} "
                    "- access token may be revoked; re-run scripts/bigsdb_auth.py")
            tok = r.json()
            self.session = OAuth1Session(keys[0], client_secret=keys[1],
                                         resource_owner_key=tok["oauth_token"],
                                         resource_owner_secret=tok["oauth_token_secret"])
            print("Authenticated to BIGSdb-Pasteur (session token obtained)",
                  file=sys.stderr)

    def get(self, path_or_url, retries=4):
        url = path_or_url if path_or_url.startswith("http") \
            else f"{self.base}/{path_or_url.lstrip('/')}"
        for attempt in range(retries):
            try:
                if self.session:
                    r = self.session.get(url, timeout=300)
                    r.raise_for_status()
                    return r.text
                req = urllib.request.Request(
                    url, headers={"User-Agent": "cdiphtheriae-typing-updater"})
                with urllib.request.urlopen(req, timeout=300) as resp:
                    return resp.read().decode("utf-8")
            except Exception as e:  # noqa: BLE001
                if attempt == retries - 1:
                    raise RuntimeError(f"GET {url} failed: {e}") from e
                time.sleep(5 * (attempt + 1))

    def json(self, path_or_url):
        return json.loads(self.get(path_or_url))


# --------------------------------------------------------------------------
# Fetching
# --------------------------------------------------------------------------

def find_mlst_scheme(client, name):
    data = client.json("schemes")
    for s in data.get("schemes", []):
        if s.get("description", "").strip().lower() == name.lower():
            return s["scheme"]
    names = [s.get("description") for s in data.get("schemes", [])]
    raise RuntimeError(f"scheme '{name}' not found; available: {names}")


def normalise_fasta(text, locus):
    """Rewrite headers to >{locus}_{allele_id}; return (fasta, ids)."""
    out, ids = [], []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            raw = line[1:].split()[0]
            allele = raw.rsplit("_", 1)[-1]
            ids.append(allele)
            out.append(f">{locus}_{allele}")
        else:
            out.append(line.upper())
    if not ids:
        raise RuntimeError(f"no alleles returned for locus {locus}")
    return "\n".join(out) + "\n", ids


def fetch(client, staging, scheme_name, tox_locus):
    scheme_url = find_mlst_scheme(client, scheme_name)
    scheme = client.json(scheme_url)
    loci = [u.rstrip("/").split("/")[-1] for u in scheme["loci"]]
    print(f"MLST scheme: {scheme_url} loci={loci}", file=sys.stderr)

    os.makedirs(os.path.join(staging, "mlst", "alleles"))
    os.makedirs(os.path.join(staging, "tox"))
    counts = {}
    for locus in loci:
        fasta, ids = normalise_fasta(
            client.get(f"loci/{locus}/alleles_fasta"), locus)
        with open(os.path.join(staging, "mlst", "alleles", f"{locus}.fasta"), "w") as fh:
            fh.write(fasta)
        counts[locus] = len(ids)

    # profiles_csv is tab-delimited despite the name
    prof = client.get(f"{scheme_url}/profiles_csv")
    rows = list(csv.reader(io.StringIO(prof), delimiter="\t"))
    header = rows[0]
    missing = [l for l in loci if l not in header]
    if header[0] != "ST" or missing:
        raise RuntimeError(f"unexpected profiles header {header[:10]}")
    with open(os.path.join(staging, "mlst", "profiles.tsv"), "w", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(
            [["ST"] + loci] + [[r[0]] + [r[header.index(l)] for l in loci]
                               for r in rows[1:] if r and r[0]])
    counts["profiles"] = len(rows) - 1

    fasta, ids = normalise_fasta(client.get(f"loci/{tox_locus}/alleles_fasta"), "tox")
    with open(os.path.join(staging, "tox", "tox_alleles.fasta"), "w") as fh:
        fh.write(fasta)
    counts["tox"] = len(ids)
    return loci, counts


# --------------------------------------------------------------------------
# Comparison
# --------------------------------------------------------------------------

def fasta_records(path):
    recs, name = {}, None
    if not os.path.exists(path):
        return recs
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                recs[name] = []
            elif name:
                recs[name].append(line)
    return {k: "".join(v) for k, v in recs.items()}


def profiles(path):
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r)
        return {row[0]: tuple(row[1:]) for row in r if row}


def compare(old_db, new_db, loci):
    changes = []
    files = [("tox", os.path.join("tox", "tox_alleles.fasta"))] + \
            [(l, os.path.join("mlst", "alleles", f"{l}.fasta")) for l in loci]
    for label, rel in files:
        old, new = fasta_records(os.path.join(old_db, rel)), \
            fasta_records(os.path.join(new_db, rel))
        added = sorted(set(new) - set(old))
        removed = sorted(set(old) - set(new))
        altered = sorted(k for k in set(old) & set(new) if old[k] != new[k])
        changes.append((label, len(old), len(new), added, removed, altered))
    op = profiles(os.path.join(old_db, "mlst", "profiles.tsv"))
    np_ = profiles(os.path.join(new_db, "mlst", "profiles.tsv"))
    changes.append(("STs", len(op), len(np_),
                    sorted(set(np_) - set(op), key=lambda x: int(x) if x.isdigit() else 0),
                    sorted(set(op) - set(np_)),
                    sorted(k for k in set(op) & set(np_) if op[k] != np_[k])))
    return changes


def short(items, n=15):
    if not items:
        return "-"
    s = ", ".join(items[:n])
    return s + (f" ... (+{len(items) - n})" if len(items) > n else "")


def write_summary(path, version, changes, authenticated, api=API):
    lines = [f"# Database update: {version}", "",
             f"Source: BIGSdb-Pasteur `{api}`  ",
             f"Retrieved: {dt.datetime.now(dt.timezone.utc):%Y-%m-%d %H:%M UTC}", ""]
    if not authenticated:
        lines += ["> **WARNING: unauthenticated download.** BIGSdb-Pasteur withholds "
                  "data curated after 1 January 2025 from anonymous API users. "
                  "This update is likely INCOMPLETE. Configure the BIGSDB_* secrets "
                  "before approving.", ""]
    lines += ["| Set | Before | After | Added | Removed | Sequence/profile changed |",
              "|---|---:|---:|---|---|---|"]
    for label, b, a, add, rem, alt in changes:
        lines.append(f"| {label} | {b} | {a} | {short(add)} | {short(rem)} | {short(alt)} |")
    lines += ["", "## Reviewer checklist", "",
              "- [ ] No unexpected removals or changed sequences/profiles",
              "- [ ] CI tests pass",
              "- [ ] Validation set re-run and results unchanged for existing types",
              "- [ ] Tag a release after merge", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


# --------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", default="database")
    p.add_argument("--summary", default="CHANGES_DB.md")
    p.add_argument("--scheme", default="MLST", help="scheme description to use")
    p.add_argument("--tox-locus", default="tox", help="BIGSdb locus name for tox")
    p.add_argument("--version-label",
                   default=f"Pasteur_{dt.date.today():%Y-%m}")
    p.add_argument("--api", default=API)
    args = p.parse_args()

    client = Client(args.api)
    if not client.authenticated:
        print("WARNING: no BIGSDB_* credentials - post-2025 data will be missing",
              file=sys.stderr)

    staging = tempfile.mkdtemp(prefix="dbupdate_")
    try:
        loci, counts = fetch(client, staging, args.scheme, args.tox_locus)
        print(f"Downloaded: {counts}", file=sys.stderr)
        changes = compare(args.db, staging, loci)
        if not any(c[3] or c[4] or c[5] for c in changes):
            print("No changes against current database.", file=sys.stderr)
            return NO_CHANGES
        with open(os.path.join(staging, "VERSION.txt"), "w") as fh:
            fh.write(f"{args.version_label}\n")
        # Swap in the new data (keep any README in the db dir)
        for sub in ("tox", "mlst"):
            shutil.rmtree(os.path.join(args.db, sub), ignore_errors=True)
            shutil.copytree(os.path.join(staging, sub), os.path.join(args.db, sub))
        shutil.copy(os.path.join(staging, "VERSION.txt"), args.db)
        write_summary(args.summary, args.version_label, changes, client.authenticated, args.api)
        print(f"Database updated to {args.version_label}", file=sys.stderr)
        return 0
    finally:
        shutil.rmtree(staging, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
