#!/usr/bin/env bash
# Download the genomes listed in validation/expected.tsv from NCBI as
# <outdir>/<sample>.fasta.gz (skips files that already exist) and write
# <outdir>/MANIFEST.tsv.
#   bash scripts/fetch_validation_genomes.sh [outdir]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$ROOT/validation/genomes}"
mkdir -p "$OUT"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

tail -n +2 "$ROOT/validation/expected.tsv" | while IFS=$'\t' read -r sample acc _rest; do
  [ -z "$acc" ] && continue
  dest="$OUT/$sample.fasta.gz"
  if [ -s "$dest" ]; then echo "have $sample"; continue; fi
  case "$acc" in
    GCF_*|GCA_*)
      curl -fsSL --retry 4 -o "$TMP/$acc.zip" \
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA"
      unzip -qo "$TMP/$acc.zip" -d "$TMP/$acc"
      cat "$TMP/$acc"/ncbi_dataset/data/"$acc"/*.fna > "$TMP/$sample.fasta" ;;
    *)
      curl -fsSL --retry 4 -o "$TMP/$sample.fasta" \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta&retmode=text"
      sleep 1 ;;
  esac
  grep -q '^>' "$TMP/$sample.fasta" || { echo "ERROR: $acc did not return FASTA" >&2; exit 1; }
  gzip -n -9 -c "$TMP/$sample.fasta" > "$dest"   # -n: reproducible output
  echo "fetched $sample ($acc)"
done

# Manifest (no dates, so it only changes when the files do)
{
  printf 'file\taccession\tsequences\tlength_bp\tsha256_uncompressed\tdescription\n'
  tail -n +2 "$ROOT/validation/expected.tsv" | while IFS=$'\t' read -r sample acc _rest; do
    f="$OUT/$sample.fasta.gz"; [ -s "$f" ] || continue
    n=$(gzip -dc "$f" | grep -c '^>')
    bp=$(gzip -dc "$f" | grep -v '^>' | tr -d '\n\r' | wc -c)
    sha=$(gzip -dc "$f" | sha256sum | cut -d' ' -f1)
    desc=$(gzip -dc "$f" | head -1 | cut -c2- | tr '\t' ' ')
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample.fasta.gz" "$acc" "$n" "$bp" "$sha" "$desc"
  done
} > "$OUT/MANIFEST.tsv"
echo "Wrote $OUT/MANIFEST.tsv"
