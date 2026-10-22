#!/usr/bin/env bash
# Assemble a Galaxy-installable tool directory with the database bundled:
#   bash scripts/bundle_tool.sh [outdir]   -> <outdir>/cdiphtheriae_typing/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$ROOT/build}"
DEST="$OUT/cdiphtheriae_typing"
[ -f "$ROOT/database/mlst/profiles.tsv" ] || {
  echo "database/ is empty - run scripts/update_database.py first" >&2; exit 1; }
rm -rf "$DEST" && mkdir -p "$DEST"
cp "$ROOT"/tool/tool.xml "$ROOT"/tool/diphtheria_typing.py "$DEST"/
cp -r "$ROOT"/tool/test-data "$DEST"/
cp -r "$ROOT"/database "$DEST"/database
echo "Bundled tool at $DEST (database $(grep -v '^#' "$ROOT/database/VERSION.txt" | head -1))"
