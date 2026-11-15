#!/usr/bin/env bash
# Populate the cohort-local fish reference resource manifest using verified
# NCBI Content-Length headers. Run on a network-enabled Seadragon transfer host.
set -euo pipefail
: "${COHORT_DIR:?Export COHORT_DIR first}"
REF="$COHORT_DIR/reference"
MANIFEST="$COHORT_DIR/manifests/fish_Pprolifica_AB_resources.tsv"
BASE='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/474/105/GCF_027474105.1_Auburn_PproV2'
PREFIX='GCF_027474105.1_Auburn_PproV2'
[[ -f "$MANIFEST" ]] || { echo "Missing $MANIFEST" >&2; exit 1; }
[[ $(wc -l < "$MANIFEST") -eq 1 ]] || { echo "Existing resource rows: refusing overwrite of $MANIFEST" >&2; exit 1; }
mkdir -p "$REF" "$COHORT_DIR/scratch"
tmp="$(mktemp "${MANIFEST}.new.XXXXXX")"
trap 'rm -f "$tmp"' EXIT
printf 'resource_id\turl\tdestination\tpartial\tresponse\tbytes\tmd5\trequired\tbuild\tdecompress_to\treadiness\n' > "$tmp"
# HEAD may redirect; use the Content-Length of the FINAL response only.
get_bytes() {
  local url=$1 h n
  h=$(mktemp)
  if ! curl -fsSIL --retry 3 --retry-delay 5 --connect-timeout 30 --max-time 120 "$url" > "$h"; then
    rm -f "$h"; echo "Failed to verify $url" >&2; return 1
  fi
  n=$(awk '/^HTTP\//{v=""} tolower($1)=="content-length:"{gsub("\r","",$2); v=$2} END{print v}' "$h")
  rm -f "$h"
  [[ $n =~ ^[1-9][0-9]*$ ]] || { echo "Missing valid Content-Length: $url" >&2; return 1; }
  printf '%s' "$n"
}
add() {
  local id=$1 file=$2 output=$3 url bytes destination
  url="$BASE/$file"; bytes=$(get_bytes "$url")
  destination="$REF/$file"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t\tTRUE\tGCF_027474105.1\t%s\tready\n' \
    "$id" "$url" "$destination" "$destination.partial" \
    "$COHORT_DIR/scratch/${id}.response" "$bytes" "$output" >> "$tmp"
  printf 'Verified %-17s %s bytes\n' "$id" "$bytes"
}
add fish_genome "${PREFIX}_genomic.fna.gz" "$REF/genome.fa"
add fish_annotation "${PREFIX}_genomic.gtf.gz" "$REF/source_annotation.gtf"
add fish_assembly_report "${PREFIX}_assembly_report.txt" ''
# Atomic manifest creation, preserving the existing header-only original until success.
mv "$tmp" "$MANIFEST"
trap - EXIT
printf 'Wrote %s (3 resources)\n' "$MANIFEST"
