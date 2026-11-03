#!/bin/bash
# Offline CyVerse tests: no MD5 catalog, network, or real VCF required.
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf -- "$work"' EXIT
export MOCK_CYVERSE_DIR="$work/mock" MOCK_CYVERSE_COUNTER="$work/calls"
mkdir -p "$MOCK_CYVERSE_DIR" "$work/out"
name=merged_flt_c1.imputed.vcf.gz
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' |
  gzip -c > "$MOCK_CYVERSE_DIR/$name"
dest="$work/out/chr1.vcf.gz"
url="https://data.cyverse.org/dav-anon/iplant/projects/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed/$name"
args=("$ROOT/tests/mock_cyverse_curl.sh" "$url" "$dest" "$dest.partial" 2 30 5 0)
bash "$ROOT/bin/fetch_cyverse.sh" "${args[@]}"
[[ -s "$dest.cyverse.tsv" && -s "$dest" ]]
gzip -t "$dest"
[[ "$(wc -l < "$MOCK_CYVERSE_COUNTER")" == 1 ]]
bash "$ROOT/bin/fetch_cyverse.sh" "${args[@]}"
[[ "$(wc -l < "$MOCK_CYVERSE_COUNTER")" == 1 ]]

# The local staging route must not attempt a network connection.
mkdir -p "$work/local/cyverse_downloads"
cp "$MOCK_CYVERSE_DIR/$name" "$work/local/cyverse_downloads/$name"
local_dest="$work/local/chr1.vcf.gz"
bash "$ROOT/bin/fetch_cyverse.sh" /usr/bin/false "$url" "$local_dest" \
  "$local_dest.partial" 1 30 5 0
[[ -s "$local_dest.cyverse.tsv" && ! -e "$work/local/cyverse_downloads/$name" ]]

# Gzip-invalid content must remain in staging rather than being promoted.
mkdir -p "$work/invalid/cyverse_downloads"
printf '<html>verification required</html>\n' > "$work/invalid/cyverse_downloads/$name"
bad_dest="$work/invalid/chr1.vcf.gz"
if bash "$ROOT/bin/fetch_cyverse.sh" /usr/bin/false "$url" "$bad_dest" \
  "$bad_dest.partial" 1 30 5 0 >/dev/null 2>&1; then
  printf 'FAIL: non-gzip staged VCF was accepted\n' >&2; exit 1
fi
[[ -s "$work/invalid/cyverse_downloads/$name" && ! -e "$bad_dest" ]]
echo 'PASS: offline CyVerse gzip integrity, receipt reuse, local staging, and invalid-response rejection'
