#!/bin/bash
set -eo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf -- "$work"' EXIT
export MOCK_CYVERSE_DIR="$work/mock"
mkdir -p "$MOCK_CYVERSE_DIR" "$work/out"
name=merged_flt_c1.imputed.vcf.gz
printf 'test compressed bytes\n' > "$MOCK_CYVERSE_DIR/$name"
md5sum "$MOCK_CYVERSE_DIR/$name" | sed "s|$MOCK_CYVERSE_DIR/||" > "$MOCK_CYVERSE_DIR/MD5.txt"
dest="$work/out/chr1.vcf.gz"
url="https://data.cyverse.org/dav-anon/iplant/projects/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed/$name"
args=("$ROOT/tests/mock_cyverse_curl.sh" "$url" "$dest" "$dest.partial" "${url%/*}/MD5.txt" 2 30 5 0)
bash "$ROOT/bin/fetch_cyverse_md5.sh" "${args[@]}"
[[ -s "$dest.cyverse.tsv" && -s "$dest" ]]
bash "$ROOT/bin/fetch_cyverse_md5.sh" "${args[@]}"
printf 'deadbeefdeadbeefdeadbeefdeadbeef  %s\n' "$name" > "$MOCK_CYVERSE_DIR/MD5.txt"
if bash "$ROOT/bin/fetch_cyverse_md5.sh" "${args[@]}" >/dev/null 2>&1; then
  printf 'FAIL: wrong publisher MD5 must be rejected\n' >&2;exit 1
fi
printf '<html>verification required</html>\n' > "$MOCK_CYVERSE_DIR/MD5.txt"
if bash "$ROOT/bin/fetch_cyverse_md5.sh" "${args[@]}" >/dev/null 2>&1; then
  printf 'FAIL: CyVerse HTML/IP-verification page must be rejected\n' >&2;exit 1
fi
echo 'PASS: offline CyVerse MD5 catalog, receipt reuse, mismatch and HTML rejection'
