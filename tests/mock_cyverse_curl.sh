#!/bin/bash
set -eo pipefail
output="";url=""
while (($#)); do
  case "$1" in
    --output) output="$2";shift 2 ;;
    --max-filesize|--max-time|--connect-timeout|--proto|--proto-redir|--continue-at) shift 2 ;;
    --*) shift ;;
    *) url="$1";shift ;;
  esac
done
[[ -n "$output" && -n "$url" ]] || exit 11
[[ "$url" != */MD5.txt ]] || { echo 'Unexpected checksum-catalog request' >&2; exit 12; }
[[ -z "${MOCK_CYVERSE_COUNTER:-}" ]] || printf '%s\n' "$url" >> "$MOCK_CYVERSE_COUNTER"
cp -- "${MOCK_CYVERSE_DIR:?}/merged_flt_c1.imputed.vcf.gz" "$output"
