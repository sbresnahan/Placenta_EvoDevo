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
if [[ "$url" == */MD5.txt ]]; then
  cp -- "${MOCK_CYVERSE_DIR:?}/MD5.txt" "$output"
else
  cp -- "${MOCK_CYVERSE_DIR:?}/merged_flt_c1.imputed.vcf.gz" "$output"
fi
