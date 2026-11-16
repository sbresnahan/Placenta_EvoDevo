#!/bin/bash
# Host-only public annotation transfer. Never source a downloaded script.
set -eo pipefail
umask 077
[[ $# == 3 ]] || { printf 'Usage: URL DEST FORMAT\n' >&2; exit 2; }
url="$1"; dest="$2"; format="$3"
[[ "$url" == https://* && "$dest" != *$'\n'* ]] || exit 2
case "$format" in rmsk|gff3|bed|sql|text|fasta_gz|h5_gz) ;; *) exit 2 ;; esac
mkdir -p "$(dirname "$dest")"
tmp="${dest}.partial"
validate() {
  local file="$1"
  [[ -s "$file" ]] || return 1
  case "$format" in
    rmsk|fasta_gz|h5_gz)
      gzip -t "$file" || return 1
      ;;
    gff3|bed)
      if [[ "$dest" == *.gz ]]; then gzip -t "$file" || return 1; fi
      ;;
  esac
  # HTTP error bodies cannot pass gzip, and non-gzip files are screened too.
  if [[ "$dest" != *.gz && "$format" != h5_gz && "$format" != fasta_gz && "$format" != rmsk ]]; then
    if head -c 1024 "$file" | grep -Eiq '<!doctype html|<html|access denied|forbidden'; then return 1; fi
  fi
  return 0
}
if [[ -e "$dest" ]]; then
  # Reuse only validated data. A malformed local/cache file is not replaced
  # implicitly because it may be an operator-supplied reference.
  validate "$dest" || { printf 'Existing resource invalid: %s\n' "$dest" >&2; exit 1; }
  exit 0
fi
trap 'rm -f "$tmp"' EXIT
curl --fail --location --silent --show-error --retry 2 --connect-timeout 30 \
  --max-time "${POE_ANNOTATION_TIMEOUT:-7200}" "$url" -o "$tmp"
validate "$tmp" || { printf 'Annotation validation failed: %s\n' "$url" >&2; exit 1; }
mv "$tmp" "$dest"
