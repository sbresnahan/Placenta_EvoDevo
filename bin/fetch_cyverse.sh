#!/bin/bash
# Host-only CyVerse VCF acquisition, without MD5 or checksum catalogs.
# Reuses an original-named file from reference/cyverse_downloads/ when present.
set -eo pipefail
umask 077
[[ $# == 8 ]] || {
  printf 'Usage: fetch_cyverse.sh CURL URL DEST PARTIAL ATTEMPTS TIMEOUT CONNECT_TIMEOUT RETRY_SECONDS\n' >&2
  exit 2
}
curl_bin="$1"; url="$2"; dest="$3"; partial="$4"
attempts="$5"; timeout="$6"; connect_timeout="$7"; retry_seconds="$8"
[[ "$curl_bin" == /* && -x "$curl_bin" &&
   "$url" == https://data.cyverse.org/dav-anon/* &&
   "$dest" == /* && "$partial" == /* && "$dest" != "$partial" &&
   "$attempts" =~ ^[1-9][0-9]*$ && "$timeout" =~ ^[1-9][0-9]*$ &&
   "$connect_timeout" =~ ^[1-9][0-9]*$ &&
   "$retry_seconds" =~ ^[0-9]+([.][0-9]+)?$ ]] || {
  printf 'Invalid CyVerse retrieval arguments\n' >&2; exit 2
}
name="${url##*/}"
[[ "$name" =~ ^[A-Za-z0-9_.-]+\.vcf\.gz$ ]] || {
  printf 'Unsafe CyVerse VCF basename\n' >&2; exit 2
}
mkdir -p -- "$(dirname "$dest")" "$(dirname "$partial")"
staged="$(dirname "$dest")/cyverse_downloads/$name"
receipt="${dest}.cyverse.tsv"
max_bytes=15000000000

valid_vcf() {
  [[ -f "$1" && ! -L "$1" ]] || return 1
  local bytes
  bytes="$(stat -c %s -- "$1")"
  (( bytes > 0 && bytes <= max_bytes )) || return 1
  gzip -t -- "$1"
}
record() {
  local bytes temp
  bytes="$(stat -c %s -- "$dest")"
  temp="$(mktemp "${receipt}.XXXXXXXX")"
  printf 'path\tbytes\turl\tmethod\n%s\t%s\t%s\tcyverse\n' \
    "$dest" "$bytes" "$url" > "$temp"
  mv -- "$temp" "$receipt"
}

if valid_vcf "$dest"; then
  record
  printf 'Reusing gzip-valid CyVerse VCF: %s\n' "$name"
  exit 0
fi
[[ ! -e "$dest" && ! -L "$dest" ]] || {
  printf 'Existing invalid destination; refusing overwrite: %s\n' "$dest" >&2
  exit 1
}
if [[ -e "$staged" || -L "$staged" ]]; then
  valid_vcf "$staged" || {
    printf 'Staged VCF is empty, too large, or gzip-invalid: %s\n' "$staged" >&2
    exit 1
  }
  mv -n -- "$staged" "$dest"
  [[ ! -e "$staged" && ! -L "$staged" && -f "$dest" ]] || {
    printf 'Unable to promote staged VCF safely: %s\n' "$dest" >&2
    exit 1
  }
  record
  printf 'Promoted locally staged CyVerse VCF: %s\n' "$name"
  exit 0
fi

for ((i=1;i<=attempts;i++)); do
  printf 'CyVerse %s: HTTPS attempt %s/%s\n' "$name" "$i" "$attempts" >&2
  if "$curl_bin" --fail --location --silent --show-error \
      --proto '=https' --proto-redir '=https' --max-filesize "$max_bytes" \
      --max-time "$timeout" --connect-timeout "$connect_timeout" \
      --continue-at - --output "$partial" "$url"; then
    if valid_vcf "$partial"; then
      mv -n -- "$partial" "$dest"
      [[ ! -e "$partial" && -f "$dest" ]] || {
        printf 'Destination appeared during VCF promotion\n' >&2; exit 1
      }
      record
      printf 'Downloaded gzip-valid CyVerse VCF: %s\n' "$name"
      exit 0
    fi
    # A complete HTTP response that is not gzip-valid cannot be resumed.
    rm -f -- "$partial"
    printf 'Downloaded VCF failed gzip integrity test: %s\n' "$name" >&2
  fi
  if (( i < attempts )); then sleep "$retry_seconds"; fi
done
printf 'CyVerse retrieval exhausted %s attempts: %s\n' "$attempts" "$name" >&2
exit 1
