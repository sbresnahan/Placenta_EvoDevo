#!/bin/bash
# Host-only CyVerse transfer. Obtain each file's published MD5 from MD5.txt,
# resume incomplete HTTPS transfers and pin the verified result in a receipt.
set -eo pipefail
umask 077
[[ $# == 9 ]] || { printf 'Usage: fetch_cyverse_md5.sh CURL URL DEST PARTIAL MD5_URL ATTEMPTS TIMEOUT CONNECT_TIMEOUT RETRY_SECONDS\n' >&2; exit 2; }
curl_bin="$1"; url="$2"; dest="$3"; partial="$4"; md5_url="$5"
attempts="$6"; timeout="$7"; connect_timeout="$8"; retry_seconds="$9"
[[ "$curl_bin" == /* && -x "$curl_bin" &&
   "$url" == https://data.cyverse.org/dav-anon/* &&
   "$md5_url" == "${url%/*}/MD5.txt" &&
   "$dest" == /* && "$partial" == /* && "$dest" != "$partial" &&
   "$attempts" =~ ^[1-9][0-9]*$ && "$timeout" =~ ^[1-9][0-9]*$ &&
   "$connect_timeout" =~ ^[1-9][0-9]*$ &&
   "$retry_seconds" =~ ^[0-9]+([.][0-9]+)?$ ]] || {
  printf 'Invalid CyVerse retrieval arguments\n' >&2; exit 2;
}
name="${url##*/}"
[[ "$name" =~ ^[A-Za-z0-9_.-]+\.vcf\.gz$ ]] || { printf 'Unsafe CyVerse VCF basename\n' >&2; exit 2; }
mkdir -p -- "$(dirname "$dest")" "$(dirname "$partial")"
receipt="${dest}.cyverse.tsv"
catalog="$(mktemp "${dest}.md5_catalog.XXXXXX")"
trap 'rm -f -- "$catalog"' EXIT
max_bytes=15000000000 # Never accept an unbounded or wrong HTTP response.

expected=""
for ((i=1;i<=attempts;i++)); do
  if "$curl_bin" --fail --location --silent --show-error \
      --proto '=https' --proto-redir '=https' --max-filesize 100000 \
      --max-time 120 --connect-timeout "$connect_timeout" \
      --output "$catalog" "$md5_url"; then
    # Accept standard md5sum, filename-first, and "MD5 (name) = hash"
    # catalogs, but never accept a hash unless its row names our VCF.
    expected="$(awk -v target="$name" '
      {sub(/\r$/, "");for(i=1;i<=NF;i++){
         file=$i;gsub(/^[*(]+|[):]+$/, "", file);sub(/^.*\//, "", file)
         if(file!=target)continue
         for(j=1;j<=NF;j++){
           hash=$j;gsub(/^[^[:xdigit:]]+|[^[:xdigit:]]+$/, "", hash)
           if(length(hash)==32&&hash~/^[[:xdigit:]]+$/){print tolower(hash);next}
         }
       }}' "$catalog")"
    if [[ "$expected" =~ ^[0-9a-f]{32}$ ]]; then break; fi
    printf 'Invalid CyVerse MD5.txt or IP-verification page for %s; refusing unverified data\n' "$name" >&2
    exit 1
  fi
  if (( i==attempts )); then
    printf 'CyVerse MD5 catalog unavailable after %s attempts: %s\n' "$attempts" "$md5_url" >&2
    exit 1
  fi
  sleep "$retry_seconds"
done

verify() {
  [[ -f "$1" ]] || return 1
  local size hash
  size="$(stat -c %s -- "$1")"
  (( size > 0 && size <= max_bytes )) || return 1
  hash="$(md5sum -- "$1")"
  [[ "${hash:0:32}" == "$expected" ]]
}
record() {
  local size temp
  size="$(stat -c %s -- "$dest")"
  temp="$(mktemp "${receipt}.XXXXXX")"
  printf 'path\tbytes\tmd5\turl\tmethod\n%s\t%s\t%s\t%s\tcyverse_md5\n' \
    "$dest" "$size" "$expected" "$url" > "$temp"
  mv -f -- "$temp" "$receipt"
}
if verify "$dest"; then
  record
  printf 'Reusing MD5-confirmed CyVerse VCF: %s\n' "$name"
  exit 0
fi

for ((i=1;i<=attempts;i++)); do
  if [[ -f "$partial" ]]; then
    size="$(stat -c %s -- "$partial")"
    if (( size >= max_bytes )); then
      if verify "$partial"; then mv -f -- "$partial" "$dest"; record; exit 0; fi
      rm -f -- "$partial"
    fi
  fi
  offset=0
  if [[ -f "$partial" ]]; then offset="$(stat -c %s -- "$partial")"; fi
  printf 'CyVerse %s: attempt %s/%s, %s bytes already staged\n' "$name" "$i" "$attempts" "$offset" >&2
  if "$curl_bin" --fail --location --silent --show-error \
      --proto '=https' --proto-redir '=https' --max-filesize "$max_bytes" \
      --max-time "$timeout" --connect-timeout "$connect_timeout" \
      --continue-at - --output "$partial" "$url"; then
    if verify "$partial"; then
      mv -f -- "$partial" "$dest"
      record
      printf 'CyVerse VCF staged with published MD5: %s\n' "$name"
      exit 0
    fi
    # A complete but wrong response (including an IP-unblocking HTML page)
    # must not be resumed or promoted on a subsequent attempt.
    rm -f -- "$partial"
    printf 'CyVerse response did not match published MD5: %s\n' "$name" >&2
  fi
  if (( i < attempts )); then sleep "$retry_seconds"; fi
done
printf 'CyVerse retrieval exhausted %s attempts: %s\n' "$attempts" "$name" >&2
exit 1
