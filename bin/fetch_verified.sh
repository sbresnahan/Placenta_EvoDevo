#!/bin/bash
# Host-only, resumable HTTPS acquisition with exact byte/MD5 verification.
set -eo pipefail
umask 077
[[ $# == 10 ]] || { printf 'Usage: fetch_verified.sh CURL URL DEST PARTIAL BYTES MD5 ATTEMPTS TIMEOUT CONNECT_TIMEOUT RETRY_SECONDS\n' >&2; exit 2; }
curl_bin="$1"; url="$2"; dest="$3"; partial="$4"; bytes="$5"; md5="$6"
attempts="$7"; timeout="$8"; connect_timeout="$9"; retry_seconds="${10}"
[[ -x "$curl_bin" && "$curl_bin" == /* && "$url" == https://* &&
   "$dest" == /* && "$partial" == /* && "$dest" != "$partial" &&
   "$bytes" =~ ^[1-9][0-9]*$ && "$md5" =~ ^[0-9a-f]{32}$ &&
   "$attempts" =~ ^[1-9][0-9]*$ && "$timeout" =~ ^[1-9][0-9]*$ &&
   "$connect_timeout" =~ ^[1-9][0-9]*$ &&
   "$retry_seconds" =~ ^[0-9]+([.][0-9]+)?$ ]] || {
  printf 'Invalid host retrieval arguments\n' >&2; exit 2;
}
verify() {
  [[ -f "$1" ]] || return 1
  local size hash
  size="$(stat -c %s -- "$1")"
  [[ "$size" == "$bytes" ]] || return 1
  hash="$(md5sum -- "$1")"
  [[ "${hash:0:32}" == "$md5" ]]
}
if verify "$dest"; then exit 0; fi
mkdir -p -- "$(dirname "$dest")" "$(dirname "$partial")"
for ((i=1;i<=attempts;i++)); do
  if [[ -f "$partial" ]]; then
    size="$(stat -c %s -- "$partial")"
    if [[ "$size" -ge "$bytes" ]]; then
      if verify "$partial"; then mv -f -- "$partial" "$dest"; exit 0; fi
      # Full-sized but invalid partials cannot be resumed safely.
      rm -f -- "$partial"
    fi
  fi
  # -C - resumes from the partial size and rejects servers ignoring Range.
  # The CLI's --max-filesize accepts sizes beyond R's 32-bit curl option.
  if "$curl_bin" --fail --location --silent --show-error \
      --proto '=https' --proto-redir '=https' \
      --max-filesize "$bytes" --max-time "$timeout" \
      --connect-timeout "$connect_timeout" \
      --continue-at - --output "$partial" "$url"; then
    if verify "$partial"; then mv -f -- "$partial" "$dest"; exit 0; fi
  fi
  if (( i < attempts )); then sleep "$retry_seconds"; fi
done
printf 'Host acquisition exhausted %s attempts: %s\n' "$attempts" "$dest" >&2
exit 1
