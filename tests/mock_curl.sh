#!/bin/bash
# Offline transport fixture for test_fetch_host.sh; never accesses the network.
set -eo pipefail
out=""; max=""; url=""
while (( $# )); do
  case "$1" in
    --output|--max-filesize|--max-time|--connect-timeout|--proto|--proto-redir|--continue-at)
      key="$1"; val="$2"; shift 2
      if [[ "$key" == --output ]]; then out="$val"; fi
      if [[ "$key" == --max-filesize ]]; then max="$val"; fi ;;
    --fail|--location|--silent|--show-error) shift ;;
    https://*) url="$1"; shift ;;
    *) exit 40 ;;
  esac
done
[[ -n "$out" && -n "$url" && "$max" == "${MOCK_EXPECT_SIZE:?}" ]] || exit 41
[[ -f "${MOCK_SOURCE:?}" ]] || exit 42
offset=0
if [[ -f "$out" ]]; then offset="$(stat -c %s -- "$out")"; fi
tail -c "+$((offset+1))" "$MOCK_SOURCE" >> "$out"
printf 'called\n' >> "${MOCK_COUNTER:?}"
