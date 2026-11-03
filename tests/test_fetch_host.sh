#!/bin/bash
# Run directly with Bash. No R, network, or LSF access is required.
set -eo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fixture="$(mktemp -d)"
trap 'rm -rf -- "$fixture"' EXIT
source_file="$fixture/source"; dest="$fixture/destination with spaces"; partial="$fixture/partial"
printf 'abcdef' > "$source_file"
export MOCK_SOURCE="$source_file" MOCK_EXPECT_SIZE=6 MOCK_COUNTER="$fixture/counter"
printf 'abc' > "$partial"
fetch=("$root/bin/fetch_verified.sh" "$root/tests/mock_curl.sh" \
  "https://example.invalid/fixture" "$dest" "$partial" 6 2 10 2 0)
"${fetch[@]}"
cmp -s "$source_file" "$dest"
[[ ! -e "$partial" && "$(wc -l < "$MOCK_COUNTER")" == 1 ]]
"${fetch[@]}" # Verified destination must be reused without another request.
[[ "$(wc -l < "$MOCK_COUNTER")" == 1 ]]
printf 'corrupt' > "$dest"
printf 'abc' > "$partial"
"${fetch[@]}" # Invalid-size destination is replaced with the resumed, correctly sized transfer.
cmp -s "$source_file" "$dest"
[[ "$(wc -l < "$MOCK_COUNTER")" == 2 ]]
if "$root/bin/fetch_verified.sh" "$root/tests/mock_curl.sh" \
     "http://example.invalid/fixture" "$dest" "$partial" 6 1 10 2 0 2>/dev/null; then
  printf 'Non-HTTPS source unexpectedly accepted\n' >&2; exit 1
fi
printf 'PASS: offline host-only retrieval, resume, size checks, and HTTPS validation\n'
