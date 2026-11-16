#!/bin/bash
# Scheduler mock only: no LSF submission, container, or RNA processing.
set -euo pipefail
ROOT="$(cd "${1:-.}" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$WORK/root with space/bin" "$WORK/logs" "$WORK/mockbin"
cp "$ROOT/bin/submit_exposure_nb.sh" "$WORK/root with space/bin/"
printf '#!/bin/bash\nprintf "%%s\\n" "$@" > "$MOCK_RECEIPT"\n' > "$WORK/root with space/bin/bsub_site"
chmod +x "$WORK/root with space/bin/bsub_site"
touch "$WORK/config with space.yaml"
export SCRIPTS_DIR="$WORK/root with space" EX_LOG_DIR="$WORK/logs"
export EX_MEMORY_GB=2 EX_WALL=00:10 EX_ARRAY_END=3 EX_MAX_CONCURRENT=1
export MOCK_RECEIPT="$WORK/receipt"
bash "$ROOT/bin/submit_exposure_nb.sh" "$WORK/config with space.yaml" fit ARRAY_INDEX
grep -Fx 'ex_fit[1-3]%1' "$MOCK_RECEIPT"
script="$(tail -1 "$MOCK_RECEIPT")"
bash -n "$script"
grep -F 'LSB_JOBINDEX' "$script"
grep -F 'config\ with\ space.yaml' "$script"
if EX_CPUS=2 bash "$ROOT/bin/submit_exposure_nb.sh" "$WORK/config with space.yaml" fit 1 > /dev/null 2>&1; then
  echo "FAIL multithread overreservation accepted"; exit 1
fi
if bash "$ROOT/bin/submit_exposure_nb.sh" "$WORK/config with space.yaml" unknown > /dev/null 2>&1; then
  echo "FAIL unknown action accepted"; exit 1
fi
echo "PASS scheduler mock, array bounds, path quoting, CPU gate, unknown action"
