#!/bin/bash
# R prepares/validates the graph inside its container; bsub runs on the HOST.
set -eo pipefail
umask 077
ROOT="${SCRIPTS_DIR:?Set SCRIPTS_DIR to the scripts directory}"
output="$("$ROOT/bin/Rscript_sif" "$ROOT/run_pipeline.R" "$@")"
printf '%s\n' "$output"
submit=no
for arg in "$@"; do [[ "$arg" != --submit ]] || submit=yes; done
if [[ "$submit" == yes ]]; then
  last="${output##*$'\n'}"
  [[ "$last" == HOST_PLAN=* ]] || { printf 'No host submission plan returned\n' >&2; exit 1; }
  plan="${last#HOST_PLAN=}"
  [[ -f "$plan" ]] || { printf 'Host plan missing\n' >&2; exit 1; }
  bash "$plan"
fi
