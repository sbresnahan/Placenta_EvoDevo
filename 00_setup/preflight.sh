#!/bin/bash
set -eo pipefail
umask 077
ROOT="${SCRIPTS_DIR:?Set SCRIPTS_DIR to the repository root}"
CONFIG="${1:-${CONFIG:-}}"
[[ "$CONFIG" == /* && -f "$CONFIG" ]] || { printf 'Provide an absolute cohort YAML\n' >&2; exit 1; }
export CONFIG
source "$ROOT/bin/seadragon_init.sh"
rs="$ROOT/bin/Rscript_sif"
if [[ "${2:-}" == --execute ]]; then
  dest="$("$rs" "$ROOT/runtime.R" "$CONFIG" preflight inventory-path shared unused)"
  export POE_NATIVE_RECORD="${dest}.pending.$$"
  "$ROOT/bin/host_inventory.sh" "$POE_NATIVE_RECORD"
  printf 'Host inventory for local review: %s\n' "$POE_NATIVE_RECORD"
  "$rs" "$ROOT/runtime.R" "$CONFIG" preflight accept-inventory shared unused
  rm -f "$POE_NATIVE_RECORD"
  unset POE_NATIVE_RECORD
fi
"$rs" "$ROOT/poe.R" "$CONFIG" preflight
