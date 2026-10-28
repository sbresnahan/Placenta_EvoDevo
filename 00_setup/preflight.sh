#!/bin/bash
set -eo pipefail
umask 077
ROOT="${SCRIPTS_DIR:?Set SCRIPTS_DIR to the scripts directory}"
CONFIG="${1:-${CONFIG:-}}"
[[ "$CONFIG" == /* && -f "$CONFIG" ]] || { printf 'Provide an absolute cohort YAML\n' >&2; exit 1; }
export CONFIG
source "$ROOT/bin/seadragon_init.sh"
for tool in bash gzip curl bsub "$POE_SINGULARITY" /usr/bin/time; do
  command -v "$tool" >/dev/null 2>&1 || {
    printf 'Required tool unavailable: %s\n' "$tool" >&2; exit 1;
  }
done
rs="$ROOT/bin/Rscript_sif"
[[ $# -le 2 && ( -z "${2:-}" || "${2:-}" == --execute ) ]] || {
  printf 'Usage: 00_setup/preflight.sh CONFIG [--execute]\n' >&2; exit 1;
}
dest="$("$rs" "$ROOT/runtime.R" "$CONFIG" preflight inventory-path shared unused)"
export POE_NATIVE_RECORD="${dest}.pending.$$"
trap 'rm -f "$POE_NATIVE_RECORD"' EXIT
"$ROOT/bin/host_inventory.sh" "$POE_NATIVE_RECORD"
"$rs" "$ROOT/runtime.R" "$CONFIG" preflight record-inventory shared unused
