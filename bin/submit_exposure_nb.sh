#!/bin/bash
# Explicit LSF submission only. No dependency chain or implicit production run.
set -euo pipefail
ROOT="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONFIG="${1:?Use CONFIG ACTION [INDEX]}"; ACTION="${2:?Specify ACTION}"; shift 2
[[ "$CONFIG" == /* && -f "$CONFIG" ]] || { echo "Use an existing absolute CONFIG" >&2; exit 1; }
case "$ACTION" in
  preflight|te-plan|te-import|prepare|pilot|fit|mediation|combine|bootstrap|enrichment|report-ready|export|status) ;;
  *) echo "Unknown exposure action" >&2; exit 1 ;;
esac
: "${EX_LOG_DIR:?Set an absolute private log directory}"
: "${EX_MEMORY_GB:?Set memory after a measured pilot}"
: "${EX_WALL:?Set wall time after a measured pilot}"
[[ "$EX_LOG_DIR" == /* ]] || { echo "Log directory must be absolute" >&2; exit 1; }
CPUS="${EX_CPUS:-1}"; MAX="${EX_MAX_CONCURRENT:-1}"; QUEUE="${EX_QUEUE:-medium}"
[[ "$EX_MEMORY_GB" =~ ^[1-9][0-9]*$ && "$CPUS" =~ ^[1-9][0-9]*$ && "$MAX" =~ ^[1-9][0-9]*$ ]] ||
  { echo "Positive integer memory, CPU and concurrency required" >&2; exit 1; }
[[ "$CPUS" == 1 ]] || { echo "NB worker requires one CPU; parallelize with arrays" >&2; exit 1; }
umask 077; mkdir -p "$EX_LOG_DIR"
script="$(mktemp "$EX_LOG_DIR/ex_${ACTION}.XXXXXX.sh")"
{
  printf '#!/bin/bash\nset -euo pipefail\numask 077\n'
  printf 'export SCRIPTS_DIR=%q\n' "$ROOT"
  printf 'exec /usr/bin/time -v %q %q %q %q' "$ROOT/bin/Rscript_sif" "$ROOT/poe_exposure.R" "$CONFIG" "$ACTION"
  for a in "$@"; do
    if [[ "$a" == ARRAY_INDEX ]]; then printf ' "${LSB_JOBINDEX:?Missing array index}"'
    else printf ' %q' "$a"; fi
  done
  printf '\n'
} > "$script"
NAME="ex_${ACTION}"
if [[ -n "${EX_ARRAY_END:-}" ]]; then
  [[ "$EX_ARRAY_END" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid EX_ARRAY_END" >&2; exit 1; }
  NAME="${NAME}[1-${EX_ARRAY_END}]%${MAX}"
fi
"$ROOT/bin/bsub_site" -q "$QUEUE" -J "$NAME" -n "$CPUS" -M "$EX_MEMORY_GB" \
  -R "rusage[mem=$EX_MEMORY_GB]" -W "$EX_WALL" \
  -oo "$EX_LOG_DIR/${ACTION}.%J.%I.out" -eo "$EX_LOG_DIR/${ACTION}.%J.%I.err" /bin/bash "$script"
