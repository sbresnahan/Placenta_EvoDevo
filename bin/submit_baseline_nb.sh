#!/bin/bash
# Host-side LSF helper. It submits only the requested job; no hidden job chain.
# Usage: submit_baseline_nb.sh CONFIG ACTION [ACTION_ARGS ...]
# Resources are explicit environment overrides, not guessed from feature count.
set -eo pipefail
umask 077
ROOT="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export SCRIPTS_DIR="$ROOT"
CONFIG="${1:?Usage: submit_baseline_nb.sh CONFIG ACTION [ARGS ...]}"
ACTION="${2:?Specify a baseline action}"
shift 2
[[ "$CONFIG" == /* && -f "$CONFIG" ]] || { printf 'Use an existing absolute CONFIG path\n' >&2; exit 1; }
case "$ACTION" in
  preflight|prepare|prepare-window|pilot|fit|combine|te-inputs|te-pilot|te-fit|h1-inputs|h1-fit|h1-retry-failed|te-combine|between-genes|context|report-ready|status) ;;
  *) printf 'Unknown baseline action\n' >&2; exit 1 ;;
esac
QUEUE="${NB_QUEUE:-medium}"
MEM="${NB_MEMORY_GB:-16}"
CPUS="${NB_CPUS:-1}"
WALL="${NB_WALL:-12:00}"
MAX="${NB_MAX_CONCURRENT:-10}"
[[ "$MEM" =~ ^[1-9][0-9]*$ && "$CPUS" =~ ^[1-9][0-9]*$ && "$MAX" =~ ^[1-9][0-9]*$ ]] || {
  printf 'Resources and concurrency must be positive integers\n' >&2; exit 1;
}
if [[ "$ACTION" == h1-inputs ]]; then
  [[ -z "${NB_ARRAY_END:-}" ]] || { printf 'h1-inputs requires one job, not an array\n' >&2; exit 1; }
  H1_PREP_WORKERS="${H1_PREP_WORKERS:-$CPUS}"
  H1_PREP_MEMORY_GB="${H1_PREP_MEMORY_GB:-$MEM}"
  [[ "$H1_PREP_WORKERS" =~ ^[1-9][0-9]*$ && "$H1_PREP_MEMORY_GB" =~ ^[0-9]+([.][0-9]+)?$ ]] || {
    printf 'Invalid H1 preparation resource controls\n' >&2; exit 1;
  }
  (( H1_PREP_WORKERS <= CPUS )) || { printf 'H1_PREP_WORKERS exceeds NB_CPUS\n' >&2; exit 1; }
  awk -v budget="$H1_PREP_MEMORY_GB" -v limit="$MEM" 'BEGIN { exit !(budget > 0 && budget <= limit) }' || {
    printf 'H1_PREP_MEMORY_GB must be positive and not exceed NB_MEMORY_GB\n' >&2; exit 1;
  }
fi
LOG_DIR="${NB_LOG_DIR:?Set NB_LOG_DIR to an absolute writable cohort-specific log directory}"
[[ "$LOG_DIR" == /* ]] || { printf 'Use absolute NB_LOG_DIR\n' >&2; exit 1; }
mkdir -p "$LOG_DIR"
script="$(mktemp "$LOG_DIR/nb_${ACTION}.XXXXXX.sh")"
{
  printf '#!/bin/bash\nset -eo pipefail\numask 077\n'
  printf 'export SCRIPTS_DIR=%q\n' "$ROOT"
  if [[ "$ACTION" == h1-inputs ]]; then
    printf 'export H1_PREP_WORKERS=%q H1_PREP_MEMORY_GB=%q\n' "$H1_PREP_WORKERS" "$H1_PREP_MEMORY_GB"
  fi
  printf 'exec %q %q %q %q' "$ROOT/bin/Rscript_sif" "$ROOT/poe_baseline.R" "$CONFIG" "$ACTION"
  for arg in "$@"; do
    if [[ "$arg" == ARRAY_INDEX ]]; then printf ' "${LSB_JOBINDEX:?Missing array index}"'
    else printf ' %q' "$arg"; fi
  done
  printf '\n'
} > "$script"
name="nb_${ACTION}"
placement=()
if [[ "$ACTION" == h1-inputs ]]; then placement=(-R 'span[hosts=1]'); fi
if [[ -n "${NB_ARRAY_END:-}" ]]; then
  [[ "$NB_ARRAY_END" =~ ^[1-9][0-9]*$ ]] || { printf 'Invalid NB_ARRAY_END\n' >&2; exit 1; }
  name="${name}[1-${NB_ARRAY_END}]%${MAX}"
fi
# This site uses GB in -M/-R, matching config/seadragon.sh and cohort YAMLs.
# Override the wrapper/site settings for installations using MB.
"$ROOT/bin/bsub_site" -q "$QUEUE" -J "$name" -n "$CPUS" -M "$MEM" \
  -R "rusage[mem=$MEM]" "${placement[@]}" -W "$WALL" \
  -oo "$LOG_DIR/${ACTION}.%J.%I.out" -eo "$LOG_DIR/${ACTION}.%J.%I.err" /bin/bash "$script"
