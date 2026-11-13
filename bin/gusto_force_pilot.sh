#!/usr/bin/env bash
# Explicit, auditable --force recovery for GUSTO partial pilot success.
set -euo pipefail
umask 077
if [[ $# -lt 1 ]]; then
  echo 'Usage: bin/gusto_force_pilot.sh CONFIG [--preview | --force | --verify] [--submit] [--max-excluded N]' >&2
  exit 2
fi
CONFIG=$1; shift
MODE=--preview; SUBMIT=0; MAX_EXCLUDED=5
while (($#)); do
  case "$1" in
    --force) MODE=--force; shift ;;
    --preview) MODE=--preview; shift ;;
    --verify) MODE=--verify; shift ;;
    --submit) SUBMIT=1; shift ;;
    --max-excluded) MAX_EXCLUDED=${2:?Missing max-excluded value}; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
if ((SUBMIT)) && [[ $MODE != --force && $MODE != --verify ]]; then
  echo '--submit requires --force or --verify' >&2; exit 2
fi
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
export SCRIPTS_DIR=$ROOT
if [[ -z "${COHORT_DIR:-}" || ! -d "${COHORT_DIR:-/nonexistent}" ]]; then
  echo 'Set COHORT_DIR to the active GUSTO cohort directory' >&2; exit 2
fi
"$ROOT/bin/Rscript_sif" "$ROOT/bin/gusto_force_pilot.R" "$CONFIG" "$MODE" --max-excluded "$MAX_EXCLUDED"
if (( ! SUBMIT )); then
  echo 'No collector submitted (add --submit to --force or --verify).'
  exit 0
fi
# Do NOT alter parent receipt by hand. Run the normal collect-pilot runtime.
# Explicit -W avoids Seadragon's 30-second default walltime.
job_dir="$COHORT_DIR/work/force-pilot-collector-$(date +%Y%m%dT%H%M%S)"
mkdir -p "$job_dir" "$COHORT_DIR/logs/lsf"
cd "$ROOT"
bsub -q short -n 2 -W 02:00 -M 8 -R 'rusage[mem=8]' \
  -J GUSTO_force_collect_pilot \
  -o "$COHORT_DIR/logs/lsf/GUSTO_force_collect_pilot.%J.out" \
  -e "$COHORT_DIR/logs/lsf/GUSTO_force_collect_pilot.%J.err" \
  "$ROOT/bin/Rscript_sif" "$ROOT/runtime.R" "$CONFIG" collect-pilot finish shared "$job_dir"
echo 'Monitor LSF. Submit scale only after status reports pilot complete and QC is checked.'
