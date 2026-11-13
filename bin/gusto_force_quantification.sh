#!/usr/bin/env bash
# Sample-agnostic, opt-in GUSTO partial quantification recovery.
set -euo pipefail
umask 077
if (( $# < 1 )); then
  echo 'Usage: bin/gusto_force_quantification.sh CONFIG [--preview|--force|--verify] [--expect-excluded N] [--max-excluded N] [--submit]' >&2
  exit 2
fi
CONFIG=$1; shift
MODE=--preview; MODE_SET=0; SUBMIT=0; ARGS=()
while (( $# )); do
  case "$1" in
    --preview|--force|--verify)
      if (( MODE_SET )); then echo 'Choose one mode' >&2; exit 2; fi
      MODE=$1; MODE_SET=1; shift ;;
    --expect-excluded|--max-excluded)
      (( $# >= 2 )) || { echo "Missing value for $1" >&2; exit 2; }
      ARGS+=("$1" "$2"); shift 2 ;;
    --submit) SUBMIT=1; shift ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
if (( SUBMIT )) && [[ $MODE != --force && $MODE != --verify ]]; then
  echo '--submit requires --force or --verify' >&2; exit 2
fi
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
export SCRIPTS_DIR="$ROOT"
: "${COHORT_DIR:?Set COHORT_DIR to data/GUSTO}"
[[ -d "$COHORT_DIR" ]] || { echo "Cohort directory not found: $COHORT_DIR" >&2; exit 2; }
[[ $(basename "$COHORT_DIR") == GUSTO ]] || { echo 'Only GUSTO is supported' >&2; exit 2; }
[[ -f "$CONFIG" ]] || { echo "Config missing: $CONFIG" >&2; exit 2; }

# Don't change the active manifest while existing quantification/assignment
# jobs can still produce outputs. A failed-dependency collector must be canceled.
if [[ $MODE == --force || $SUBMIT == 1 ]]; then
  if command -v bjobs >/dev/null 2>&1; then
    active=$(bjobs -w -u "$USER" 2>/dev/null | awk \
      '$3=="RUN" && $7 ~ /^GUSTO_(quantify|assign|collect-quantification|force_collect_quantification)/ {print}' || true)
    if [[ -n "$active" ]]; then
      echo 'Active GUSTO quantification/assignment/collector jobs still RUN:' >&2
      echo "$active" >&2
      exit 1
    fi
    pending_collector=$(bjobs -w -u "$USER" 2>/dev/null | awk \
      '$3=="PEND" && $7 ~ /^GUSTO_(collect-quantification|force_collect_quantification)/ {print}' || true)
    if [[ -n "$pending_collector" ]]; then
      echo 'Cancel the obsolete pending GUSTO collector before recovery:' >&2
      echo "$pending_collector" >&2
      exit 1
    fi
  fi
fi
"$ROOT/bin/Rscript_sif" "$ROOT/bin/gusto_force_quantification.R" "$CONFIG" "$MODE" "${ARGS[@]}"
if (( ! SUBMIT )); then
  echo 'No collector submitted. Use --force --submit (or --verify --submit after filtering).'
  exit 0
fi
# Ensure the filtered manifest truly contains only complete, valid units.
"$ROOT/bin/Rscript_sif" "$ROOT/bin/gusto_force_quantification.R" "$CONFIG" --verify
mkdir -p "$COHORT_DIR/logs/lsf" "$COHORT_DIR/work"
job_dir="$COHORT_DIR/work/force-quant-collector-$(date +%Y%m%dT%H%M%S)-$$"
mkdir -p "$job_dir"
cd "$ROOT"
bsub -q medium -n 2 -W 04:00 -M 8 -R 'rusage[mem=8]' \
  -J GUSTO_force_collect_quantification \
  -o "$COHORT_DIR/logs/lsf/GUSTO_force_collect_quantification.%J.out" \
  -e "$COHORT_DIR/logs/lsf/GUSTO_force_collect_quantification.%J.err" \
  "$ROOT/bin/Rscript_sif" "$ROOT/runtime.R" "$CONFIG" collect-quantification finish shared "$job_dir"
