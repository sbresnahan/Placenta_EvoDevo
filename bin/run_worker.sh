#!/bin/bash
# Every checked-in LSF stage enters here; native plans execute only on the host.
set -euo pipefail
umask 077
ROOT="${SCRIPTS_DIR:?Set SCRIPTS_DIR explicitly before bsub}"
CONFIG="${CONFIG:-${2:-}}"
[[ "$CONFIG" == /* && -f "$CONFIG" ]] || { printf 'CONFIG must be an existing absolute YAML path\n' >&2; exit 1; }
stage="${1:?Worker stage required}"
sample="${POE_SAMPLE:-${3:-shared}}"
export SCRIPTS_DIR="$ROOT" CONFIG
source "$ROOT/bin/seadragon_init.sh"
rs="$ROOT/bin/Rscript_sif"
bootstrap="$("$rs" "$ROOT/runtime.R" "$CONFIG" "$stage" bootstrap "$sample" unused)"
IFS=$'\t' read -r job lock <<< "$bootstrap"
[[ -d "$job" && -n "$lock" ]] || { printf 'Invalid worker bootstrap\n' >&2; exit 1; }
if ! mkdir "$lock" 2>/dev/null; then
  printf 'Stage/library locked: %s. Inspect owner and active jobs; do not auto-delete.\n' "$lock" >&2; exit 1
fi
printf 'pid=%s job=%s host=%s\n' "$$" "${LSB_JOBID:-local}" "$(hostname)" > "$lock/owner"
cleanup() { rm -f "$lock/owner"; rmdir "$lock"; }
trap cleanup EXIT
export POE_NATIVE_RECORD="$job/native.tsv"
"$ROOT/bin/host_inventory.sh" "$POE_NATIVE_RECORD"
run_r() {
  local entry=""
  case "$stage" in
    reference) entry=02_reference/prepare_reference.R ;;
    markers) entry=03_genotypes/prepare_markers.R ;;
    child) entry=03_genotypes/prepare_child.R ;;
    assign) entry=05_parental_expression/assign_parents.R ;;
    models) entry=06_promoter_models/run_models.R ;;
    retrieve) entry=01_retrieval/retrieve_public.R ;;
  esac
  if [[ -n "$entry" ]]; then "$rs" "$ROOT/$entry" "$CONFIG" "$1" "$sample" "$job"
  else "$rs" "$ROOT/runtime.R" "$CONFIG" "$stage" "$1" "$sample" "$job"; fi
}
case "$stage" in
  reference)
    run_r prepare
    bash "$job/pre.sh"
    run_r middle
    bash "$job/post.sh"
    run_r finish ;;
  child|index|pilot|quantify|retrieve)
    run_r prepare
    bash "$job/pre.sh"
    run_r finish ;;
  markers|collect-pilot|assign|collect-quantification|models|preflight)
    run_r finish ;;
  *) printf 'Unknown worker stage: %s\n' "$stage" >&2; exit 1 ;;
esac
