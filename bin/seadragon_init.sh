#!/bin/bash
# Source in a host subprocess, never inside the R container.
POE_ROOT="${SCRIPTS_DIR:?Set SCRIPTS_DIR to the scripts directory}"
POE_SITE_FILE="${POE_SITE_FILE:-${POE_ROOT}/config/seadragon.sh}"
export POE_SITE_FILE SCRIPTS_DIR
[[ -f "$POE_SITE_FILE" ]] || { printf 'Missing site file: %s\n' "$POE_SITE_FILE" >&2; return 1; }
source "$POE_SITE_FILE"
poe_activate() {
  case "$1" in
    samtools|bcftools)
      [[ -f "$POE_MODULE_INIT" && -x "$POE_CONDA" ]] || { printf 'Modules/conda initialization missing\n' >&2; return 1; }
      source "$POE_MODULE_INIT"
      eval "$("$POE_CONDA" shell.bash hook)"
      if [[ "$1" == samtools ]]; then conda activate "$POE_SAMTOOLS_ENV"
      else conda activate "$POE_BCFTOOLS_ENV"; fi
      # Native htslib tools may need this on Seadragon; only in this subprocess.
      if [[ -n "${CONDA_PREFIX:-}" ]]; then export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"; fi
      ;;
    ornaments|personalizer|python)
      [[ -f "$POE_ORNAMENTS_ACTIVATE" ]] || { printf 'Ornaments activation missing: %s\n' "$POE_ORNAMENTS_ACTIVATE" >&2; return 1; }
      source "$POE_ORNAMENTS_ACTIVATE"
      ;;
    *) printf 'Unknown native environment: %s\n' "$1" >&2; return 1 ;;
  esac
}
