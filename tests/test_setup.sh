#!/bin/bash
# Setup test: confirm every tool the pipeline calls is reachable on this host.
# Usage: SCRIPTS_DIR=/path/to/scripts tests/test_setup.sh
set -eo pipefail
umask 077
ROOT="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export SCRIPTS_DIR="$ROOT"
source "$ROOT/bin/seadragon_init.sh"

fail=0
ok()  { printf 'ok    %-14s %s\n' "$1" "$2"; }
bad() { printf 'FAIL  %-14s %s\n' "$1" "$2" >&2; fail=1; }

# Host commands used by the submit wrappers and LSF workers.
for tool in bash gzip curl bsub "$POE_SINGULARITY" /usr/bin/time; do
  if path="$(command -v "$tool" 2>/dev/null)"; then ok "$(basename "$tool")" "$path"
  else bad "$(basename "$tool")" "not on PATH"; fi
done

# Native tools, activated exactly as the workers activate them.
for tool in ornaments personalizer samtools bcftools python; do
  if line="$("$ROOT/bin/native_tool.sh" "$tool" --identity 2>&1)"; then
    line="${line##*$'\n'}"   # drop any activation chatter
    IFS=$'\t' read -r _ path version <<< "$line"
    ok "$tool" "$path (${version:0:60})"
  else bad "$tool" "$line"; fi
done

# SRA toolkit (fish and any retrieval.method: sra cohort).
if line="$(source "$POE_MODULE_INIT" && module load sratoolkit 2>&1 &&
            command -v prefetch && command -v fasterq-dump)"; then
  ok sratoolkit "$(printf '%s' "$line" | tail -n 1)"
else bad sratoolkit "module load sratoolkit / prefetch / fasterq-dump failed"; fi

# R container, library and packages.
r_ok=1
if [[ -f "$POE_R_SIF" ]]; then ok r_image "$POE_R_SIF"; else bad r_image "missing: $POE_R_SIF"; r_ok=0; fi
if [[ -d "$POE_R_LIBRARY" ]]; then ok r_library "$POE_R_LIBRARY"; else bad r_library "missing: $POE_R_LIBRARY"; r_ok=0; fi
if (( r_ok )); then
  if out="$("$ROOT/bin/Rscript_sif" -e '
      pk <- c("data.table","yaml","curl","digest","survival","GenomicRanges","IRanges",
              "S4Vectors","Rsamtools","Biostrings")
      miss <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
      if (getRversion() < "4.3.1") stop("R >= 4.3.1 required; found ", getRversion())
      if (length(miss)) stop("missing R packages: ", paste(miss, collapse = ", "))
      cat("R", as.character(getRversion()), "with", length(pk), "packages\n")' 2>&1)"; then
    ok r_packages "$(printf '%s' "$out" | tail -n 1)"
  else bad r_packages "$(printf '%s' "$out" | tail -n 1)"; fi
fi

if (( fail )); then printf '\nSetup test FAILED\n' >&2; exit 1; fi
printf '\nSetup test passed: all required tools are available.\n'
