#!/usr/bin/env bash
# Stage 11 candidate placental P. prolifica RNA-seq runs using Placenta_EvoDevo's
# EXISTING fetch_sra.sh backend. Does not alter cohort manifests or submit the
# complete pipeline; reciprocal-cross parentage is curated separately.
set -euo pipefail

PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
SCRIPTS_DIR="$PROJECT/scripts"
COHORT_DIR="$PROJECT/data/fish_Pprolifica_AB"
RUNS="$SCRIPTS_DIR/manifests/fish_Pprolifica_AB_placenta_SRA_runs.tsv"

if [[ ${1:-} == worker ]]; then
  : "${LSB_JOBINDEX:?This worker must run inside an LSF array task}"
  [[ "$LSB_JOBINDEX" =~ ^[0-9]+$ ]] || exit 2
  record=$(sed -n "$((LSB_JOBINDEX+1))p" "$RUNS")
  [[ -n $record ]] || { echo "No SRA run for array index $LSB_JOBINDEX" >&2; exit 2; }
  IFS=$'\t' read -r RUN GSM FRAGMENTS BASES SRA_BYTES <<< "$record"
  [[ "$RUN" =~ ^SRR[0-9]+$ && "$FRAGMENTS" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid SRA row: $record" >&2; exit 2; }
  export SCRIPTS_DIR COHORT_DIR
  # Mate output paths and .sra.tsv receipt match the existing pipeline templates.
  exec "$SCRIPTS_DIR/bin/fetch_sra.sh" \
      "$RUN" PE \
      "$COHORT_DIR/raw/${RUN}_1.fastq.gz" \
      "$COHORT_DIR/raw/${RUN}_2.fastq.gz" \
      "$COHORT_DIR/scratch" 2 "$FRAGMENTS" \
      "$COHORT_DIR/raw/${RUN}_1.fastq.gz.sra.tsv" \
      0 0 3 5
fi

[[ ${1:-} == --submit || ${1:-} == --preview || ${1:-} == --verify ]] || {
  echo "Usage: $0 --preview|--submit|--verify" >&2; exit 2;
}
[[ -s "$RUNS" ]] || { echo "Missing $RUNS; copy the provided TSV there" >&2; exit 2; }
[[ -x "$SCRIPTS_DIR/bin/fetch_sra.sh" ]] || {
  echo "Missing executable $SCRIPTS_DIR/bin/fetch_sra.sh" >&2; exit 2;
}
[[ $(($(wc -l < "$RUNS") - 1)) -eq 11 ]] || {
  echo "Expected 11 placental accession rows; refusing submission" >&2; exit 2;
}
mkdir -p "$COHORT_DIR/raw" "$COHORT_DIR/scratch" "$COHORT_DIR/logs/lsf" "$COHORT_DIR/private"

echo "Fish placental SRA accession list:"
column -t -s $'\t' "$RUNS" 2>/dev/null || cat "$RUNS"
echo
printf 'SRA Toolkit backend: %s\nMax concurrent transfers: 2\n' "$SCRIPTS_DIR/bin/fetch_sra.sh"

if [[ $1 == --verify ]]; then
  ok=0
  missing=0
  while IFS=$'\t' read -r run gsm fragments bases sra_bytes; do
    [[ $run == run_accession ]] && continue
    f1="$COHORT_DIR/raw/${run}_1.fastq.gz"
    f2="$COHORT_DIR/raw/${run}_2.fastq.gz"
    receipt="${f1}.sra.tsv"
    if [[ -s $f1 && -s $f2 && -s $receipt ]]; then
      echo "OK $run"
      ((ok+=1))
    else
      echo "INCOMPLETE $run"
      ((missing+=1))
    fi
  done < "$RUNS"
  echo "$ok complete, $missing incomplete (presence-only check; inspect receipts and LSF logs)"
  ((missing==0))
  exit $?
fi

if [[ $1 == --preview ]]; then
  echo 'Preview only. Re-run --submit to schedule the 11 SRA retrieval jobs.'
  exit 0
fi

# Limit to TWO simultaneously running array elements, same backend used by maize/mouse.
# Memory resource conventions match the existing Seadragon scripts (GB).
bsub -q transfer -n 2 -W 08:00 -M 8 -R 'rusage[mem=8]' \
    -J 'fish_sra_placenta[1-11]%2' \
    -o "$COHORT_DIR/logs/lsf/fish_sra_placenta.%J.%I.out" \
    -e "$COHORT_DIR/logs/lsf/fish_sra_placenta.%J.%I.err" \
    bash "$0" worker
