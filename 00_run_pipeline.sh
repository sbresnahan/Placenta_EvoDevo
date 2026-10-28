#!/bin/bash
#BSUB -J poe_submit
#BSUB -q short
#BSUB -n 1
#BSUB -M 4
#BSUB -R "rusage[mem=4]"
#BSUB -W 0:30
#BSUB -o /rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/data/logs/submit.%J.out
#BSUB -e /rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/data/logs/submit.%J.err
set -euo pipefail
umask 077
ROOT="${SCRIPTS_DIR:?Set SCRIPTS_DIR explicitly}"
CONFIG="${CONFIG:?Set CONFIG to an absolute cohort YAML}"
STAGE="${STAGE:?Set STAGE to reference, pilot, scale, models, or retrieve}"
args=("$CONFIG" "$STAGE")
[[ "${SUBMIT:-no}" != yes ]] || args+=(--submit)
[[ -z "${ACCEPT:-}" ]] || args+=(--accept "$ACCEPT")
[[ -z "${AFTER:-}" ]] || args+=(--after "$AFTER")
exec "$ROOT/bin/submit_pipeline.sh" "${args[@]}"
