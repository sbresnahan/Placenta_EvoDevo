#!/bin/bash
set -euo pipefail
# xQTL-style loop interface: driver submits one checked-in worker per library.
ROOT="${SCRIPTS_DIR:?Set SCRIPTS_DIR}"
CONFIG="${1:?Absolute CONFIG required}"; STAGE="${2:-pilot}"; shift; [[ $# -eq 0 ]] || shift
[[ "$STAGE" == pilot || "$STAGE" == scale ]] || { printf "Use pilot or scale\n" >&2; exit 1; }
exec "$ROOT/bin/submit_pipeline.sh" "$CONFIG" "$STAGE" "$@"
