#!/bin/bash
# Input test for one cohort, run after retrieval (GUSTO: after private inputs are staged).
# Usage: tests/test_inputs.sh "$CONFIG"
set -eo pipefail
umask 077
ROOT="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CONFIG="${1:?Usage: tests/test_inputs.sh CONFIG}"
[[ "$CONFIG" == /* && -f "$CONFIG" ]] || { printf 'CONFIG must be an existing absolute YAML path\n' >&2; exit 1; }
export SCRIPTS_DIR="$ROOT"
exec "$ROOT/bin/Rscript_sif" "$ROOT/tests/test_inputs.R" "$CONFIG"
