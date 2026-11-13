#!/bin/bash
set -eo pipefail
umask 077
ROOT="${SCRIPTS_DIR:?Set SCRIPTS_DIR}"
source "$ROOT/bin/seadragon_init.sh"
tool="${1:?Native tool required}"; shift
poe_activate "$tool"
case "$tool" in personalizer) command_name=create_personalized_transcriptome ;; python) command_name=python ;; *) command_name="$tool" ;; esac
exe="$(command -v "$command_name")"
[[ -n "$exe" && -x "$exe" && "$exe" == /* ]] || { printf 'Native executable not resolved: %s\n' "$command_name" >&2; exit 1; }
if [[ "${1:-}" == --identity ]]; then
  exe="$(readlink -f "$exe")"
  case "$tool" in
    personalizer) version="installed console entry point" ;;
    python|samtools|bcftools) version="$("$exe" --version 2>&1)" ;;
    ornaments) version="$("$exe" version 2>&1 || "$exe" --version 2>&1 || true)" ;;
  esac
  version="${version//$'\t'/ }"; version="${version//$'\r'/ }"; version="${version//$'\n'/; }"
  printf '%s\t%s\t%s\n' "$tool" "$exe" "${version:0:500}"
else
  exec "$exe" "$@"
fi
