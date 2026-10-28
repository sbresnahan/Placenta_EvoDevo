#!/bin/bash
set -eo pipefail
umask 077
ROOT="${SCRIPTS_DIR:?Set SCRIPTS_DIR}"
source "$ROOT/bin/seadragon_init.sh"
dest="${1:?Inventory destination required}"
mkdir -p "$(dirname "$dest")"
tmp="$(mktemp "${dest}.XXXXXX")"
trap 'rm -f "$tmp"' EXIT
{
  printf 'component\tpath\tversion\n'
  for tool in ornaments personalizer samtools bcftools python; do
    "$ROOT/bin/native_tool.sh" "$tool" --identity
  done
  printf 'site\t%s\tSeadragon site configuration\n' "$POE_SITE_FILE"
  printf 'r_image\t%s\tConfigured R image\n' "$POE_R_SIF"
  freeze="$("$ROOT/bin/native_tool.sh" python -m pip freeze)"
  freeze="${freeze//$'\t'/ }"; freeze="${freeze//$'\r'/ }"; freeze="${freeze//$'\n'/; }"
  printf 'python_dependencies\t%s\t%s\n' "$POE_ORNAMENTS_ACTIVATE" "$freeze"
} > "$tmp"
mv -f "$tmp" "$dest"
