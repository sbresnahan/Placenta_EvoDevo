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
  printf 'component\tpath\tversion\tsha256\n'
  for tool in ornaments personalizer samtools bcftools python; do
    "$ROOT/bin/native_tool.sh" "$tool" --identity
  done
  revision=unresolved
  if [[ -n "$POE_ORNAMENTS_SOURCE" && -d "$POE_ORNAMENTS_SOURCE/.git" ]]; then
    revision="$(git -C "$POE_ORNAMENTS_SOURCE" rev-parse HEAD)"
    # Dirty tracked sources cannot establish correspondence to a clean revision.
    [[ -z "$(git -C "$POE_ORNAMENTS_SOURCE" status --porcelain --untracked-files=no)" ]] || revision=dirty_source
  fi
  printf 'source_revision\t%s\t%s\t-\n' "${POE_ORNAMENTS_SOURCE:-unresolved}" "$revision"
  hash="$(sha256sum "$POE_SITE_FILE")"; hash="${hash%% *}"
  printf 'site\t%s\tSeadragon site configuration\t%s\n' "$POE_SITE_FILE" "$hash"
  hash="$(sha256sum "$POE_R_SIF")"; hash="${hash%% *}"
  printf 'r_image\t%s\tR image SHA256\t%s\n' "$POE_R_SIF" "$hash"
  freeze="$("$ROOT/bin/native_tool.sh" python -m pip freeze)"
  hash="$(printf '%s' "$freeze" | sha256sum)"; hash="${hash%% *}"
  printf 'python_dependencies\t%s\tpip-freeze digest (not printed)\t%s\n' "$POE_ORNAMENTS_ACTIVATE" "$hash"
} > "$tmp"
mv -f "$tmp" "$dest"
