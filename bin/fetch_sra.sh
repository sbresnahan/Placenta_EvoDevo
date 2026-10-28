#!/bin/bash
# Host-only NCBI SRA retrieval. Never run from the R container.
set -eo pipefail
umask 077
[[ $# == 12 ]] || {
  echo 'Usage: fetch_sra.sh RUN LAYOUT DEST1 DEST2 SCRATCH THREADS READ_COUNT RECEIPT ENA_BYTES1 ENA_MD5_1 ENA_BYTES2 ENA_MD5_2' >&2
  exit 2
}
run="$1"; layout="$2"; dest1="$3"; dest2="$4"; scratch="$5"
threads="$6"; read_count="$7"; receipt="$8"
ena_bytes1="$9"; ena_md51="${10}"; ena_bytes2="${11}"; ena_md52="${12}"
[[ "$run" =~ ^(SRR|ERR|DRR)[0-9]+$ && "$layout" =~ ^(SE|PE)$ &&
   "$threads" =~ ^[1-9][0-9]*$ && "$read_count" =~ ^[1-9][0-9]*$ &&
   "$dest1" == /* && "$scratch" == /* && "$receipt" == /* ]] || {
  echo 'Unsafe SRA argument or accession' >&2; exit 2;
}
[[ "$layout" == PE && "$dest2" == /* || "$layout" == SE && "$dest2" == - ]] || exit 2

verify_exact() {
  local path="$1" bytes="$2" md5="$3" actual
  [[ -f "$path" && "$bytes" =~ ^[1-9][0-9]*$ && "$md5" =~ ^[0-9a-f]{32}$ ]] || return 1
  [[ "$(stat -c %s -- "$path")" == "$bytes" ]] || return 1
  actual="$(md5sum -- "$path")"
  [[ "${actual:0:32}" == "$md5" ]]
}

receipt_valid() {
  [[ -f "$receipt" ]] || return 1
  local mate path values bytes md5 n
  for mate in 1 2; do
    [[ "$mate" == 1 || "$layout" == PE ]] || continue
    path="$dest1"; [[ "$mate" == 1 ]] || path="$dest2"
    values="$(awk -F '\t' -v r="$run" -v m="$mate" -v p="$path" \
      '$1==r && $2==m && $3==p {print $4, $5}' "$receipt")"
    n="$(awk -F '\t' -v r="$run" -v m="$mate" -v p="$path" \
      '$1==r && $2==m && $3==p {n++} END {print n+0}' "$receipt")"
    [[ "$n" == 1 ]] || return 1
    read -r bytes md5 <<< "$values"
    verify_exact "$path" "$bytes" "$md5" || return 1
  done
}

# A prior verified ENA transfer is also acceptable. Never overwrite it.
if verify_exact "$dest1" "$ena_bytes1" "$ena_md51" &&
   { [[ "$layout" == SE ]] || verify_exact "$dest2" "$ena_bytes2" "$ena_md52"; }; then
  echo "Reusing ENA-verified FASTQs: $run"
  exit 0
fi
if receipt_valid; then
  echo "Reusing locally verified SRA-derived FASTQs: $run"
  exit 0
fi
if [[ -e "$dest1" || ( "$layout" == PE && -e "$dest2" ) || -e "$receipt" ]]; then
  echo "Unverified existing SRA destination/receipt; refusing to overwrite $run" >&2
  exit 1
fi

# Site modules are available on LSF transfer hosts, not inside R.
source "${POE_MODULE_INIT:-/etc/profile.d/modules.sh}"
module load sratoolkit
for tool in prefetch fasterq-dump vdb-validate gzip md5sum; do
  command -v "$tool" >/dev/null || { echo "Missing host tool: $tool" >&2; exit 1; }
done

cache="$scratch/sra_cache/$run"
temp_parent="$scratch/sra_tmp"
[[ ! -L "$temp_parent" && ! -L "$cache" ]] || { echo 'SRA scratch symlink rejected' >&2; exit 1; }
mkdir -p "$cache" "$temp_parent" "$(dirname "$dest1")" "$(dirname "$receipt")"
[[ "$layout" == SE ]] || mkdir -p "$(dirname "$dest2")"
temp="$(mktemp -d "$temp_parent/${run}.XXXXXXXX")"
trap 'rm -r -- "$temp"' EXIT
mkdir -p "$temp/output" "$temp/work"

# Preserve the accession as assigned by ENA (ERR), NCBI (SRR), or DDBJ (DRR).
# The NCBI run must actually be available; do not silently substitute another run.
prefetch --max-size u --output-directory "$cache" "$run"
sra="$cache/$run.sra"
if [[ ! -f "$sra" && -f "$cache/$run/$run.sra" ]]; then sra="$cache/$run/$run.sra"; fi
[[ -f "$sra" ]] || { echo "NCBI did not provide $sra" >&2; exit 1; }
vdb-validate "$sra"
fasterq-dump --split-files --threads "$threads" --temp "$temp/work" \
  --outdir "$temp/output" "$sra"

if [[ "$layout" == SE && -s "$temp/output/${run}.fastq" &&
      ! -e "$temp/output/${run}_1.fastq" ]]; then
  mv -- "$temp/output/${run}.fastq" "$temp/output/${run}_1.fastq"
fi
[[ -s "$temp/output/${run}_1.fastq" ]] || { echo 'SRA mate 1 missing' >&2; exit 1; }
if [[ "$layout" == PE ]]; then
  [[ -s "$temp/output/${run}_2.fastq" ]] || { echo 'SRA mate 2 missing' >&2; exit 1; }
else
  [[ ! -e "$temp/output/${run}_2.fastq" ]] || { echo 'Unexpected SRA mate 2' >&2; exit 1; }
fi
[[ ! -e "$temp/output/${run}.fastq" ]] || { echo 'Unpaired SRA reads detected; refusing data loss' >&2; exit 1; }
shopt -s nullglob
produced=("$temp/output/"*.fastq)
expected_count=1; [[ "$layout" == SE ]] || expected_count=2
[[ ${#produced[@]} == "$expected_count" ]] || {
  echo 'Unexpected extra SRA FASTQ output (possible technical read); refusing data loss' >&2
  exit 1
}

gzip -n -c -- "$temp/output/${run}_1.fastq" > "$temp/output/1.fastq.gz"
gzip -t -- "$temp/output/1.fastq.gz"
if [[ "$layout" == PE ]]; then
  gzip -n -c -- "$temp/output/${run}_2.fastq" > "$temp/output/2.fastq.gz"
  gzip -t -- "$temp/output/2.fastq.gz"
fi

# File promotion is fail-closed; an interrupted run without a receipt is not
# silently reused. Keep the downloaded .sra cache for a retry.
[[ ! -e "$dest1" && ( "$layout" == SE || ! -e "$dest2" ) ]] || exit 1
mv -n -- "$temp/output/1.fastq.gz" "$dest1"
[[ ! -e "$temp/output/1.fastq.gz" ]] || { echo 'Destination appeared during SRA promotion' >&2; exit 1; }
if [[ "$layout" == PE ]]; then mv -n -- "$temp/output/2.fastq.gz" "$dest2"; fi
[[ "$layout" == SE || ! -e "$temp/output/2.fastq.gz" ]] || { echo 'Destination appeared during SRA promotion' >&2; exit 1; }
tmp_receipt="$(mktemp "${receipt}.XXXXXXXX")"
{
  printf 'run_id\tmate\tpath\tbytes\tmd5\ttransport\n'
  for mate in 1 2; do
    [[ "$mate" == 1 || "$layout" == PE ]] || continue
    path="$dest1"; [[ "$mate" == 1 ]] || path="$dest2"
    bytes="$(stat -c %s -- "$path")"
    hash="$(md5sum -- "$path")"
    printf '%s\t%s\t%s\t%s\t%s\tsra\n' "$run" "$mate" "$path" "$bytes" "${hash:0:32}"
  done
} > "$tmp_receipt"
mv -- "$tmp_receipt" "$receipt"
echo "SRA-derived FASTQs staged: $run (expected fragments: $read_count; checked by R finish gate)"
