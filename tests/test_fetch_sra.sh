#!/bin/bash
# Offline SRA transfer-node test: no real NCBI connection or scheduler required.
set -eo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fixture="$(mktemp -d)"
trap 'rm -r -- "$fixture"' EXIT
mkdir -p "$fixture/mockbin" "$fixture/raw with spaces"
export MOCK_PREFETCH_COUNT="$fixture/prefetch.count"
export MOCK_SRA_SCRATCH="$fixture/scratch with spaces"
export POE_MODULE_INIT="$fixture/modules.sh"
printf '%s\n' 'module() { [[ "$1" == load && "$2" == sratoolkit ]]; }' > "$POE_MODULE_INIT"

cat > "$fixture/mockbin/prefetch" <<'SCRIPT'
#!/bin/bash
set -eo pipefail
[[ "$1" == --max-size && "$2" == u && "$3" == --output-directory ]] || exit 71
dir="$4"; run="$5"
[[ "$run" =~ ^(SRR|ERR|DRR)[0-9]+$ ]] || exit 72
mkdir -p "$dir"
printf 'MOCK SRA\n' > "$dir/$run.sra"
printf '%s\n' "$run" >> "$MOCK_PREFETCH_COUNT"
SCRIPT
cat > "$fixture/mockbin/vdb-validate" <<'SCRIPT'
#!/bin/bash
[[ -s "$1" ]]
SCRIPT
cat > "$fixture/mockbin/fasterq-dump" <<'SCRIPT'
#!/bin/bash
set -eo pipefail
[[ "$1" == --split-files && "$2" == --threads ]] || exit 73
threads="$3"; [[ "$4" == --temp && "$6" == --outdir ]] || exit 74
tmp="$5"; dest="$7"; sra="$8"
[[ "$threads" == 2 && "$tmp" == "$MOCK_SRA_SCRATCH"/sra_tmp/* ]] || exit 75
mkdir -p "$dest"
run="$(basename "$sra" .sra)"
if [[ "$run" == ERR992 ]]; then
  printf '@%s\nACGT\n+\n!!!!\n' "$run" > "$dest/$run.fastq"
else
  printf '@%s/1\nACGT\n+\n!!!!\n' "$run" > "$dest/${run}_1.fastq"
  printf '@%s/2\nTGCA\n+\n!!!!\n' "$run" > "$dest/${run}_2.fastq"
fi
SCRIPT
chmod +x "$fixture/mockbin/"*
export PATH="$fixture/mockbin:$PATH"

dest1="$fixture/raw with spaces/ERR991_1.fastq.gz"
dest2="$fixture/raw with spaces/ERR991_2.fastq.gz"
receipt="$dest1.sra.tsv"
args=(ERR991 PE "$dest1" "$dest2" "$MOCK_SRA_SCRATCH" 2 1 "$receipt" \
  123 0123456789abcdef0123456789abcdef 123 0123456789abcdef0123456789abcdef)
bash "$root/bin/fetch_sra.sh" "${args[@]}"
gzip -t "$dest1" "$dest2"
[[ -s "$receipt" && "$(wc -l < "$MOCK_PREFETCH_COUNT")" == 1 ]]
[[ "$(wc -l < "$receipt")" == 3 ]]
bash "$root/bin/fetch_sra.sh" "${args[@]}" # Valid receipt permits reuse.
[[ "$(wc -l < "$MOCK_PREFETCH_COUNT")" == 1 ]]
printf 'tampered' > "$dest1"
if bash "$root/bin/fetch_sra.sh" "${args[@]}" 2>/dev/null; then
  echo 'Tampered SRA destination incorrectly accepted' >&2; exit 1
fi
[[ "$(wc -l < "$MOCK_PREFETCH_COUNT")" == 1 ]]

single="$fixture/raw with spaces/ERR992_1.fastq.gz"
bash "$root/bin/fetch_sra.sh" ERR992 SE "$single" - "$MOCK_SRA_SCRATCH" 2 1 \
  "$single.sra.tsv" 123 0123456789abcdef0123456789abcdef 0 -
gzip -t "$single"
[[ ! -e "$fixture/raw with spaces/ERR992_2.fastq.gz" ]]
echo 'PASS: offline SRA prefetch, scratch isolation, PE/SE splitting, receipt reuse and corruption rejection'
