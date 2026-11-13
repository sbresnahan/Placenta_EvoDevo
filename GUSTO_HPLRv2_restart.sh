# GUSTO restart on HPLRv2.0.annotated.PANTRY.gtf
# ---------------------------------------------------------------------------
# Run block by block from a Seadragon login shell. DO NOT run this file
# end-to-end: a few steps wait for LSF jobs, and step 3 is a go/no-go gate
# that must pass before anything is deleted.
#
# Before you start, copy bin/build_transcript_fasta.R into $SCRIPTS_DIR/bin/.
# ---------------------------------------------------------------------------

## 0. Environment ------------------------------------------------------------
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export DATA_ROOT="$PROJECT/data"
export SCRIPTS_DIR="$PROJECT/scripts"
export POE_SITE_FILE="$SCRIPTS_DIR/config/seadragon.sh"
umask 077
cd "$SCRIPTS_DIR"

use_cohort() {
  export COHORT="${1:?Specify cohort}"
  export CONFIG="$SCRIPTS_DIR/config/$COHORT.private.yaml"
  export COHORT_DIR="$DATA_ROOT/$COHORT"
  mkdir -p "$COHORT_DIR"/{manifests,raw,reference,private,work,scratch,logs,receipts}
  for SOURCE in "$SCRIPTS_DIR/manifests/${COHORT}_"*.tsv; do
    DEST="$COHORT_DIR/manifests/$(basename "$SOURCE")"
    [[ -e "$DEST" ]] || cp -- "$SOURCE" "$DEST"
  done
  [[ -e "$CONFIG" ]] || cp -- "$SCRIPTS_DIR/config/$COHORT.example.yaml" "$CONFIG"
}
use_cohort GUSTO

export NEW_GTF=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/reference/HPLRv2.0.annotated.PANTRY.gtf
export NEW_REF_DIR="$DATA_ROOT/reference/HPLRv2.0"          # outside the GUSTO cohort dir
export NEW_TX_FA="$NEW_REF_DIR/HPLRv2.0.annotated.PANTRY.transcripts.fa"
export ARCHIVE="$DATA_ROOT/GUSTO_ESPRESSO_archive_$(date +%Y%m%d)"
export BCFTOOLS="$SCRIPTS_DIR/bin/native/bcftools"
export SAMTOOLS="$SCRIPTS_DIR/bin/native/samtools"
GENOME_SRC=$(awk '/^reference:/{f=1;next} f&&/^  genome_source:/{print $2;exit}' "$CONFIG")
GUSTO_VCF=$(awk '/^gusto_precomputed:/{f=1;next} f&&/^  vcf:/{print $2;exit}' "$CONFIG")
echo "genome_source: $GENOME_SRC"; echo "GUSTO VCF:     $GUSTO_VCF"
ls -l "$NEW_GTF" "$GENOME_SRC" "$GUSTO_VCF"


## 1. Stop all GUSTO jobs (pipeline + regulatory; names start with GUSTO_) ----
bjobs -u "$USER" -J 'GUSTO_*' -noheader -o 'jobid stat queue job_name'
bkill -u "$USER" -J 'GUSTO_*' 0
# Re-run the bjobs line until it prints "No unfinished job found".


## 2. Build the HPLRv2.0 transcript FASTA from the genome ---------------------
mkdir -p "$NEW_REF_DIR"
GENOME_FAI="$GENOME_SRC.fai"
if [[ ! -s "$GENOME_FAI" ]]; then          # source dir may be read-only
  GENOME_FAI="$NEW_REF_DIR/genome_source.fa.fai"
  "$SAMTOOLS" faidx --fai-idx "$GENOME_FAI" "$GENOME_SRC"
fi

bsub -q short -n 2 -M 16 -R "rusage[mem=16]" -W 3:00 -J HPLRv2_txfa \
  -o "$NEW_REF_DIR/build_txfa.%J.out" -e "$NEW_REF_DIR/build_txfa.%J.err" \
  "cd '$SCRIPTS_DIR' && bin/Rscript_sif bin/build_transcript_fasta.R \
     '$GENOME_SRC' '$NEW_GTF' '$NEW_TX_FA' '$GENOME_FAI'"
# WAIT for DONE, then:
"$SAMTOOLS" faidx "$NEW_TX_FA"
cut -f1 "$NEW_TX_FA.fai" | wc -l          # should equal the transcript count in the .out log


## 3. GO/NO-GO: is the precomputed GUSTO VCF in HPLRv2.0 transcript space? ----
# The adapter takes the transcript ID from VCF CHROM and the transcript position
# from POS. The current file was built on ESPRESSO transcripts, so every CHROM
# must exist in the new GTF or child preparation will fail.
# After running GUSTO_vcf_liftover_HPLRv2.sh, check the lifted file instead:
#   export NEW_VCF=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_POE/ornaments/tx_vcf_HPLRv2/biallelic.filtered.txomic.merged.N200.HPLRv2.vcf.gz
#   GUSTO_VCF="$NEW_VCF"
CHK="$NEW_REF_DIR/preflight"; mkdir -p "$CHK"
awk -F'\t' '$3=="exon"{if(match($9,/transcript_id "[^"]+"/))print substr($9,RSTART+15,RLENGTH-16)}' \
  "$NEW_GTF" | sort -u > "$CHK/hplr_tx.txt"
"$BCFTOOLS" query -f '%CHROM\n' "$GUSTO_VCF" | uniq | sort -u > "$CHK/vcf_tx.txt"   # full pass over the VCF; can take a while
cut -f1 "$NEW_TX_FA.fai" | sort -u > "$CHK/fa_tx.txt"

echo "HPLRv2.0 transcripts:            $(wc -l < "$CHK/hplr_tx.txt")"
echo "FASTA vs GTF ID differences:     $(comm -3 "$CHK/fa_tx.txt" "$CHK/hplr_tx.txt" | wc -l)   (must be 0)"
echo "VCF transcripts (CHROM):         $(wc -l < "$CHK/vcf_tx.txt")"
echo "  found in HPLRv2.0 GTF:         $(comm -12 "$CHK/vcf_tx.txt" "$CHK/hplr_tx.txt" | wc -l)"
echo "  MISSING from HPLRv2.0 GTF:     $(comm -23 "$CHK/vcf_tx.txt" "$CHK/hplr_tx.txt" | wc -l)   (must be 0)"
comm -23 "$CHK/vcf_tx.txt" "$CHK/hplr_tx.txt" | head

# If MISSING > 0: STOP. You need a GUSTO transcriptomic VCF projected onto
# HPLRv2.0 transcripts, either regenerated upstream (Placenta_POE tx_vcf) or
# re-projected from the hg38 locus in the VCF ID field. Set it as NEW_VCF and
# use it in step 6. Even when every ID matches, a shared ID with a different
# exon structure in HPLRv2.0 will fail the REF check (bcftools norm -c e
# against the new FASTA) and the genomic-ID cross-check, so a pass here is
# necessary but not sufficient.
# The lifted VCF comes from GUSTO_vcf_liftover_HPLRv2.sh (see the NEW_VCF lines above).


## 4. Archive the small, hand-edited inputs ----------------------------------
mkdir -p "$ARCHIVE"
cp -p -- "$CONFIG" "$ARCHIVE/"
cp -rp -- "$COHORT_DIR/manifests" "$ARCHIVE/manifests"                     # may be force-filtered
cp -rp -- "$COHORT_DIR"/private/force_pilot_* "$COHORT_DIR"/private/force_quant_* "$ARCHIVE/" 2>/dev/null
# The regulatory download cache doesn't depend on the GTF; keep it so nothing is re-downloaded.
[[ -d "$COHORT_DIR/reference/regulatory_resources" ]] && \
  mv -- "$COHORT_DIR/reference/regulatory_resources" "$ARCHIVE/regulatory_resources"
# Check the regulatory block for inputs stored inside $COHORT_DIR (e.g. ENCODE
# cCRE BEDs under private/). Move any GTF-independent ones into $ARCHIVE now.
awk '/^regulatory:/,0' "$CONFIG" | grep -n "$COHORT_DIR" || echo "No regulatory paths inside COHORT_DIR"
ls -la "$ARCHIVE"


## 5. Delete every ESPRESSO-derived GUSTO output -----------------------------
# Keeps only manifests/ (regenerated in step 7) and raw/. Removes reference/,
# private/ (child VCFs, personalized FASTAs/indices, snapshot, job registry,
# force_* audits, qtl_context_v1), work/, scratch/, logs/, receipts/ (including
# locks), additional_v1/, additional_v2/, approved_exports/, and every top-level
# table (exons, tx_map, event_map, events, pairs, pilot_report, model_*, status ...).
: "${COHORT_DIR:?}"; [[ "$COHORT_DIR" == "$DATA_ROOT/GUSTO" ]] || { echo "Wrong COHORT_DIR"; false; }
ls -la "$COHORT_DIR/raw" | head                       # expect empty/symlinks only for GUSTO
find "$COHORT_DIR" -mindepth 1 -maxdepth 1 ! -name manifests ! -name raw -printf '%y  %p\n'   # dry run
# After reviewing the list:
find "$COHORT_DIR" -mindepth 1 -maxdepth 1 ! -name manifests ! -name raw -exec rm -rf -- {} +


## 6. Point the private config at HPLRv2.0 -----------------------------------
sed -i -E \
  -e "s#^(  annotation: ).*#\1$NEW_GTF#" \
  -e "s#^(  transcript_source: ).*#\1$NEW_TX_FA#" \
  -e "s#^(  custom_gtf: ).*#\1$NEW_GTF#" \
  "$CONFIG"
# Point at the lifted VCF (GUSTO_vcf_liftover_HPLRv2.sh):
: "${NEW_VCF:?Set NEW_VCF to the lifted VCF}"
sed -i -E "s#^(  vcf: ).*#\1$NEW_VCF#" "$CONFIG"
grep -nE '^  (annotation|transcript_source|custom_gtf|vcf|genome_source):' "$CONFIG"
diff "$ARCHIVE/$(basename "$CONFIG")" "$CONFIG"


## 7. Recreate directories and regenerate unfiltered manifests ----------------
use_cohort GUSTO
[[ -d "$ARCHIVE/regulatory_resources" ]] && \
  mv -- "$ARCHIVE/regulatory_resources" "$COHORT_DIR/reference/regulatory_resources"

STRAND=$(awk -F'\t' 'NR>1{print $5}' "$ARCHIVE/manifests/GUSTO_libraries.tsv" | sort -u)
STAGE=$(awk -F'\t' 'NR>1{print $5}' "$ARCHIVE/manifests/GUSTO_units.tsv" | sort -u)
echo "strand=$STRAND stage=$STAGE"                   # each must be a single value
GUSTO_VCF=$(awk '/^gusto_precomputed:/{f=1;next} f&&/^  vcf:/{print $2;exit}' "$CONFIG")
python3 bin/prepare_gusto_manifests.py \
  --vcf "$GUSTO_VCF" \
  --origin-matrix /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_POE/ornaments/tx_vcf_SEAD/hap_origin_matrix.tsv \
  --id-map /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/mapqtl/GUSTO/genotypes/trios/infant/infant_samples_n200_fastqIDs.txt \
  --fastq-dir /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/mapqtl/GUSTO/rnaseq/fastq \
  --out-dir "$COHORT_DIR/manifests" --strand-mode "$STRAND" --stage "$STAGE" --overwrite
# Bring back hand edits (shared family_id, covariates). Differences that are only
# rows dropped by the old force-pilot/force-quant step are expected; don't restore those.
for M in trios covariates units libraries sources; do
  echo "== $M"; diff "$ARCHIVE/manifests/GUSTO_$M.tsv" "$COHORT_DIR/manifests/GUSTO_$M.tsv" | head -20
done


## 8. Restart GUSTO ----------------------------------------------------------
bin/submit_pipeline.sh "$CONFIG" reference            # dry run: check the plan uses the new GTF/FASTA
bin/submit_pipeline.sh "$CONFIG" reference --submit   # full build, NOT --resume
# WAIT
bin/Rscript_sif poe.R "$CONFIG" status
column -t "$COHORT_DIR/annotation_audit.tsv"
ls -lh "$COHORT_DIR/receipts/J1001/child.rds" "$COHORT_DIR/receipts/J1001/human.rds"
awk -F '\t' '!/^#/ && length($4)!=length($5){n++} END{print "Personalized VCF indels:",n+0}' \
  "$COHORT_DIR/private/J1010/tx.vcf"

bin/submit_pipeline.sh "$CONFIG" pilot --submit
# WAIT; review assignment_qc.tsv (precomputed_origin / origin_undefined / no_precomputed_variant)
bin/Rscript_sif poe.R "$CONFIG" status

bin/submit_pipeline.sh "$CONFIG" scale --submit
# WAIT
bin/Rscript_sif poe.R "$CONFIG" status
bin/submit_pipeline.sh "$CONFIG" models --submit
# WAIT
bin/submit_pipeline.sh "$CONFIG" additional --submit
# WAIT
bin/Rscript_sif poe.R "$CONFIG" additional-status

# Regulatory (#1/#5). The QTL crosswalk is rebuilt against the new event_map;
# PANTRY HPLRT_*/HPLRG_* IDs should now link directly rather than falling back
# to gene context.
python3 bin/prepare_gusto_qtl_context.py \
  --pantry-dir /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/qtl_results \
  --event-map "$(grep -m1 '^  event_map:' "$CONFIG" | awk '{print $2}')" \
  --out-dir "$COHORT_DIR/private/qtl_context_v1"      # same path the YAML already references
bin/submit_pipeline.sh "$CONFIG" regulatory-resources --submit   # reuses the restored cache
# WAIT
bin/submit_pipeline.sh "$CONFIG" te-tss --submit
# WAIT
bin/submit_pipeline.sh "$CONFIG" te-genes --submit
bin/submit_pipeline.sh "$CONFIG" human-context --submit
bin/Rscript_sif poe_regulatory.R "$CONFIG" regulatory-status
