# Lift biallelic.filtered.txomic.merged.N200.fixed.vcf.gz from ESPRESSO to HPLRv2.0
# ---------------------------------------------------------------------------
# Run block by block on Seadragon. Prerequisite: the HPLRv2.0 transcript FASTA
# from GUSTO_HPLRv2_restart.sh step 2 ($NEW_TX_FA plus its .fai).
# Copy bin/liftover_gusto_vcf.R into $SCRIPTS_DIR/bin/ first.
#
# Linking: gffcompare v0.12.10 (local build) -r HPLRv2.0 (reference) vs ESPRESSO (query) gives each
# ESPRESSO transcript its best HPLRv2.0 match (.tmap ref_id + class code).
# Shared ENST IDs are ignored; only structure counts. Variants are re-projected
# from the hg38 locus in VCF ID onto the matched transcript and REF-checked
# against the HPLRv2.0 transcript FASTA. GT/DS/GP columns are copied unchanged,
# and hap_origin_matrix.tsv is unaffected.
# ---------------------------------------------------------------------------

## 0. Environment ------------------------------------------------------------
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export DATA_ROOT="$PROJECT/data"
export SCRIPTS_DIR="$PROJECT/scripts"
export POE_SITE_FILE="$SCRIPTS_DIR/config/seadragon.sh"
umask 077
cd "$SCRIPTS_DIR"
source "$POE_SITE_FILE"

BCF="$SCRIPTS_DIR/bin/native/bcftools"
SRC_VCF=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_POE/ornaments/tx_vcf_SEAD/biallelic.filtered.txomic.merged.N200.fixed.vcf.gz
ESP_GTF=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered/ESPRESSO_corrected_SC_filtered.gtf
HPLR_GTF=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/reference/HPLRv2.0.annotated.PANTRY.gtf
NEW_TX_FA="$DATA_ROOT/reference/HPLRv2.0/HPLRv2.0.annotated.PANTRY.transcripts.fa"
OUT_DIR=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_POE/ornaments/tx_vcf_HPLRv2   # private, beside the original
OUT_PREFIX="$OUT_DIR/biallelic.filtered.txomic.merged.N200.HPLRv2"
FINAL_VCF="$OUT_PREFIX.vcf.gz"
GFFCOMPARE=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/gffcompare_src/gffcompare/gffcompare   # v0.12.10
mkdir -p "$OUT_DIR"
ls -l "$SRC_VCF" "$ESP_GTF" "$HPLR_GTF" "$NEW_TX_FA" "$NEW_TX_FA.fai"


## 1. Check gffcompare ------------------------------------------------------------
"$GFFCOMPARE" --version


## 2. Best-match ESPRESSO -> HPLRv2.0 with gffcompare ---------------------------
# gffcompare writes .tmap/.refmap beside the QUERY file, so work on a copy.
cp -- "$ESP_GTF" "$OUT_DIR/ESPRESSO_corrected_SC_filtered.gtf"
bsub -q short -n 1 -M 8 -R "rusage[mem=8]" -W 1:00 -J HPLRv2_gffcompare \
  -o "$OUT_DIR/gffcompare.%J.out" -e "$OUT_DIR/gffcompare.%J.err" \
  "cd '$OUT_DIR' && '$GFFCOMPARE' -r '$HPLR_GTF' -o esp_vs_hplr ESPRESSO_corrected_SC_filtered.gtf"
# WAIT for DONE
TMAP=$(ls "$OUT_DIR"/esp_vs_hplr*.tmap); echo "$TMAP"
# Class-code mix for ESPRESSO transcripts (= identical intron chain; c contained;
# k containment of ref; m/n retained intron(s); j novel junction; u no match ...)
awk -F'\t' 'NR>1{n[$3]++} END{for(k in n) print n[k]"\t"k}' "$TMAP" | sort -rn
# The same counts restricted to transcripts that actually carry VCF records:
"$BCF" query -f '%CHROM\n' "$SRC_VCF" | sort -u > "$OUT_DIR/vcf_tx.txt"
awk -F'\t' 'NR==FNR{keep[$1];next} FNR>1 && ($5 in keep){n[$3]++} END{for(k in n) print n[k]"\t"k}' \
  "$OUT_DIR/vcf_tx.txt" "$TMAP" | sort -rn
# Default accepted codes: =,c,k,m,n,j. Change --class-codes in step 3 if this mix argues for it.


## 3. Lift the VCF -------------------------------------------------------------
"$BCF" view -Ov -o "$OUT_PREFIX.source.vcf" "$SRC_VCF"       # uncompressed working copy (private)
bsub -q short -n 2 -M 16 -R "rusage[mem=16]" -W 3:00 -J HPLRv2_liftover \
  -o "$OUT_DIR/liftover.%J.out" -e "$OUT_DIR/liftover.%J.err" \
  "cd '$SCRIPTS_DIR' && bin/Rscript_sif bin/liftover_gusto_vcf.R \
     --vcf '$OUT_PREFIX.source.vcf' --tmap '$TMAP' \
     --target-gtf '$HPLR_GTF' --target-fasta '$NEW_TX_FA' \
     --out-prefix '$OUT_PREFIX' --class-codes '=,c,k,m,n,j'"
#   Add --project-all to also personalize HPLRv2.0 isoforms with no ESPRESSO best match.
# WAIT for DONE, then review:
column -t -s $'\t' "$OUT_PREFIX.summary.tsv"
column -t -s $'\t' "$OUT_PREFIX.tmap_class_codes.tsv" | head -20


## 4. Sort, left-align in transcript space, REF-check, compress, index ----------
"$BCF" sort -T "$OUT_DIR/sort_tmp" -Ov -o "$OUT_PREFIX.sorted.vcf" "$OUT_PREFIX.unsorted.vcf"
"$BCF" norm -f "$NEW_TX_FA" -c e -Ov -o "$OUT_PREFIX.norm.vcf" "$OUT_PREFIX.sorted.vcf"
"$BCF" sort -T "$OUT_DIR/sort_tmp" -Oz -o "$FINAL_VCF" "$OUT_PREFIX.norm.vcf"
"$BCF" index -t "$FINAL_VCF"


## 5. Verify ----------------------------------------------------------------------
"$BCF" view -H "$FINAL_VCF" | wc -l                               # = lifted_records in summary
"$BCF" view -H "$FINAL_VCF" | head -10 | cut -f1-9
"$BCF" query -l "$FINAL_VCF" | cmp - <("$BCF" query -l "$SRC_VCF") && echo "Sample columns identical"
# Every CHROM is an HPLRv2.0 transcript (the restart runbook step-3 gate):
"$BCF" query -f '%CHROM\n' "$FINAL_VCF" | sort -u > "$OUT_DIR/lifted_tx.txt"
cut -f1 "$NEW_TX_FA.fai" | sort -u > "$OUT_DIR/hplr_tx.txt"
echo "Lifted CHROMs missing from HPLRv2.0: $(comm -23 "$OUT_DIR/lifted_tx.txt" "$OUT_DIR/hplr_tx.txt" | wc -l)   (must be 0)"
# No repeated transcript position/allele records:
"$BCF" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$FINAL_VCF" | sort | uniq -d | wc -l   # expect 0
# IDs still satisfy the GUSTO adapter's locus pattern:
"$BCF" query -f '%ID\n' "$FINAL_VCF" | grep -vcE '^chr([0-9]+|X|Y|M|MT):[0-9]+:[ACGT]+:[ACGT]+v[0-9]+$'   # expect 0
# Where were records dropped?
awk -F'\t' 'NR>1{n[$9]++} END{for(k in n) print n[k]"\t"k}' "$OUT_PREFIX.source_records.tsv" | sort -rn
# Spot-check the loci shown earlier (minus-strand indel + multi-transcript SNPs):
"$BCF" view -H "$FINAL_VCF" | awk -F'\t' '$3 ~ /^chr10:(1014278|1015827|1040857):/' | cut -f1-5


## 6. Remove uncompressed genotype intermediates ----------------------------------
rm -f -- "$OUT_PREFIX".source.vcf "$OUT_PREFIX".unsorted.vcf "$OUT_PREFIX".sorted.vcf "$OUT_PREFIX".norm.vcf
rm -rf -- "$OUT_DIR/sort_tmp"


## 7. Hand off to the restart ------------------------------------------------------
# In GUSTO_HPLRv2_restart.sh use:
#   export NEW_VCF="$FINAL_VCF"
# then step 3 (gate) should report 0 missing, and step 6 sets gusto_precomputed.vcf.
echo "$FINAL_VCF"
