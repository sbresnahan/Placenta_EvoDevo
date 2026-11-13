# Placenta_EvoDevo — revised analysis runbook (additive v1)

Applies to repository revision **after applying the supplied patch**, on MD Anderson Seadragon. The original reference / pilot / scale / models stages and their prior outputs are preserved. New work is in the **`additional`** stage (`data/<cohort>/additional_v1/`) and the aggregate-only cross-cohort comparator.

## 0. What is being tested

| Question | Cohorts | Endpoint | Meaning and limitations |
|---|---|---|---|
| Q1. Is dosage-adjusted allele-specific expression systematically maternal or paternal in reciprocal crosses? | Maize, mouse | Baseline `model_gene.tsv`; exploratory `additional_v1/gene_POE.tsv` | Reciprocal-cross intercept, conditional on genotype markers and Ornaments' estimated counts. In maize, subtract expected 2:1 maternal:paternal dosage. |
| Q2. Does the *TSS choice* change parental imbalance, conditional on an identical first splice donor? | Maize, mouse | `additional_v1/separated_TSS_donor_POE.tsv` with `TSS_conditional_on_first_donor` | Distinct TSS regions in the same gene and same first donor; not proof of causal promoter regulation or a molecule-spanning read. |
| Q3. Does the *first-donor choice* change parental imbalance, conditional on nearby TSS? | Maize, mouse | Same file with `first_donor_conditional_on_TSS` | Different first donor within 50 bp TSS matching; not a purely isolated splicing perturbation. |
| Q4. Is there gene-level parent-of-origin imbalance beyond cis orientation? | GUSTO | New `additional_v1/gene_POE.tsv` | Trio-specific chromosome parental assignment; adjusts for cis allele and configured covariates; one child/family and one fixed SNP/gene. Private analysis. |
| Q5. Do two transcript events differ in their parental expression? | GUSTO | Frozen baseline `model_promoter.tsv` | HC3/within-family bootstrap model conditional on cis orientation. GUSTO event pairs still follow the old structural definition; no new GUSTO two-event refit is implied. |
| Q6. Do TSS-overlapping TEs associate with altered transcript parental imbalance? | Cohorts with matching repeat annotation | New `event_TE_annotation.tsv`; the between-gene TE contrast is now `te-genes` (`additional_v2/te_genes_result.tsv`, see `REGULATORY_MODULES.md` #1b) | Annotation and matched contrasts run only with a verified build-matched UCSC RepeatMasker table; otherwise explicitly `not_assessed`. Exploratory, not causal. |
| Q7. Are gene-level effects shared across cohorts/species? | All with comparable aggregates | `within_stratum_meta.tsv`, `directional_conservation.tsv` | Fixed-effect meta **only within same species × tissue × developmental stage and independent biological samples**. Cross-species direction comparison is descriptive, with curated one-to-one orthologs. |

The additive stage is **more permissive** than the original two-event/full-callability gate: default ≥10 expected counts per event/locus per unit, ≥2 independent reciprocal units **in each direction**, same-locus comparison for structural contrasts, up to 50 structurally supported candidate pairs per gene/type screened for callability and the top 10 selected by callable unit count (never effect size). The result records `direct_coordinate_both`, `nominal_p_lt_0_05`, BH `q_value`, 90% CI equivalence to a 1.25-fold margin, and explicit ineligibility states. This **does not** establish read-level molecule phasing, independent transcription starts, or statistical replication with small n. `p < 0.05` is nominal, not FDR control; lack of significance is not proof of equivalence.

## 1. One-time environment and cohort selector

Run at a Seadragon **login** shell:

```bash
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
```

Make sure the updated code has reached the **same** `$SCRIPTS_DIR` used for prior runs. Private YAML and manifests are *not* replaced by code updates. To inspect a stage:

```bash
use_cohort maize_B73_Ki11  # change cohort as needed
bin/Rscript_sif poe.R "$CONFIG" status
bin/Rscript_sif poe.R "$CONFIG" additional-status
bjobs -a
ls -lt "$COHORT_DIR/logs/lsf/" | head
```

Submitting `bin/submit_pipeline.sh "$CONFIG" <stage> --submit` creates and submits LSF jobs. Its return means **submitted**, not completed. Wait for `DONE` and the pipeline's receipt/status before starting the dependent stage. A dry run omits `--submit`. Do not rerun `retrieve`, `reference`, `pilot`, `scale`, or `models` merely because this patch was applied.

### One-time optional analysis configuration

Defaults are automatically applied by `additional`; existing private YAMLs need not be edited. To specify different settings, add this *top-level* section to a cohort private YAML **before its first `additional` run**:

```yaml
additional_analysis:
  min_depth: 10
  min_per_direction: 2
  max_pairs_per_gene: 10
  # te_rmsk: /absolute/path/assembly-matched-ucsc-rmsk.txt.gz
```

`te_rmsk` must be the 17-column UCSC RepeatMasker `rmsk` table, gzipped; verify it uses the **exact reference genome build and matching chromosome naming**. Do not add this key until verified. Without it the TE outcome is explicitly not tested. For a changed exploratory parameter set, create a new versioned analysis implementation/output path rather than overwriting `additional_v1`.

## 2. Maize B73 × Ki11 — 20 DAP endosperm

This cohort's original five-stage run has already reported **`reference complete`, `pilot complete`, `scale requested`, `quantification complete`, `models complete`**. Keep those preliminary results.

```bash
use_cohort maize_B73_Ki11
bin/Rscript_sif poe.R "$CONFIG" status

# New analyses ONLY; reads the existing complete parent-assigned counts.
bin/submit_pipeline.sh "$CONFIG" additional          # optional dry run
bin/submit_pipeline.sh "$CONFIG" additional --submit
```

After LSF `DONE`:

```bash
bin/Rscript_sif poe.R "$CONFIG" additional-status
ls -lh "$COHORT_DIR/additional_v1/"
less -S "$COHORT_DIR/additional_v1/gene_POE.tsv"
less -S "$COHORT_DIR/additional_v1/separated_TSS_donor_POE.tsv"
```

Expect **`additional_v1 complete`**. Inspect full status and direct-coordinate support, and report nominal p < 0.05 separately from `q_value <= 0.05`. No 29-file retrieval recheck and no new HapMap3 download is needed for this step.

## 3. Maize B73 × W22 — 14 DAP endosperm

Last verified status was **`reference complete`**, with **pilot blocked** because it had not run. Confirm current status first—work may have progressed. Existing reference reuses the *Ki11* HapMap3 AGPv3 VCFs, not Ki11-specific B73/Ki11 marker calls. The working W22 marker manifest must have **11 lines** and columns `source_id path build processed_path`.

```bash
use_cohort maize_B73_W22
bin/Rscript_sif poe.R "$CONFIG" status
wc -l "$COHORT_DIR/manifests/maize_B73_W22_marker_sources.tsv"
```

If and only if `pilot` is incomplete and `reference` complete:

```bash
bin/submit_pipeline.sh "$CONFIG" pilot --submit
# WAIT for pilot DONE, then:
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/pilot_report.tsv"
```

If and only if `quantification` is incomplete after pilot completes:

```bash
bin/submit_pipeline.sh "$CONFIG" scale --submit
# WAIT for quantify + assign + collector jobs DONE, then:
bin/Rscript_sif poe.R "$CONFIG" status
```

If `models` is incomplete after quantification completes:

```bash
bin/submit_pipeline.sh "$CONFIG" models --submit
# WAIT for models DONE, then:
bin/Rscript_sif poe.R "$CONFIG" status
```

Finally, once `models complete`:

```bash
bin/submit_pipeline.sh "$CONFIG" additional --submit
# WAIT for additional DONE, then:
bin/Rscript_sif poe.R "$CONFIG" additional-status
less -S "$COHORT_DIR/additional_v1/separated_TSS_donor_POE.tsv"
```

Do not mix W22 14 DAP and Ki11 20 DAP estimates in an unstratified meta-analysis. The new comparator explicitly separates these stages.

## 4. Mouse Babr × CAST — E7.5 extraembryonic ectoderm

```bash
use_cohort mouse_Babr_CAST
bin/Rscript_sif poe.R "$CONFIG" status
```

For a new cohort, review `$COHORT_DIR/manifests/` sources, libraries, units, resources, and marker sources. Confirm `cluster_id` distinguishes biological units; inspect maternal/paternal assignment and the configured single-end fragment mean/SD. This is **extraembryonic ectoderm**, not placental tissue. The MGP `REFERENCE_C57BL_6J` genotype is a proxy for Babr and must remain labeled as such.

If raw data/resources are not yet fetched, submit retrieval via **transfer**:

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
# WAIT for transfer jobs and resource receipt
less -S "$COHORT_DIR/resource_receipt.tsv"
```

Supply/verify the exact genome-matched source GTF and parental-marker files specified in `mouse_Babr_CAST_marker_sources.tsv`; inspect detailed acquisition instructions in `docs/Placenta_EvoDevo_step_by_step.md`, Part 4. Then advance **one completed stage at a time**:

```bash
bin/submit_pipeline.sh "$CONFIG" reference --submit
# WAIT; verify `reference complete`
bin/Rscript_sif poe.R "$CONFIG" status

bin/submit_pipeline.sh "$CONFIG" pilot --submit
# WAIT; verify `pilot complete`
bin/Rscript_sif poe.R "$CONFIG" status

bin/submit_pipeline.sh "$CONFIG" scale --submit
# WAIT; verify `quantification complete`
bin/Rscript_sif poe.R "$CONFIG" status

bin/submit_pipeline.sh "$CONFIG" models --submit
# WAIT; verify `models complete`
bin/Rscript_sif poe.R "$CONFIG" status

bin/submit_pipeline.sh "$CONFIG" additional --submit
# WAIT; verify
bin/Rscript_sif poe.R "$CONFIG" additional-status
```

The E7.5 single-end fragment settings are assumptions unless externally measured; evaluate that uncertainty when interpreting transcript events. `marker_proxy` is recorded in the results.

## 5. Mouse FVB × CAST — E12.5 placenta

```bash
use_cohort mouse_FVB_CAST_E12
bin/Rscript_sif poe.R "$CONFIG" status
```

Review the paired-end sources, unit clusters, reciprocal orientations, GTF and parental-marker manifest. Confirm the **FVB** and **CAST** parental marker identities and build (not the Babr genotype proxy). If new, retrieve into the cohort via transfer:

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
# WAIT; inspect resource_receipt.tsv; validate reference/marker sources
```

Then, advancing only after preceding stages complete:

```bash
bin/submit_pipeline.sh "$CONFIG" reference --submit
# WAIT; bin/Rscript_sif poe.R "$CONFIG" status
bin/submit_pipeline.sh "$CONFIG" pilot --submit
# WAIT; bin/Rscript_sif poe.R "$CONFIG" status
bin/submit_pipeline.sh "$CONFIG" scale --submit
# WAIT; bin/Rscript_sif poe.R "$CONFIG" status
bin/submit_pipeline.sh "$CONFIG" models --submit
# WAIT; bin/Rscript_sif poe.R "$CONFIG" status
bin/submit_pipeline.sh "$CONFIG" additional --submit
# WAIT; bin/Rscript_sif poe.R "$CONFIG" additional-status
```

Use `docs/Placenta_EvoDevo_step_by_step.md`, Part 5, for source files, manifest names and parent-marker quality gates. E12.5 placenta and E7.5 extraembryonic ectoderm are distinct strata even though both are mouse.

## 6. Human GUSTO — precomputed phased transcript variants and chromosome parental origin (private)

**No existing GUSTO preliminary outputs need preserving.** The new GUSTO adapter reuses completed infant phasing/imputation and parental-origin assignment. The gene-level analysis remains additive, after the primary models. Other cohorts and the legacy trio adapter are unchanged.

```bash
use_cohort GUSTO
```

1. Copy the updated `config/GUSTO.example.yaml` into your private GUSTO configuration (preserve any local site settings). The active input fields are:

```yaml
phase_origin:
  enabled: false
gusto_precomputed:
  enabled: true
  vcf: /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_POE/ornaments/tx_vcf_SEAD/biallelic.filtered.txomic.merged.N200.fixed.vcf.gz
  origin_matrix: /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_POE/ornaments/tx_vcf_SEAD/hap_origin_matrix.tsv
  fastq_id_map: /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/mapqtl/GUSTO/genotypes/trios/infant/infant_samples_n200_fastqIDs.txt
reference:
  build: GRCh38
  contig_aliases: []
  genome_source: /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon.fasta
  transcript_source: /rsrch5/home/epi/bhattacharya_lab/data/Placenta_scRNAseq/salmon/ESPRESSO_corrected_SC_filtered.fasta
```

Also set `paths.annotation` to:

```text
/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered/ESPRESSO_corrected_SC_filtered.gtf
```

`paths.typed_vcf` and `paths.personalization_vcf` are legacy schema fields and are ignored by this adapter. The uploaded **sample** VCF is for local format tests; keep the **full** `.fixed.vcf.gz` as the live input.

2. Create the private GUSTO manifests from the actual full VCF and RNA identity table. Specify the **experimentally verified** library strand mode (`fr`, `rf` or `unstranded`); example:

```bash
python3 bin/prepare_gusto_manifests.py \
  --vcf /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_POE/ornaments/tx_vcf_SEAD/biallelic.filtered.txomic.merged.N200.fixed.vcf.gz \
  --origin-matrix /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_POE/ornaments/tx_vcf_SEAD/hap_origin_matrix.tsv \
  --id-map /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/mapqtl/GUSTO/genotypes/trios/infant/infant_samples_n200_fastqIDs.txt \
  --fastq-dir /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/mapqtl/GUSTO/rnaseq/fastq \
  --out-dir "$COHORT_DIR/manifests" --strand-mode unstranded
```

This command refuses to replace existing manifests; after reviewing existing files, add `--overwrite` **only if replacement is intended**. It creates source, library, unit, child/family and covariate manifests for the exact VCF/matrix/FASTQ overlap. The underscore in a VCF column such as `010-04002_B010-04002` separates an upstream label from the actual infant genotype identity, **`B010-04002`**. The adapter selects the literal composite column with bcftools but uses only the `B...` identity for biological linkage. No fabricated parental IDs are required. By default, `family_id` is the infant's unique ID; audit family independence and update shared family IDs before modeling. Confirm cohort stage metadata. The helper deliberately leaves total FASTQ `read_count` blank rather than scan hundreds of private compressed reads on a login node; GUSTO pilots then report projected full-run time as `NA`, while final quantification counts/validates all input reads. Other cohort pilot gates remain strict.

The uploaded VCF sample contains 197 genotypes. Three mapped infants (`B010-20753`, `B010-21035`, `B020-66040`) are absent from that sample's VCF header. The helper checks the **full** VCF header at runtime and adapts the manifest accordingly.

3. GUSTO's additional gene bootstrap can exceed the original short model stage. Confirm this optional resource block **under `lsf.stages`**:

```yaml
    additional:
      memory_value: 16.0
      reservation_value: 16.0
      cpus: 2
      wall_time: '12:00'
      queue: medium
```

The source VCF uses transcript identifiers in `CHROM` and transcript-relative `POS`. Its `ID` stores an hg38 locus, e.g. `chr10:1014278:G:Av1`; the `vN` suffix is stripped to link the same genomic variant across transcripts. **Both biallelic SNPs and indels are retained for transcriptome personalization and indexing**, with the transcript FASTA used for native REF validation and indel normalization. The phased fetal `GT` first allele is Left, second Right; `hap_origin_matrix.tsv` has one chromosome row per sample, not SNP rows. For downstream parental-count orientation, only SNPs have a direct genomic-coordinate mapping; indels remain in the personalized reference but are excluded from parent-of-origin locus tests until allele-aware indel coordinate mapping is implemented. Only phased fetal heterozygous SNPs with a defined maternal Left/Right chromosome assignment are orientable. Undefined chromosomes and homozygous genotypes retain missing maternal/paternal counts. Genomic position and allele orientation must match the GTF exon mapping or SNP assignment fails; the source FASTA's transcript IDs and spliced lengths are compared against the GTF at reference completion. This adapter **does not infer or repair** phasing/phase-set switches; that accuracy belongs to the supplied parent-of-origin pipeline.

4. Inspect real VCF sample IDs; the configured native bcftools wrapper is activated on the host:

```bash
# Inspect VCF sample IDs and ensure every trio ID is present.
# Example using the available site bcftools wrapper, as configured:
"$SCRIPTS_DIR/bin/native/bcftools" query -l /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_POE/ornaments/tx_vcf_SEAD/biallelic.filtered.txomic.merged.N200.fixed.vcf.gz | head
```

5. Build the ESPRESSO reference and per-child indices:

```bash
bin/submit_pipeline.sh "$CONFIG" reference --submit
# WAIT for reference DONE, verify
bin/Rscript_sif poe.R "$CONFIG" status
```

**GUSTO transcript VCF boundary records and recovery.** The original private
VCF can include nonstandard clipped deletions with nonpositive transcript
`POS`, `REF=.` and a numeric length in `ALT` (91 such records were found in a
2026-11-11 audit). These are not sequence-resolved VCF alleles. The GUSTO
adapter excludes nonpositive positions and missing REF/ALT values with
bcftools comparisons verified on the full GUSTO file, while retaining ordinary
sequence-resolved SNPs **and indels**. The separate audit reports all
non-nucleotide alleles (including numeric deletion placeholders). bcftools
1.16 on Seadragon returned zero matches for regex checks on REF/ALT, so those
predicates must not be used for this source. Retained records still undergo strict FASTA REF validation;
never replace `-c e` with an automatic reference-allele rewrite.

**Per-child genotype filtering:** bcftools evaluates `-i 'GT="alt"'` before
`-s SAMPLE` if they appear in the same `view` call. A first `view -s SAMPLE`
must therefore create the one-infant subset; a second `view -i` retains only
that infant's alternate-bearing calls. A historical J1001 query incorrectly
retained 94,739 `0|0` and 32,370 `./.` genotypes because the filter considered
other infants. A successful child query must contain only `0|1`, `1|0`,
`1|1` (or unphased equivalent) genotypes. The pipeline now enforces this
ordering and removes the temporary one-infant subset after normalization.

**Ambiguous multiallelic positions:** a small number of genomic sites have
multiple alternate alleles for one child, sometimes on the same nominal
haplotype. J1001 had 29 repeated SNP transcript positions (58 records),
spanning fewer genomic loci. Every child is screened independently after
sample extraction; the filter removes *all* transcript representations of
each implicated genomic site from the child VCF, query, SNP map, and resulting
personalized transcriptome/index. The same rule covers conflicting indels,
while retaining all other sequence-resolved SNPs and indels. Full locus and
genotype details are audited privately in `tx.vcf.excluded_ambiguous.tsv`.
The original cohort VCF is preserved. This exclusion is deliberately
conservative until multiallelic phased genotypes and allele counts are supported.

If the **shared reference receipt is already complete**, do not resubmit the
entire reference stage. The following resume form requires a valid completed
shared reference and submits only missing/incomplete human child preparation
and personalized-index receipts. First submit one child, inspect its logs and
receipts, then submit the rest. Completed child/index pairs are skipped:

```bash
bin/submit_pipeline.sh "$CONFIG" reference --resume --sample J1001
bin/submit_pipeline.sh "$CONFIG" reference --resume --sample J1001 --submit
# WAIT for J1001 child and index DONE; verify both receipts and variant QC
ls -lh "$COHORT_DIR/receipts/J1001/child.rds" "$COHORT_DIR/receipts/J1001/human.rds"

bin/submit_pipeline.sh "$CONFIG" reference --resume --submit
# WAIT for all remaining child/index jobs before attempting pilot
```

`--sample` requires `reference --resume`, and selects the RNA `sample_id`
(e.g. `J1001`), not the infant's `B...` genotype ID. This recovery path is
also available for legacy human trio cohorts; it does not alter public-cross
reference workflows. A normal `reference --submit` still prepares the shared
reference and all child/index jobs from scratch.

After the child/index jobs complete, verify that indels survived selection and
were included in the VCF passed to Personalizer (choose a child with indels in
the full source VCF; a zero count in one child is not by itself a failure):

```bash
awk -F '\t' '!/^#/ && length($4)!=length($5) {n++} END {print "Personalized VCF indels:",n+0}' \
  "$COHORT_DIR/private/J1010/tx.vcf"
```

6. Run pilot, which performs the first direct matrix-based parental assignment:

```bash
bin/submit_pipeline.sh "$CONFIG" pilot --submit
# WAIT; verify pilot receipt and inspect family-specific QC
bin/Rscript_sif poe.R "$CONFIG" status
read -r -p "Select one RNA sample_id (e.g. J1010): " SAMPLE_ID
less -S "$COHORT_DIR/work/$SAMPLE_ID/pilot/"*/assignment_qc.tsv
```

Review `precomputed_origin`, `origin_undefined`, `child_not_phased_heterozygous` and `no_precomputed_variant` statuses and confirm the transcript overlap before full quantification. Never share per-child QC or counts as public cohort aggregates.

7. Once the pilot and parental-origin QC are acceptable, quantify and model:

```bash
bin/submit_pipeline.sh "$CONFIG" scale --submit
# WAIT for quantification / parental assignment DONE
bin/Rscript_sif poe.R "$CONFIG" status
bin/submit_pipeline.sh "$CONFIG" models --submit
# WAIT; verify `models complete`
bin/Rscript_sif poe.R "$CONFIG" status
bin/submit_pipeline.sh "$CONFIG" additional --submit
# WAIT; verify additive gene-level analysis
bin/Rscript_sif poe.R "$CONFIG" additional-status
less -S "$COHORT_DIR/additional_v1/gene_POE.tsv"
```

The new human gene model selects a single locus per gene based on family callability/depth, one offspring/family, controls for fetal cis orientation and configured covariates, uses HC3 uncertainty and ≥2,000 family bootstraps. It requires ≥10 independent families and both cis orientations (at least three each). It reports nominal and BH-adjusted results. The original `model_gene.tsv` remains a `deferred_unimplemented` placeholder as part of the frozen baseline contract; use the new additive gene table instead.

**Privacy:** all GUSTO child/parent genotype files, phase QC, and per-family counts remain under the restricted project data root. Existing `review`/`export` promotes only explicitly permitted **promoter** aggregate columns; it does **not** authorize exporting the new gene estimates. Obtain IRB/data-governance approval and apply the same minimum-cell and disclosure checks to any new gene-level export before cross-cohort use.

## 7. Cross-cohort comparison (after analyses finish)

This is a separate, **aggregate-only** Python analysis. A curated one-to-one ortholog mapping is required; identical gene symbols are not sufficient. No combined individual-level data is constructed. You can begin with only the maize cohorts; add mouse and human when ready. Run in the authorized analysis environment.

1. Create a tab-delimited `cross_manifest.tsv` with fields:

```text
cohort_id  species  tissue  stage  results_path  effect_scale  independence_group  approved_release  minimum_cell_size
```

The actual file must use **tabs**, not spaces. For each cohort, `results_path` is its versioned `additional_v1/gene_POE.tsv`. Use `log2` for public maize/mouse and `ln` for GUSTO. `independence_group` must be a distinct identifier for independent sets of biological units; shared individuals/families must use the same ID. Example rows (replace paths with actual validated paths):

```text
maize_B73_Ki11  maize  endosperm  20_DAP  /.../maize_B73_Ki11/additional_v1/gene_POE.tsv  log2  ki11_20dap  false  0
maize_B73_W22  maize  endosperm  14_DAP  /.../maize_B73_W22/additional_v1/gene_POE.tsv  log2  w22_14dap  false  0
mouse_Babr_CAST  mouse  extraembryonic_ectoderm  E7_5  /.../mouse_Babr_CAST/additional_v1/gene_POE.tsv  log2  babr_e7  false  0
mouse_FVB_CAST_E12  mouse  placenta  E12_5  /.../mouse_FVB_CAST_E12/additional_v1/gene_POE.tsv  log2  fvb_e12  false  0
```

For GUSTO add a human row **only** after gene-aggregate approval: `approved_release=true`, `minimum_cell_size=<approved integer >=10>` and a disclosure-reviewed `results_path`. Never point this manifest to private per-family files. The comparator rechecks `n_families` row-wise and refuses unapproved human input.

2. Create a *tab-delimited* `orthologs.tsv` with columns `species`, `gene_id`, `orthogroup_id`, curated to one-to-one orthologs (no within-species paralogs). It is fine to start with mouse+human or a maize-containing orthology resource; genes without one-to-one mapping are excluded with a reason.
3. Run the comparator with a new, empty versioned output directory:

```bash
python3 "$SCRIPTS_DIR/bin/compare_cohorts.py" \
  --manifest /absolute/path/cross_manifest.tsv \
  --orthologs /absolute/path/orthologs.tsv \
  --out "$DATA_ROOT/cross_cohort/analysis_v1"
```

4. Review `included_gene_effects.tsv`, `excluded_gene_effects.tsv`, `within_stratum_meta.tsv` and `directional_conservation.tsv` under the output directory. A fixed-effect p-value/q-value is reported **only** where multiple independent cohorts measure the *same tissue × stage × species*. The two maize cohorts **will be descriptive separate strata** (14 vs 20 DAP). Directional comparisons across species are descriptive and should be interpreted against 2:1 versus 1:1 genomic dosage, timepoint/tissue differences, proxy genotypes and orthology uncertainty. The comparator does not claim evolutionary conservation significance.

## 8. Verification and troubleshooting

```bash
# 1. Check primary stage status, unchanged by additional analysis:
bin/Rscript_sif poe.R "$CONFIG" status

# 2. Check additive receipt and expected output files:
bin/Rscript_sif poe.R "$CONFIG" additional-status
ls -lh "$COHORT_DIR/additional_v1/"

# 3. Inspect LSF failures, output and dependency conditions:
bjobs -a
ls -lt "$COHORT_DIR/logs/lsf/" | head

```

If `additional_v1` is partially written or the new receipt is `incomplete`, **do not delete or overwrite the preliminary outputs**. Inspect the LSF error, work out why it failed and decide whether to repair the incomplete additional-only files or start a clearly versioned new analysis. No automatic rerun of baseline results is performed. Do not interpret `models complete` as indicating that `additional` was run; it is a separate stage.

For detailed acquisition URLs, repository setup, read-layout specifics and original sample manifests, see `docs/Placenta_EvoDevo_step_by_step.md`. This revised document describes the new post-model analytical layer and the current user-specific status without superseding those source-specific details.
