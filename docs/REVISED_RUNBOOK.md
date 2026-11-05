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
| Q6. Do TSS-overlapping TEs associate with altered transcript parental imbalance? | Cohorts with matching repeat annotation | New `event_TE_annotation.tsv`, `TE_comparison.tsv` | Annotation and matched contrasts run only with a verified build-matched UCSC RepeatMasker table; otherwise explicitly `not_assessed`. Exploratory, not causal. |
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
00_setup/preflight.sh "$CONFIG"
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
less -S "$COHORT_DIR/additional_v1/TE_comparison.tsv"
```

Expect **`additional_v1 complete`**. Inspect full status and direct-coordinate support, and report nominal p < 0.05 separately from `q_value <= 0.05`. No 29-file retrieval recheck and no new HapMap3 download is needed for this step.

## 3. Maize B73 × W22 — 14 DAP endosperm

Last verified status was **`reference complete`**, with **pilot blocked** because it had not run. Confirm current status first—work may have progressed. Existing reference reuses the *Ki11* HapMap3 AGPv3 VCFs, not Ki11-specific B73/Ki11 marker calls. The working W22 marker manifest must have **11 lines** and columns `source_id path build processed_path`.

```bash
use_cohort maize_B73_W22
bin/Rscript_sif poe.R "$CONFIG" status
wc -l "$COHORT_DIR/manifests/maize_B73_W22_marker_sources.tsv"
00_setup/preflight.sh "$CONFIG"
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
00_setup/preflight.sh "$CONFIG"
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
00_setup/preflight.sh "$CONFIG"
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

## 6. Human GUSTO — phased/imputed trios and placenta (private)

**No existing GUSTO preliminary outputs need preserving.** GUSTO now includes chromosome-wise origin assignment during `assign_sample`; gene-level models are added only as an analysis *after* primary models.

```bash
use_cohort GUSTO
00_setup/preflight.sh "$CONFIG"
```

1. Edit the private copies of `GUSTO_sources.tsv`, `GUSTO_libraries.tsv`, `GUSTO_units.tsv`, `GUSTO_trios.tsv`, `GUSTO_covariates.tsv` in `$COHORT_DIR/manifests/`; verify child/sample/family linkage, strand, local RNA paths and read counts. Use one primary library per child and independent families for inference.
2. Set `paths.typed_vcf` to the **phased/imputed trio VCF** containing child, mother and father in the *same GRCh38 coordinate/REF/ALT representation*. Set `paths.personalization_vcf` to the suitable child imputed VCF (it can be the same normalized imputation resource, if appropriate). Do not assume imputed GP/DS alone is a hard-called phased genotype; `GT` values must be present.
3. In **`$CONFIG` (the private YAML)** add or confirm the top-level block below. It is already present in updated `config/GUSTO.example.yaml`; copying it is necessary if the private YAML predates this patch:

```yaml
phase_origin:
  enabled: true
  min_informative: 20
  min_agreement: 0.95
  window_bp: 5000000
  min_window_votes: 5
  min_covered_windows: 3
```

GUSTO's additional gene-wide family bootstrap can take longer than the original short model stage. If your private YAML predates the patch, add this optional block **under `lsf.stages`** before submission (not at the top level):

```yaml
    additional:
      memory_value: 16.0
      reservation_value: 16.0
      cpus: 2
      wall_time: '12:00'
      queue: medium
```

The method uses **all** eligible heterozygous phased fetal sites with one or both homozygous informative parents (not random 10k SNP panels). For each chromosome, it votes for whether fetal haplotype 1 or 2 is maternal, requires ≥20 informative sites spanning ≥3 separately supported 5-Mb windows and ≥95% agreement, and excludes conflicting windows. An uncallable chromosome falls back **only to unambiguous direct Mendelian sites**. Suspected switches cannot be repaired by assigning a single label to the whole chromosome. **Important:** the current query captures `GT`, not FORMAT `PS`; if the imputation represents separate phase sets or switches, validate phase sets upstream or implement a per-phase-set adapter before treating a chromosome as callable. Do not assume this QC proves chromosome-long phasing.

4. Validate VCF samples and build; the `bcftools` command must be available on a node where the native environment is initialized:

```bash
# Inspect VCF sample IDs and ensure every trio ID is present.
# Example using the available site bcftools wrapper, as configured:
"$SCRIPTS_DIR/bin/native/bcftools" query -l "$COHORT_DIR/private/typed.vcf.gz" | head
```

5. Build the reference and per-child indices; use the existing GRCh38 genome and HPLRv2.0 transcript GTF paths from the GUSTO example YAML:

```bash
bin/submit_pipeline.sh "$CONFIG" reference --submit
# WAIT for reference DONE, verify
bin/Rscript_sif poe.R "$CONFIG" status
```

6. Run pilot, which also performs the **first chromosome-origin assignment**:

```bash
bin/submit_pipeline.sh "$CONFIG" pilot --submit
# WAIT; verify pilot receipt and inspect family-specific QC
bin/Rscript_sif poe.R "$CONFIG" status
read -r -p "Select one RNA sample_id: " SAMPLE_ID
less -S "$COHORT_DIR/private/$SAMPLE_ID/typed.query.tsv.phase_origin.tsv"
```

Review the `status`, `n_informative`, `agreement`, `n_covered_windows` and `n_conflicting_windows` per chromosome for *all* children. Investigate failed chromosomes, parent mismatches, sex chromosomes, and phase-set boundaries **before** committing to full quantification. These files are private and must never be shared as a cohort aggregate.

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

# 4. Run synthetic R model/phase tests (only on a host with the R container):
bin/Rscript_sif tests/test_additional.R

# 5. Run aggregate-only Python tests:
python3 -m unittest discover -s tests -p test_compare_cohorts.py -v
```

If `additional_v1` is partially written or the new receipt is `incomplete`, **do not delete or overwrite the preliminary outputs**. Inspect the LSF error, work out why it failed and decide whether to repair the incomplete additional-only files or start a clearly versioned new analysis. No automatic rerun of baseline results is performed. Do not interpret `models complete` as indicating that `additional` was run; it is a separate stage.

For detailed acquisition URLs, repository setup, read-layout specifics and original sample manifests, see `docs/Placenta_EvoDevo_step_by_step.md`. This revised document describes the new post-model analytical layer and the current user-specific status without superseding those source-specific details.
