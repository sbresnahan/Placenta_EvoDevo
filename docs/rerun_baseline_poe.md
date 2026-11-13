# All-SNP baseline NB rerun

## Scope and safety

This is a new, versioned analysis. It does not replace or delete existing
quantification, legacy POE tables, or their receipts. It reuses completed raw
Ornaments products, reconstructs parental transcript-SNP counts, aggregates
gene-SNP counts using transcript usage weights, and fits each gene and each
individual transcript across its own informative SNPs.

Within-gene TE-TSS comparisons require neither a shared callable SNP nor the
same first donor. Reciprocal crosses include `allele * cross`. Interaction
`p < 0.05` is flagged `not_strictly_POE`; a nonsignificant interaction is not
evidence of strict imprinting. Public cohorts remain exploratory.

Private raw products, parent assignments, metadata, manifests, and model input
objects must remain on the deployment. Only the strict final archive should
be returned for reanalysis. The old selected-SNP results are not validation of
this model.

## 1. Apply the patch to the pinned source

Use a clean checkout, or a new worktree if the deployment has local edits.
Do not discard local changes.

```bash
git clone https://github.com/sbresnahan/Placenta_EvoDevo.git Placenta_EvoDevo_baseline_nb
cd Placenta_EvoDevo_baseline_nb
git checkout ccfa3c825d941b08c6941b92fd7990996bb4cd74
git apply --check /absolute/downloads/baseline_poe_model_patch.diff
git apply /absolute/downloads/baseline_poe_model_patch.diff
git diff --check
export SCRIPTS_DIR="$PWD"
```

The comparison reference is `sbresnahan/placenta-poe-pfas` at
`c285f131071becec7347c830d263e78a88ac45de`. Do not apply this patch blindly
to a later upstream revision. Resolve and test any upstream conflicts first.
Keep this checkout unchanged throughout the rerun.

## 2. R library and container

Use the existing site configuration and `bin/Rscript_sif`, not a different
host R. Inspect `config/seadragon.sh`, `POE_SITE_FILE` if used, the image path,
and the writable, container-visible `POE_R_LIBRARY`. Preserve the existing
deployment library, or select a new compatible library in the site file.
All compute nodes must see the same library and source.

Install in an allocated interactive/LSF session, not on a restricted login
node. Compiling glmmTMB/TMB can require several GB of RAM.

```bash
bin/Rscript_sif -e 'install.packages(c("digest","data.table","yaml","lme4","glmmTMB"), lib=.libPaths()[1], repos="https://cloud.r-project.org", Ncpus=1)'
bin/Rscript_sif -e 'library(glmmTMB); library(lme4); library(GenomicRanges); library(Rsamtools); library(rtracklayer); library(digest); print(sessionInfo())'
bin/Rscript_sif tests/test_baseline_nb.R "$SCRIPTS_DIR"
bin/Rscript_sif tests/test_baseline_workflow.R "$SCRIPTS_DIR"
bin/Rscript_sif tests/test_baseline_adapter.R "$SCRIPTS_DIR"
bin/Rscript_sif tests/test_baseline_between.R "$SCRIPTS_DIR"
bin/Rscript_sif tests/test_baseline_deployment.R "$SCRIPTS_DIR"
```

The tests themselves do not source the deployment bootstrap. The wrapper
passes `R_LIBS_USER` explicitly. If Matrix ABI warnings occur, rebuild lme4,
TMB and glmmTMB against the Matrix installed inside this same R environment.
Do not mix packages built by different R/container installations. Stop if
the required numerical tests are skipped or fail. See `validation_summary.md`
in the supplied bundle for the tested package versions and limitations.

## 3. Create new full cohort configurations and private metadata

Copy each of the six actual full schema-2 YAMLs, retaining resolved private
paths, acquisition/reference settings, parent definitions, dosage, regulatory
resources, units and libraries. Do not run an example YAML or the baseline
fragment as a complete configuration. Add one `baseline_nb` section to each
copy, following `config/baseline_nb.example.yaml`.

Choose a NEW absolute `baseline_nb.output_dir` for each cohort, for example
`/absolute/cohort/results/baseline_nb_v1_run1`. Do not change the old
`paths.model_gene` or the raw quantification paths. Never share one output
directory between cohorts. Keep bootstrap at 2000 unless deliberately
documenting a lower-resolution exploratory run (minimum 200).

For GUSTO the literal model requires complete sex and ancestry metadata.
The archived covariate placeholder is insufficient. The helper uses the
configured infant-to-FASTQ mapping and exact `J<fastq_id>` identifiers.
The source CSV must have `SubjectID`, `sex`, `mother_ethnicity`.

```bash
export CONFIG_GUSTO=/absolute/private/config/GUSTO.baseline_nb.yaml
export GUSTO_COVARIATES=/absolute/private/actual_covariates.csv
export GUSTO_NB_METADATA=/absolute/private/manifests/GUSTO_baseline_metadata.tsv
umask 077
bin/Rscript_sif bin/prepare_baseline_gusto_metadata.R \
  "$CONFIG_GUSTO" "$GUSTO_COVARIATES" "$GUSTO_NB_METADATA"
```

Set `baseline_nb.metadata` in that YAML to the generated TSV. This file is
private. The helper refuses missing/duplicate mappings and refuses overwrite.
Inspect covariate coding and its biological meaning. The field called
`ancestry` here is the supplied maternal-ethnicity variable, not inferred
genetic ancestry. No label imputation is performed.

For public cohorts, preserve verified `unit_id`, parent and `cluster_id`
assignments. Unknown clusters are explicitly flagged independence-unverified.
Do not fill them with invented biological IDs. Sex/ancestry covariates, if
available and variable, can be supplied in a `sample_id` keyed TSV.

Set the remaining five absolute YAML paths:

```bash
export CONFIG_KI11=/absolute/private/config/maize_B73_Ki11.baseline_nb.yaml
export CONFIG_W22=/absolute/private/config/maize_B73_W22.baseline_nb.yaml
export CONFIG_BABR=/absolute/private/config/mouse_Babr_CAST.baseline_nb.yaml
export CONFIG_FVB=/absolute/private/config/mouse_FVB_CAST_E12.baseline_nb.yaml
export CONFIG_FISH=/absolute/private/config/fish_Pprolifica_AB.baseline_nb.yaml
CONFIGS=("$CONFIG_GUSTO" "$CONFIG_KI11" "$CONFIG_W22" "$CONFIG_BABR" "$CONFIG_FVB" "$CONFIG_FISH")
```

The rest of this runbook is executed **one cohort at a time**. Start with
GUSTO, then repeat for each other config. Do not infer comparability across
different tissues, stages, builds or species.

## 4. Preflight, prepare and measure a pilot

The scheduler helper submits only the requested action. It does not create
dependencies or wait for jobs. Wait for DONE, inspect stderr/stdout and
receipts after EVERY command before proceeding to its dependent command.
Scheduler success alone does not establish model convergence.

```bash
export CONFIG="$CONFIG_GUSTO"
export NB_LOG_DIR=/absolute/private/logs/GUSTO_baseline_nb_run1
export NB_QUEUE=medium NB_CPUS=1 NB_MEMORY_GB=16 NB_WALL=12:00
export NB_MAX_CONCURRENT=10
unset NB_ARRAY_END
bin/submit_baseline_nb.sh "$CONFIG" preflight
# After successful completion:
bin/submit_baseline_nb.sh "$CONFIG" prepare
# After successful completion:
bin/submit_baseline_nb.sh "$CONFIG" pilot 10
```

These memory/time settings are starting requests, not measured guarantees.
Inspect prepare's normalization matrix estimate and `pilot_resources.tsv`
in the selected baseline output directory. Measure scheduler peak memory
with `bjobs -l JOB_ID` or the site's completed-job accounting. Matrix memory
alone excludes raw/merged tables, copies and fitting overhead. A 4 GiB
matrix limit is enforced by default. If exceeded, allocate sufficient memory,
choose a NEW output directory and explicitly raise `max_dense_gb`; do not
remove the guard or run on a login node. Changing batch counts also requires
a new prepared directory.

Preflight requires raw `allele_counts.txt`, `tpms.txt`, valid complete
quantification receipts, parent assignment and reference resources.
`prepare` checks all quantification products/receipts. Existing gene/event
aggregate tables alone cannot reconstruct individual transcripts.
The new orientation path deliberately reuses this repository's REF/ALT
parent-assignment implementation rather than the reference's L/R schema.

If raw products are absent or incompatible, stop this route. Use the existing
pipeline's `reference`, `pilot` and `scale` stages as documented in
`docs/OPERATIONS.md` and `docs/Placenta_EvoDevo_step_by_step.md`, with public
retrieval first where needed. For example,
`bin/submit_pipeline.sh "$CONFIG" scale --submit` requests upstream scale
quantification after its prerequisites pass. Review fish provenance and any
failed-library exclusions first. Do not rerun retrieval/quantification merely
because old sentinel model outputs are unsuitable.

## 5. Fit all gene and transcript batches

Read `gene_batches.tsv` and `transcript_batches.tsv` from the NEW output
directory. They contain the batch indices and expected feature counts.
Use the exact row counts, including any empty batches. The default configuration
creates 20 gene batches and 100 transcript batches.

```bash
# Replace these values if the prepared manifests differ.
export NB_ARRAY_END=20
bin/submit_baseline_nb.sh "$CONFIG" fit gene ARRAY_INDEX
export NB_ARRAY_END=100
bin/submit_baseline_nb.sh "$CONFIG" fit transcript ARRAY_INDEX
unset NB_ARRAY_END
# Only after EVERY array element succeeds:
bin/submit_baseline_nb.sh "$CONFIG" combine
```

Choose memory/wall time/concurrency from the pilot and available site capacity.
Fits are single-threaded. More requested CPUs do not accelerate one fit.
Each batch checks raw/config/source fingerprints, which adds shared-filesystem
I/O. Avoid aggressive concurrency, especially for GUSTO.

Inspect `gene_POE.tsv`, `transcript_POE.tsv` and the two availability tables.
Statuses distinguish insufficient samples/SNPs/counts, missing or unidentifiable
covariates, rank deficiency, failed optimization/Hessian and estimable fits.
The gate is >=3 informative SNPs, >=50 rounded counts, >=15 independent
families for human, and >=4 units/clusters with >=2 each reciprocal direction
for public exploratory models. NA uncertainty is never a significant result.

## 6. TE-TSS annotation and within-gene comparisons

Use the prepared, provenance-checked repeat resources from the configured
regulatory analysis. `te-inputs` reads those resources; it does not download
them or silently infer a complete repeat reference.

```bash
bin/submit_baseline_nb.sh "$CONFIG" te-inputs
# After success inspect TE_availability.tsv and within_gene_jobs.tsv.
# If at least one mixed gene is eligible:
bin/submit_baseline_nb.sh "$CONFIG" te-pilot 5
```

Inspect `TE_pilot_resources.tsv` before requesting the full bootstrap arrays.
The pilot reports projected time for one gene, not all loci; transcript-rich
genes can take much longer. Each production job jointly resamples biological
clusters (within reciprocal direction for public cohorts), refits every
eligible isoform using its own SNPs and retains signed and absolute-effect
contrasts. The default is 2000 draws, checkpointed in chunks of 25. >=90%
valid draws are required for inferential summaries.

```bash
# Set to the number of data rows in within_gene_jobs.tsv, not 2000.
export NB_ARRAY_END=NUMBER_OF_MIXED_GENE_JOBS
# Increase NB_WALL / NB_MEMORY_GB if the measured TE pilot requires it.
bin/submit_baseline_nb.sh "$CONFIG" te-fit ARRAY_INDEX
unset NB_ARRAY_END
# After EVERY element succeeds, or directly if there are zero mixed genes:
bin/submit_baseline_nb.sh "$CONFIG" te-combine
```

If no mixed genes exist, skip `te-pilot` and `te-fit`, but still run
`te-combine`. Missing fish repeats mean TE analysis is unavailable, not that
fish has no TE-TSS transcripts. Partial references cannot certify non-TE TSSs.
Unknown repeat classes remain distinct.

## 7. Between-gene comparisons and report-ready outputs

Run these sequentially, waiting and inspecting each job:

```bash
bin/submit_baseline_nb.sh "$CONFIG" between-genes
bin/submit_baseline_nb.sh "$CONFIG" context
bin/submit_baseline_nb.sh "$CONFIG" report-ready
bin/submit_baseline_nb.sh "$CONFIG" status
```

`between-genes` reuses the existing propensity/balance, block-bootstrap,
permutation and placebo machinery, but the POE outcome and depth inventory
come exclusively from the new all-SNP gene model. Insufficient exposure
groups and unavailable full repeat references are explicit states.
The transcript-class between-gene comparison uses equal gene weights and is
separate from aggregated gene POE.

For compatibility the adjusted gene-level exposure is still any TE-overlapping
multi-exon TSS event versus all such TSSs non-TE. Its magnitude outcome is
`beta_adjusted^2 - SE^2`, not `abs(beta_adjusted)`. The all-transcript
class comparison and within-gene contrasts instead use absolute beta.
Keep these estimands separate.

Human xQTL/chromatin annotations are context, not proof of transcript-specific
POE, colocalization or causality. Unconfigured resources remain not assessed.
Nonhuman context is explicitly not applicable.

Repeat Sections 4-7 for the other five configs, changing CONFIG and NB_LOG_DIR.
All seven statuses must be `complete_verified` before export, including
explicit unavailable/not-applicable results. Never remove receipts or run
legacy `models` to fill missing baseline results. Source/R-package/config/metadata/raw
changes require a NEW baseline output directory. Repeat/context resource
changes require a new directory and rerun so stale annotations cannot be
mixed with finished contrasts.

## 8. Create the return archive

Run in an allocated session using the same checkout, image and library:

```bash
unset NB_ARRAY_END
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="/absolute/private/exports/results_baseline_nb_${STAMP}.tar.gz"
mkdir -p "$(dirname "$OUT")"
bin/Rscript_sif bin/export_results.R --baseline-nb --out "$OUT" "${CONFIGS[@]}"
tar -tzf "$OUT" | head -n 40
sha256sum "$OUT"
```

Do not use `--force`. Export refuses incomplete/stale stages and refuses
individual-identifier columns in the allowlisted result tables. It includes
all-SNP models, annotation/availability, within- and between-gene outputs,
cross-dependent flags, nominal transcript candidates, context, package/code
provenance, cohort/settings/QC summaries, separate cohort summary,
documentation and `MANIFEST.tsv` checksums. It excludes raw counts, private
metadata/manifests, fit objects, bootstrap records and all legacy sentinel
POE endpoints.

Before sharing, extract to a NEW private directory and verify:

```bash
CHECK_DIR="$(mktemp -d)"
tar -xzf "$OUT" -C "$CHECK_DIR"
bin/Rscript_sif -e '
args <- commandArgs(TRUE)
top <- list.dirs(args[1], recursive=FALSE, full.names=TRUE)
stopifnot(length(top)==1L)
m <- data.table::fread(file.path(top,"MANIFEST.tsv"))
actual <- vapply(file.path(top,m$archive_path), digest::digest, "",
                 algo="sha256", file=TRUE)
stopifnot(identical(unname(actual), as.character(m$sha256)))
cat("Verified",nrow(m),"result checksums\n")
' "$CHECK_DIR"
tar -tzf "$OUT" | grep -E '/private/|/bootstraps/|/fits/|model_promoter.tsv|model_gene.tsv' \
  && printf 'STOP: inspect unexpected private/legacy entries\n' || true
```

Also inspect the tables for free-text identifiers and paths. The automated
privacy guard is a safeguard, not a substitute for local privacy review.
Return only `results_baseline_nb_<timestamp>.tar.gz` and its SHA256. Do not
return the source covariate CSV, parent counts or prepared input directories.

## Interpretation on return

Nominal raw `p < 0.05` is the requested selection rule. BH q values are
exported for auditing, not used as an undisclosed filter. Report gene and
transcript effects separately, and do not call TE causality or validated
imprinting from association. The literal human model lacks explicit cis
orientation adjustment. Small reciprocal cohorts, uncertain biological
independence, mixed fish stages and limited interaction power remain major
constraints even when the optimizer converges.
