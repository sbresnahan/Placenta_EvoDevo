# All-SNP baseline NB analysis

## Scope and safety

The baseline analysis uses completed Ornaments quantification products,
reconstructs parental transcript-SNP counts, aggregates gene-SNP counts using
transcript usage weights, and fits each gene and each individual transcript
across its own informative SNPs. Quantification outputs and their receipts are
read, never modified.

Within-gene TE-TSS comparisons require neither a shared callable SNP nor the
same first donor. Reciprocal crosses include `allele * cross`. Interaction
`p < 0.05` is flagged `not_strictly_POE`; a nonsignificant interaction is not
evidence of strict imprinting. Public cohorts remain exploratory.

Private raw products, parent assignments, metadata, manifests, and model input
objects must remain on the deployment. Only the final archive is returned for
reanalysis.

Stages are reused while their receipt and outputs exist. The workflow does not
verify file-content or source-code identity, so it does not detect changed inputs.
Run a known code version on known inputs, and rerun affected stages explicitly
after changing code, configuration, metadata, raw products or resources (see
"Rerunning a stage" below).

## 1. Source, R library and container

Run from the deployment scripts directory (`$SCRIPTS_DIR`) and use the site
configuration and `bin/Rscript_sif`, not a different host R. Inspect
`config/seadragon.sh`, `POE_SITE_FILE` if used, the image path, and the
writable, container-visible `POE_R_LIBRARY`. All compute nodes must see the
same library and source. Keep the source unchanged while cohorts are running.

Install in an allocated interactive/LSF session, not on a restricted login
node. Compiling glmmTMB/TMB can require several GB of RAM.

```bash
export SCRIPTS_DIR=/absolute/path/to/scripts
cd "$SCRIPTS_DIR"
bin/Rscript_sif -e 'install.packages(c("data.table","yaml","lme4","glmmTMB"), lib=.libPaths()[1], repos="https://cloud.r-project.org", Ncpus=1)'
bin/Rscript_sif -e 'library(glmmTMB); library(lme4); library(GenomicRanges); library(Rsamtools); library(rtracklayer); print(sessionInfo())'
bin/Rscript_sif tests/test_baseline_nb.R "$SCRIPTS_DIR"
bin/Rscript_sif tests/test_baseline_workflow.R "$SCRIPTS_DIR"
bin/Rscript_sif tests/test_baseline_adapter.R "$SCRIPTS_DIR"
bin/Rscript_sif tests/test_baseline_between.R "$SCRIPTS_DIR"
bin/Rscript_sif tests/test_baseline_deployment.R "$SCRIPTS_DIR"
```

Each test sources `R/library_paths.R`, so it uses the configured library. If
Matrix ABI warnings occur, rebuild lme4, TMB and glmmTMB against the Matrix
installed inside this same R environment. Do not mix packages built by
different R/container installations. Stop if the required numerical tests are
skipped or fail.

## 2. Cohort configurations and private metadata

Copy each of the six full schema-2 YAMLs, retaining resolved private paths,
acquisition/reference settings, parent definitions, dosage, regulatory
resources, units and libraries. Do not run an example YAML or the baseline
fragment as a complete configuration. Add one `baseline_nb` section to each
copy, following `config/baseline_nb.example.yaml`.

Give each cohort its own absolute `baseline_nb.output_dir`, for example
`/absolute/cohort/results/baseline_nb_v1_run1`. Never share one output
directory between cohorts. `bootstrap` sets the within-gene refits per mixed
gene (default 2000, minimum 200); record the value used.

For GUSTO the model requires complete sex and ancestry metadata. The helper
uses the configured infant-to-FASTQ mapping and exact `J<fastq_id>`
identifiers. The source CSV must have `SubjectID`, `sex`, `mother_ethnicity`.

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

For public cohorts, keep verified `unit_id`, parent and `cluster_id`
assignments. Unknown clusters are explicitly flagged independence-unverified.
Do not fill them with invented biological IDs; set `cluster_id = unit_id` only
when every unit is known to be an independent biological sample.
Sex/ancestry covariates, if available and variable, can be supplied in a
`sample_id` keyed TSV.

Set the remaining five absolute YAML paths:

```bash
export CONFIG_KI11=/absolute/private/config/maize_B73_Ki11.baseline_nb.yaml
export CONFIG_W22=/absolute/private/config/maize_B73_W22.baseline_nb.yaml
export CONFIG_BABR=/absolute/private/config/mouse_Babr_CAST.baseline_nb.yaml
export CONFIG_FVB=/absolute/private/config/mouse_FVB_CAST_E12.baseline_nb.yaml
export CONFIG_FISH=/absolute/private/config/fish_Pprolifica_AB.baseline_nb.yaml
CONFIGS=("$CONFIG_GUSTO" "$CONFIG_KI11" "$CONFIG_W22" "$CONFIG_BABR" "$CONFIG_FVB" "$CONFIG_FISH")
```

Sections 3-6 are run **one cohort at a time**. Do not infer comparability
across different tissues, stages, builds or species.

## 3. Preflight, prepare and measure a pilot

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
Inspect prepare's normalization matrix estimate and `pilot_resources.tsv` in
the baseline output directory. Measure scheduler peak memory with
`bjobs -l JOB_ID` or the site's completed-job accounting. Matrix memory alone
excludes raw/merged tables, copies and fitting overhead. A 4 GiB matrix limit
is enforced by default. If it is exceeded, allocate sufficient memory,
explicitly raise `max_dense_gb` and rerun prepare; do not remove the guard or
run on a login node. Changing batch counts requires rerunning prepare in an
empty output directory.

Preflight requires raw `allele_counts.txt`, `tpms.txt`, valid complete
quantification receipts, parent assignment and reference resources. `prepare`
checks all quantification products/receipts. Gene/event aggregate tables alone
cannot reconstruct individual transcripts. Orientation uses this repository's
REF/ALT parent-assignment implementation.

If raw products are absent or incompatible, run the pipeline's `reference`,
`pilot` and `scale` stages as documented in `docs/OPERATIONS.md` and
`docs/Placenta_EvoDevo_step_by_step.md`, with public retrieval first where
needed. For example, `bin/submit_pipeline.sh "$CONFIG" scale --submit`
requests scale quantification after its prerequisites pass. Review fish
provenance and any failed-library exclusions first. The baseline stages do not
use `models` outputs.

## 4. Fit all gene and transcript batches

Read `gene_batches.tsv` and `transcript_batches.tsv` from the output
directory. They contain the batch indices and expected feature counts. Use the
exact row counts, including any empty batches. The default configuration
creates 20 gene batches and 100 transcript batches.

```bash
# Set these to the prepared manifests' row counts.
export NB_ARRAY_END=20
bin/submit_baseline_nb.sh "$CONFIG" fit gene ARRAY_INDEX
export NB_ARRAY_END=100
bin/submit_baseline_nb.sh "$CONFIG" fit transcript ARRAY_INDEX
unset NB_ARRAY_END
# Only after EVERY array element succeeds:
bin/submit_baseline_nb.sh "$CONFIG" combine
```

Choose memory/wall time/concurrency from the pilot and available site capacity.
Fits are single-threaded. More requested CPUs do not accelerate one fit. A
failed element can be resubmitted alone with its index (`fit gene 7`);
finished batches are skipped.

Inspect `gene_POE.tsv`, `transcript_POE.tsv` and the two availability tables.
Statuses distinguish insufficient samples/SNPs/counts, missing or unidentifiable
covariates, rank deficiency, failed optimization/Hessian and estimable fits.
The gate is >=3 informative SNPs, >=50 rounded counts, >=15 independent
families for human, and >=4 units/clusters with >=2 each reciprocal direction
for public exploratory models. NA uncertainty is never a significant result.

## 5. TE-TSS annotation and within-gene comparisons

Use the prepared, provenance-checked repeat resources from the configured
regulatory analysis. `te-inputs` reads those resources; it does not download
them or silently infer a complete repeat reference.

```bash
bin/submit_baseline_nb.sh "$CONFIG" te-inputs
# After success inspect TE_availability.tsv and within_gene_jobs.tsv.
# If at least one mixed gene is eligible:
bin/submit_baseline_nb.sh "$CONFIG" te-pilot 5
```

The `te-pilot`/`te-fit` bootstrap commands apply to GUSTO descriptive H1
diagnostics. Inspect `TE_pilot_resources.tsv` before requesting the arrays.
Each job resamples biological clusters, refits the eligible isoforms with
their own SNPs and records signed and absolute-effect contrasts.
Draws (`baseline_nb.bootstrap`) are checkpointed in chunks of 25.
Diagnostic intervals and centered-tail probabilities require ≥90% valid refits.
Reciprocal cohorts use `h1-inputs` and `h1-fit` as specified in
`Seadragon_analysis_runbook.md` and `H1_BAYES_MODELS.md`.

```bash
# Set to the number of data rows in within_gene_jobs.tsv.
export NB_ARRAY_END=NUMBER_OF_MIXED_GENE_JOBS
# Increase NB_WALL / NB_MEMORY_GB if the measured TE pilot requires it.
bin/submit_baseline_nb.sh "$CONFIG" te-fit ARRAY_INDEX
unset NB_ARRAY_END
# After EVERY element succeeds, or directly if there are zero mixed genes:
bin/submit_baseline_nb.sh "$CONFIG" te-combine
```

A gene killed at the wall-time limit resumes from its last finished chunk when
resubmitted alone (`te-fit INDEX`). If no mixed genes exist, skip `te-pilot`
and `te-fit`, but still run `te-combine`. Missing fish repeats mean TE analysis
is unavailable, not that fish has no TE-TSS transcripts. Partial references
cannot certify non-TE TSSs. Unknown repeat classes remain distinct.

## 6. Between-gene comparisons and report-ready outputs

Run these sequentially, waiting and inspecting each job:

```bash
bin/submit_baseline_nb.sh "$CONFIG" between-genes
bin/submit_baseline_nb.sh "$CONFIG" context
bin/submit_baseline_nb.sh "$CONFIG" report-ready
bin/submit_baseline_nb.sh "$CONFIG" status
```

`between-genes` uses the propensity/balance, block-bootstrap, permutation and
placebo machinery of the regulatory module, with the POE outcome and depth
inventory taken from the all-SNP gene model. Insufficient exposure groups and
unavailable full repeat references are explicit states. The transcript-class
between-gene comparison uses equal gene weights and is separate from aggregated
gene POE.

The adjusted gene-level exposure uses the active native TSS window.
Its magnitude outcome is `beta_adjusted^2 - SE^2`. The all-transcript
class comparison and within-gene contrasts use absolute beta.

Human xQTL/chromatin annotations provide association context. Transcript-specific
POE-QTL colocalization and causal TE effects require additional measurements.
`context` requires the configured QTL tables
when `regulatory.human_context.require_qtl` is true. Unconfigured resources
remain not assessed. Nonhuman context is explicitly not applicable.

Repeat Sections 3-6 for the other five configs, changing CONFIG and
NB_LOG_DIR. A cohort is complete when `status` reports all seven stages as
`complete_verified`, including explicit unavailable/not-applicable results.

### Rerunning a stage

A stage whose receipt and outputs exist is skipped. To rerun one, move its
receipt (`receipts/<stage>.rds`) and its outputs aside, then submit it again.
After changing code, configuration, metadata or raw products, rerun every
stage that depends on them, or start the cohort in an empty output directory.
After changing repeat or context resources, rerun `te-inputs` and every later
stage so annotations and contrasts stay consistent.

## 7. Create the return archive

Run in an allocated session using the same source, image and library:

```bash
unset NB_ARRAY_END
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="/absolute/private/exports/results_baseline_nb_${STAMP}.tar.gz"
mkdir -p "$(dirname "$OUT")"
bin/Rscript_sif bin/export_results.R --baseline-nb --out "$OUT" "${CONFIGS[@]}"
tar -tzf "$OUT" | head -n 40
```

Every cohort must have all seven stages `complete_verified`. A cohort whose
later stages are still pending (for example fish while repeat discovery runs)
can be exported with `--partial COHORT_ID`: its POE layer is exported once
`prepare` and `combine` verify, plus any later stage that verifies; pending
stages are listed as `pending_not_exported` in `cohorts/<id>/STAGE_STATUS.tsv`
and as not run (not negative) in `AVAILABILITY.tsv`.

Export refuses identifier columns in the result tables. It includes all-SNP
models, annotation/availability, within- and between-gene outputs,
cross-dependent flags, nominal transcript candidates, context, package
versions, cohort/settings/QC summaries, separate cohort summary, documentation
and `MANIFEST.tsv` filenames/sizes. It excludes raw counts, private
metadata/manifests, fit objects, bootstrap records and promoter-model endpoints.

Before sharing, extract to an empty private directory and verify:

```bash
CHECK_DIR="$(mktemp -d)"
tar -xzf "$OUT" -C "$CHECK_DIR"
bin/Rscript_sif -e '
args <- commandArgs(TRUE)
top <- list.dirs(args[1], recursive=FALSE, full.names=TRUE)
stopifnot(length(top)==1L)
m <- data.table::fread(file.path(top,"MANIFEST.tsv"))
actual <- file.info(file.path(top,m$archive_path))$size
stopifnot(all(is.finite(actual)), all(actual == m$bytes))
cat("Inspected",nrow(m),"result paths and sizes; content identity is not verified\\n")
' "$CHECK_DIR"
tar -tzf "$OUT" | grep -E '/private/|/bootstraps/|/fits/|model_promoter.tsv|model_gene.tsv' \
  && printf 'STOP: inspect unexpected private/promoter-model entries\n' || true
```

Also inspect the tables for free-text identifiers and paths. The automated
privacy guard is a safeguard, not a substitute for local privacy review.
Return only `results_baseline_nb_<timestamp>.tar.gz`. Do not
return the source covariate CSV, parent counts or prepared input directories.

## Interpretation on return

Nominal raw `p < 0.05` is the requested selection rule. BH q values are
exported for auditing, not used as an undisclosed filter. Report gene and
transcript effects separately, and do not call TE causality or validated
imprinting from association. The human model lacks explicit cis orientation
adjustment. Small reciprocal cohorts, uncertain biological independence, mixed
fish stages and limited interaction power remain major constraints even when
the optimizer converges.
