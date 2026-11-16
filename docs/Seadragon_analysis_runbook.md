# Native promoter-window H1/H5/H6/H7 execution

## Scope and prerequisites

Use completed upstream gene/transcript POE models, preparation metadata and
immutable original private count batches. Recompute window-dependent analyses.
Keep all existing external results intact. Commands below are operator actions;
this source package performs no remote execution or automatic submission.

Roots:

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export SCRIPTS_DIR="$PROJECT/scripts"
export GXE=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/qtl_results/gxe
cd "$SCRIPTS_DIR"
```

The six cohort IDs are GUSTO, mouse_Babr_CAST, mouse_FVB_CAST_E12,
maize_B73_Ki11, maize_B73_W22 and fish_Pprolifica_AB.
Upstream inputs are `data/<COHORT>/results/baseline_nb_v1_run1`;
the corresponding full input configuration is
`data/<COHORT>/private/<COHORT>.baseline_nb.yaml`.
The helper verifies the actual configured locations and compatibility.

Use the cluster `bin/Rscript_sif` and configured R library for NB stages.
Bayesian H1 additionally requires brms, cmdstanr, posterior, a C++ toolchain
and compiled CmdStan available within that same runtime. Set `CMDSTAN` to
the CmdStan installation accessible in the container. Check packages, compiler
and a compiled sampler on an allocated compute node before fitting private data.
Keep compiler output and model objects in private storage.

## 1. Configure independent window runs

Primary radius is 3000 bp; sensitivities are 0 and 1000 bp. All use inclusive
interval intersections clipped to native contigs. `tss_te_label` and its TE IDs
retain exact-overlap definitions. `promoter_te_label` and `promoter_*` fields
record the active window. The event −500/+100 annotation stays separate.

Create a private operator settings YAML with:

- `h1_bayes`: reviewed count provenance and producing pipeline, plus the
  settings in `config/promoter_h1.example.yaml`.
- `gusto_human_context`: `gxe_tables`, `exposure_manifest` and the verified
  `phenotype_map`, retaining the original baseline QTL and BED inputs.
- Optional `exposure_by_cohort`: complete per-cohort H6/H7 configuration
  fragments from `config/exposure_nb.example.yaml`.

Count provenance choices are `fragment_unique_assignment`,
`upstream_duplicate_corrected` or `assignment_covariance_unknown`.
The last records unresolved assignment covariance. Known fragment duplication
requires upstream correction. Review independently verified clusters per
reciprocal direction. Missing fish cluster support yields explicit ineligibility.

```bash
OPERATOR="$PROJECT/data/private/promoter_operator_settings.yaml"
CONFIG_DIR="$PROJECT/data/private/promoter_configs_run1"
bin/Rscript_sif bin/prepare_promoter_configs.R "$PROJECT" "$OPERATOR" "$CONFIG_DIR"
```

The helper requires absolute upstream paths and an empty config destination.
Review `CONFIGS.tsv`, all 18 YAMLs, repeat-reference scope, dosage, builds,
parental assignments and H6/H7 private input paths. Preparation never changes
the source count batches.

## 2. Reuse the completed POE layer

For one cohort/radius:

```bash
COHORT=mouse_Babr_CAST
RADIUS=3000
CONFIG="$CONFIG_DIR/$COHORT.r$RADIUS.yaml"
POE_INPUT="$PROJECT/data/$COHORT/private/$COHORT.baseline_nb.yaml"
export NB_LOG_DIR="$PROJECT/data/$COHORT/logs/promoter_r$RADIUS"
bin/Rscript_sif poe_baseline.R "$CONFIG" prepare-window "$POE_INPUT"
bin/Rscript_sif poe_baseline.R "$CONFIG" te-inputs
```

`prepare-window` refuses nonempty destinations. It requires completed upstream
prepare/combine receipts, compatible biological/reference/model settings and
available absolute private batch pointers. It copies allowlisted small POE
tables, inventories, manifests and preparation metadata; the destination
receipts describe this reuse operation. It excludes downstream context,
TE comparisons, indexed jobs, fits and bootstrap chunks. No upstream POE
refitting is required.

Review `window_preparation.tsv`, `TE_availability.tsv` and
`transcript_TE_annotation.tsv` before fitting. Missing repeats leave TE
comparisons unavailable. A partial repeat reference cannot certify nonTE.

## 3. Reciprocal-cohort Bayesian H1

```bash
bin/Rscript_sif poe_baseline.R "$CONFIG" h1-inputs
```

Review `h1_transcript_coverage.tsv`, `h1_gene_coverage.tsv` and
`h1_bayes_jobs.tsv`. Count-eligible transcripts can include upstream GLMM
failures. The job index is the `job_index` column.

Run one representative gene on an allocated compute node, measuring elapsed
time, peak RSS and task storage. The normal job runs five prior profiles and
all support-preserving leave-one-cluster-out tasks. Its private task objects
permit restart after a completed task. Include compiler startup in the pilot.
Use four CPUs for four simultaneous one-thread chains. Reduce
`parallel_chains` explicitly to match a smaller allocation.

```bash
bin/Rscript_sif poe_baseline.R "$CONFIG" h1-fit 1
```

Inspect all diagnostics before setting production resources:

```bash
export NB_CPUS=4
export NB_MEMORY_GB=16
export NB_WALL=24:00
export NB_MAX_CONCURRENT=2
# Set after inspecting the actual job inventory and measured pilot.
export NB_ARRAY_END=NUMBER_OF_H1_JOBS
bin/submit_baseline_nb.sh "$CONFIG" h1-fit ARRAY_INDEX
unset NB_ARRAY_END
# Run after every scheduled gene job has produced its coverage outputs.
bin/Rscript_sif poe_baseline.R "$CONFIG" te-combine
```

The displayed memory/time are starting allocations requiring pilot review.
If there are zero eligible genes, skip the fit array and run `te-combine`
to write complete unavailable coverage.

Ready means Rhat <1.01, bulk/tail ESS ≥400, no divergence/tree-depth hits,
and per-chain E-BFMI ≥0.3 for model parameters and derived contrasts.
`fit_failure`, `diagnostic_failure` and `unavailable_LOCO_support` remain
visible. Diagnose the task and configure longer runs explicitly in a fresh
output directory; thresholds are fixed. See `H1_BAYES_MODELS.md` for priors,
estimands, cancellation diagnostics and interpretation limits.

## 4. GUSTO descriptive H1

Use the GUSTO destination configuration for the active radius.

```bash
COHORT=GUSTO
CONFIG="$CONFIG_DIR/$COHORT.r$RADIUS.yaml"
export NB_LOG_DIR="$PROJECT/data/$COHORT/logs/promoter_r$RADIUS"
bin/Rscript_sif poe_baseline.R "$CONFIG" te-pilot 5
# Inspect TE_pilot_resources.tsv; size the array from within_gene_jobs.tsv.
export NB_CPUS=1
export NB_ARRAY_END=NUMBER_OF_GUSTO_MIXED_GENES
bin/submit_baseline_nb.sh "$CONFIG" te-fit ARRAY_INDEX
unset NB_ARRAY_END
bin/Rscript_sif poe_baseline.R "$CONFIG" te-combine
```

Reuse supplied transcript POE for class contrasts. Bootstrap diagnostics
are recomputed from private counts for this radius. Their field names include
`bootstrap_centered_tail_prob`; ≥90% valid refits are required for diagnostic
uncertainty. H1 discovery thresholds are not assigned.

## 5. GUSTO baseline and ppBMI H5

Copy and complete `manifests/gxe_tables.template.tsv` and
`manifests/ppBMI_exposure.template.tsv` in private storage. Declare the actual
producing config/provenance, enabled ppBMI row, continuous/binary definition,
transformation, scaling and range. A filename cannot establish these properties.
The inspected producer scales continuous exposures pooled across ancestries;
the production manifest must declare its actual settings. Participant-level
exposure matrices, genotype files and phenotype BEDs are unnecessary here.

The adapter reads only these authoritative aggregates:

```bash
test -s "$GXE/aggregate/gxe_tier1_results.tsv.gz"
test -s "$GXE/aggregate/gxe_tier1_meta.tsv.gz"
test -s "$GXE/aggregate/gxe_tier2_results.tsv.gz"
bin/Rscript_sif poe_baseline.R "$CONFIG" context
```

Remove an unavailable optional tier from the private manifest and document
the missing source; never point it to a significant-only table or variant shard.
Tier 1 preserves all tested features, b_gi/b_gi_se and nominal versus
feature-calibrated statistics. Primary source support is qval ≤0.05 within
ancestry × modality × exposure. Main b_g/b_e effects remain distinct.
Tier1 ACAT synthesis retains ancestry-only and discordant-lead features;
an IVW variant effect requires concordant leads. Tier 2 q_meta is an
exploratory source flag because the producer applies BH after selecting
the minimum cis meta-p without within-cis calibration.

Review `gxe_ppBMI_context.tsv.gz`, `gxe_ppBMI_links.tsv.gz`,
mapping/resource coverage and exposure provenance. Baseline and GxE outputs,
ancestries, tiers, modalities, exact-transcript/event/gene contexts and
exact/proximal same-insertion evidence remain separate. Same-insertion
requires a variant inside that insertion. Untested and unavailable evidence
stays missing. GUSTO cohort overlap with the atlas requires independent
sample-provenance review; H5 provides contextual association evidence.

## 6. Between-gene and H6/H7 stages

```bash
bin/Rscript_sif poe_baseline.R "$CONFIG" between-genes
bin/Rscript_sif poe_baseline.R "$CONFIG" report-ready
bin/Rscript_sif poe_baseline.R "$CONFIG" status
```

H6/H7 `baseline_dir` must point to this radius-specific POE/annotation layer.
Keep native family RNA quantification immutable where its catalogue is
unchanged. Recompute feature-family targets, target exclusions, predictor
scaling, jobs, fits and multiplicity corrections in fresh H6/H7 output folders.

```bash
bin/Rscript_sif poe_exposure.R "$CONFIG" preflight
bin/Rscript_sif poe_exposure.R "$CONFIG" prepare
bin/Rscript_sif poe_exposure.R "$CONFIG" pilot 9
```

Review `target_family_links.tsv` for insertion IDs, distances, radius and
exact/proximal-only flags. Same-family insertions produce one transcript-family
target; different families produce separate targets. RNA predictors are
family-wide. The target gene span ±1 kb exclusion remains fixed; an intronic
insertion can identify the family while its own reads are excluded.

Use the measured H6/H7 array commands in `rerun_exposure_nb.md` for `fit`,
GUSTO `mediation`, `combine`, selected bootstrap diagnostics and `enrichment`.
H7 is maternal PFAS → target-excluded family RNA → signed POE. ppBMI stays
in H5 context. H5 follow-up retains gene-cluster uncertainty and separate
support/multiplicity families. No H1/H5/H6 significance gate selects H7 jobs.

## 7. Completion and export

Verify all planned rows, failures, missing resources and sensitivity
unavailability before interpreting estimates. Repeat across all three radii
without selecting a radius using significance. Export only approved aggregate
tables and review all free text for human-data privacy.

```bash
bin/Rscript_sif poe_exposure.R "$CONFIG" export
```

Use `bin/export_results.R --baseline-nb` for POE/H1/H5 aggregates and
`--exposure-nb` for completed H6/H7 aggregate snapshots, following their
documented required files. Raw counts, posterior model objects, batch manifests
with private pointers, participant identifiers and exposure values remain private.
