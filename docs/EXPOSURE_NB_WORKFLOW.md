# H6/H7 workflow

## 1. Runtime and configuration

The entry point is `poe_exposure.R`. Upstream gene/transcript POE commands use
`poe_baseline.R`. Use a fresh output directory for each promoter radius.
See `Seadragon_analysis_runbook.md` for the integrated H1/H5/H6/H7 execution order.

Required R packages: data.table, yaml, glmmTMB, sandwich,
GenomicRanges, IRanges, S4Vectors, Rsamtools and ggplot2.
Use the existing cluster `bin/Rscript_sif` wrapper and site library. Direct local
tests may set `POE_R_LIBRARY` to an existing library. Match compiled Matrix,
TMB, glmmTMB and lme4 dependencies; resolve ABI warnings
in production. Tests here used R 4.4.2 and glmmTMB 1.1.15.2.

Merge `config/exposure_nb.example.yaml` into a complete schema-2 cohort YAML.
The fragment requires a complete cohort configuration. Choose a fresh private output
directory outside `baseline_dir`. Use absolute paths for deployment. All files
under this output, including aggregate results, initially remain private.

## 2. Input-ready checklist

- [ ] Existing gene/transcript private count batches and their `*_batches.tsv`
      and `*_feature_batches.tsv` manifests are available at declared paths.
- [ ] `transcript_TE_annotation.tsv` includes transcript_id, gene_id, chrom,
      start, end, tss, tss_te_label, promoter_te_label, promoter_start,
      promoter_end and te_tss_window_bp; build/contigs match the reference.
- [ ] Explicit sample/unit/family/SubjectID crosswalk agrees with count metadata.
      One RNA library per biological unit; independently verified families.
      Public cross labels are `focal_mother` / `alternate_mother`.
- [ ] Complete native TE catalogue, locus counts and completed library QC exist.
      Missing count rows mean zero ONLY for fully quantified catalogue libraries.
      A missing library is not zero.
- [ ] QC includes mapped **fragments** (not paired-end alignment/read counts),
      layout, strand, retained multimappers, quantifier version and declared
      `annotation_id` (TE catalogue ID/release). If YAML `annotation_id` is set,
      QC must declare that same ID. Keep the catalogue immutable.
- [ ] GUSTO maternal eight-column PFAS table and full private dictionary exist.
      Declare pre-z transform, z-center/scale/reference, original units,
      biospecimen/timing, LOD handling and prior imputation. Do not re-z-score.
      Center and scale refer to the declared pre-z transformed variable (for
      example log concentration if logging preceded standardization), not
      necessarily the raw concentration scale.
- [ ] Freeze annotation memberships before fitting. Fill species, source release,
      native gene namespace and orthology provenance. Use the verified selection
      list in `manifests/exposure_functional_selection.tsv`; it does not itself
      contain gene memberships or pretend a release has been downloaded.
- [ ] GUSTO H5 context is available; the default is the existing
      `context/human_transcript_context.tsv` under baseline_dir.
- [ ] Run a measured library/model pilot, review failure/support distributions,
      and obtain production approval for memory, storage, wall time and concurrency.

No private allele counts, exposure values or completed preprocessing dictionary
were available for this implementation. A public aggregate archive is insufficient
to refit these models. Fish H6 remains blocked until repeat catalogue and exact
TE–TSS annotation are complete. Use GRCh38 for GUSTO, GRCm38 for mouse, B73
RefGen v4 for maize, and the configured GCF_027474105.1 fish reference.

## 3. TE quantification: generate, inspect, pilot

Import already validated locus-resolved counts when possible. Otherwise:

```bash
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml te-plan
```

This **only writes scripts** and a disjoint annotation/QC. No alignment or
quantification is submitted. The libraries template requires absolute BAM or
FASTQ/index paths, layout and strandedness. If `te_loci` is absent, `te-plan`
can construct the native catalogue from `repeat_intervals`; point the YAML
at the resulting catalogue before `prepare`.

The source-compatible strand-aware Telescope pin is
`26a6c38df6d143387799edc2f3e9f976505f093d`
(declared version `1.0.3.1`).
Install this exact checkout into an isolated environment and set
`EX_TELESCOPE_SOURCE` to that checkout. The source requires a compatible
environment, e.g. Python 3.8 and NumPy 1.23.x (it uses `np.int` and `open(...,'rU')`);
validate installation and quantification in the deployment environment.
Freeze all resolved package versions after the real-library
pilot. Do not silently substitute current HEAD: its report/CLI differ.
STAR is pinned to 2.7.11b. See model documentation references for source verification.

Set `EX_MODULE_ROOT` to this checkout and `EX_ALIGN_CPUS` only after pilot sizing.
Generated STAR jobs retain up to 100 mappings by default and NH/HI/AS/nM tags.
They quick-check the BAM, run a 10,000-alignment tag/multimapper probe, and collate
by read name before Telescope `assign --reassign_mode average`. Strand mapping:
PE forward FR / reverse RF; SE forward F / reverse R; unstranded omits the flag.
The probe is not a full-library guarantee: verify header/reference, retained
secondary alignments, mapping ceiling and strand on 2–3 representative libraries.
Review GTF ambiguous-segment exclusions and nested loci before acceptance.

Record `/usr/bin/time -v` peak memory, elapsed time, input/output/temporary disk
bytes, fragment denominator, assignment proportions and mapping/strand QC.
Generated scripts are not a measured resource request and do not submit themselves.

```bash
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml te-import
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml preflight
```

The adapter imports named `transcript`, `final_count`, `unique_count` columns,
skipping the RunInfo comment and `__no_feature`. Retain fractional TE counts.
Telescope's annotation-unique count is stored as `annotation_unique_count`, not
mistaken for genomic unique mapping. To request strict unique-mapper sensitivity,
declare `unique_telescope_report` in the library manifest; generated scripts
then filter genomic alignments to NH=1 and quantify them separately. Import that
report's final_count as `unique_count`, with QC `unique_count_definition:
genomic_NH1`. Use the original mapped-fragment denominator for both predictors.
If a library has no NH1-overlapping fragments and Telescope writes no report,
supply a reviewed header-only completed zero-count report and retain the
no-overlap log/QC; never infer completion from an absent report.
No missing report is imputed. Populate QC after inspection, never by declaring
success from a generated script alone. For counts-only sensitivity imports,
set `family_only: true` and `predictors: [inclusive_family_only]`.

## 4. Prepare, pilot, fit and combine

```bash
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml prepare
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml pilot 9
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml status
```

The run record stores structured settings, dosage and input size/mtime. Changed
inputs/settings require a new output directory. Metadata equality is not
file-content identity: same-size edits with unchanged timestamps can go undetected.
Frozen annotations retain memberships, release metadata and freeze time.
Input batches must not change in place.
The new module does not apply baseline success/significance filters.

Preparation assigns sequential job IDs. Keep each run's input manifests,
settings, predictor cache and indexed model outputs together.

Review job_inventory, preflight, exposure_coverage, TE_scaling, pilot_results and
pilot_resources. Nine pilot jobs spread across input size/order may not cover
every module/level; missing projections remain NA and require extra pilot jobs.
Pilot memory projections report serialized input-object size, not peak RSS.
Use external `/usr/bin/time -v` on representative stages to choose memory.
Inspect the prepared predictor-cache size and temporary disk too.

For each requested batch, invoke applicable stages:

```bash
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml fit 1
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml mediation 1
# Repeat applicable nonempty stages for all batch IDs in private/jobs.tsv.
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml combine
```

`fit` runs H6 and standalone H7 totals; `mediation` runs primary H7 and optional
modification. Not every batch contains both stages. Use private/jobs.tsv to submit
only applicable IDs; a stage with no requested jobs exits with an error.
`combine` blocks until every planned stage file exists and every
job ID occurs exactly once. Failed/unavailable fits remain explicit rows.
An H7 dictionary blocker allows H6-only preparation; adding PFAS later needs a
new run directory and new preparation.

Explicit LSF helper, only after input-ready production approval:

```bash
export EX_LOG_DIR=/absolute/private/logs
export EX_MEMORY_GB=<measured_GB> EX_WALL=<measured_HH:MM>
export EX_CPUS=1 EX_MAX_CONCURRENT=1
bash bin/submit_exposure_nb.sh /absolute/cohort.yaml fit 1
```

The site's existing wrapper interprets `-M`/`rusage[mem]` in GB. Adapt on sites
using MB. For a contiguous applicable range, `EX_ARRAY_END=N` and argument
`ARRAY_INDEX` use LSB_JOBINDEX; otherwise submit explicit IDs. There is no
automatic job chain, no hidden bsub call in R, and no automatic bootstrap.
The new worker is single-threaded. Increase array concurrency only after pilot.

## 5. Bootstrap, biological follow-up, figures, export

```bash
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml bootstrap 1
# Repeat selection_index from bootstrap_selection.tsv, not arbitrary model IDs.
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml enrichment
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml report-ready
bin/Rscript_sif poe_exposure.R /absolute/cohort.yaml export
```

Bootstrap is secondary and top-20 selected by default; 2,000 draws are real
paired model refits and require a cost pilot. Empty discoveries are not a reason
to fabricate a bootstrap selection. `report-ready` writes counts only, not a
scientific report claiming results. Missing H5 context blocks H5 follow-up
explicitly, even if functional results were saved.

Functional import: obtain species-native GO gene associations for the selected
terms and descendants using a pinned GO ontology/annotation release; omit NOT
annotations, record evidence-code policy, and deduplicate set/gene membership.
Map source gene IDs to the exact count-table namespace using a versioned map;
report unmapped and ambiguous IDs. Use Reactome all-level gene memberships for
the verified human pathways. Broad categories >2,000 tested members are flagged,
not silently split after seeing results. Geneimprint's status `Imprinted` is the
prespecified independent set; do not mix `Predicted`, `Conflicting Data`,
`Tissue Dependent` or `Not Imprinted` into it. Species/tissue applicability
still requires review. For maize/fish use curated native resources or an explicit,
versioned orthology mapping; do not assume human placenta/imprinting annotations
apply to endosperm or fish. Keep ambiguous mappings in a coverage audit; no
automatic best-hit assignment is implemented.

PNG helpers in `R/exposure_plots.R` return without drawing an empty forest and
support POE curves, effect forests, signed/magnitude comparisons, PFAS correlation
and joint support. Source the new modules in listed order (see tests), then call
helpers on results or an explicitly refitted supported model. Curves for
mediation hold M at the value supplied to the helper; natural-effect scenarios
are separate outputs. Fonts use Liberation Sans. No biological figure was
generated from unavailable private data.

For pooled H6 correction, source exposure_inputs and exposure_followup/workflow,
then call `ex_pool_h6(c("cohortA/model_results.tsv","cohortB/model_results.tsv"),
"pooled_H6_results.tsv")`. Use all intended cohorts, including unavailable rows.
This computes cross-cohort BH, not a meta-analysis or a shared TE-family mapping.

Only the `aggregate_export` file allowlist is intended for sharing. It removes
participant identifiers/raw predictors/errors and records filenames and byte sizes.
Review source/provenance text and cohort-specific privacy rules before sharing:
it is not anonymization or automatic human-data release approval. Never share
private/prepared.rds, jobs, predictor matrices, raw batch fits or bootstrap draws.

## 6. Independent report archive (separate from the initial report)

H6/H7 are **not added as required parts** of the initial H1/H5 or baseline
report export. Export each cohort's aggregate snapshot first, then package one
or more snapshots independently. Neither command submits production jobs.

Single cohort, after `combine` and any desired optional stages:

```bash
CONFIG=/absolute/private/GUSTO/cohort.yaml
EXPORTS=/absolute/private/exports
mkdir -p "$EXPORTS"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT="$EXPORTS/results_H6_H7_${STAMP}.tar.gz"
bin/Rscript_sif poe_exposure.R "$CONFIG" export
bin/Rscript_sif bin/export_results.R --exposure-nb --out "$OUT" "$CONFIG"
tar -tzf "$OUT" | head -n 40
```

Multiple cohorts, exported separately from their own runs:

```bash
CONFIGS=(/absolute/private/GUSTO/cohort.yaml /absolute/private/mouse/cohort.yaml)
for C in "${CONFIGS[@]}"; do
  bin/Rscript_sif poe_exposure.R "$C" export
done
OUT="/absolute/private/exports/results_H6_H7_$(date +%Y%m%d_%H%M%S).tar.gz"
bin/Rscript_sif bin/export_results.R --exposure-nb --out "$OUT" "${CONFIGS[@]}"
```

The archive mode reads only cohort ID and `exposure_nb.output_dir` from each
YAML (relative paths resolve beside the YAML; literal absolute paths recommended).
It consumes existing `aggregate_export/` files; it does **not** open private
count batches, sample mappings, PFAS inputs, baseline receipts or the initial
report archive. Invoke the archive command directly when a complete aggregate
snapshot is available. Inventories retain allowlisted metadata fields.
Do not put private values into aggregate provenance/free-text fields.

Required snapshot files:

- `MANIFEST.tsv`: unique relative `file` and positive `bytes`;
- `model_results.tsv`, including failed/unavailable rows or a valid header only;
- `model_coverage.tsv`, `job_inventory.tsv`, `provenance.tsv`.

Optional files include pilot resources, TE scaling, exposure coverage, frozen
annotation manifest, bootstrap selection and aggregate summaries, gene evidence,
functional/H5 enrichment and report summary. They must be allowlisted and listed
in the snapshot manifest. Unlisted stale files are not collected. Optional
presence is **not** evidence of a successful analysis, and absence is **not**
a null result. Do not interpret an archive as a `report-ready` completion receipt.
`family_id` is a TE family and is retained; participant IDs, private predictors,
counts and raw errors are rejected.

To include later bootstrap/enrichment stages, first move the existing snapshot
aside inside the private run, then export again. Set `RUN` to that YAML's actual
`exposure_nb.output_dir`; do not move the whole run or edit its private inputs:

```bash
RUN=/absolute/private/GUSTO/results/exposure_nb_v2_run1
mv "$RUN/aggregate_export" "$RUN/aggregate_export_$(date +%Y%m%d_%H%M%S)"
bin/Rscript_sif poe_exposure.R "$CONFIG" export
# Create a NEW archive filename and run --exposure-nb again.
```

Archive layout:

- `cohorts/<ID>/exposure_nb_v1/`: the module's versioned aggregate tables and
  descriptive manifest;
- `AVAILABILITY.tsv`: per-cohort optional snapshot file presence;
- `MANIFEST.tsv`: relative archived filenames and byte sizes (not its own entry);
- `SCRIPTS.tsv`: descriptive code filename/size inventory;
- `docs/` and `README.md`: model, schema and runbook context.

The mode refuses conflicting modes, `--force`, `--partial`, duplicate cohorts,
unsafe paths, symlinks, forbidden columns, missing core tables, size mismatches,
and existing output filenames. It lists and extracts the new archive to confirm
the expected files/sizes and table schemas. Inventories provide **no file-content
identity or snapshot-freshness verification**. Git storage and
upstream package/tool internals are not modified. Keep inputs immutable and
record declared software/annotation releases.

Before local release, extract into an empty private directory and inspect:

```bash
CHECK_DIR="$(mktemp -d)"
tar -xzf "$OUT" -C "$CHECK_DIR"
bin/Rscript_sif -e '
args <- commandArgs(TRUE)
top <- list.dirs(args[1], recursive=FALSE, full.names=TRUE)
stopifnot(length(top)==1L)
m <- data.table::fread(file.path(top, "MANIFEST.tsv"))
sizes <- file.info(file.path(top, m$archive_path))$size
stopifnot(all(is.finite(sizes)), all(sizes == m$bytes))
print(data.table::fread(file.path(top, "AVAILABILITY.tsv")))
' "$CHECK_DIR"
```

Review the tables and free text under local human-data release rules. Automated
screening is not anonymization or release authorization. Return the independent
`results_H6_H7_<timestamp>.tar.gz`, not private run folders.

## 7. Reproduce software validation

```bash
POE_R_LIBRARY=/existing/R/library Rscript tests/test_exposure_modules.R "$PWD" /scratch/validation
POE_R_LIBRARY=/existing/R/library Rscript tests/test_exposure_workflow.R "$PWD" /scratch/validation
POE_R_LIBRARY=/existing/R/library Rscript tests/test_exposure_edges.R "$PWD" /scratch/validation
bash tests/test_exposure_deployment.sh "$PWD"
POE_R_LIBRARY=/existing/R/library Rscript tests/test_exposure_archive.R "$PWD" /scratch/validation
POE_R_LIBRARY=/existing/R/library Rscript tests/test_repository_contract.R "$PWD" /scratch/validation
```

The numerical tests use observed glmmTMB Salamanders counts/environmental values
with artificial labels strictly as software fixtures. No simulated cohort or
biological H6/H7 finding is supplied. The short two-draw bootstrap in the test
is not a production precision recommendation.

Run all `tests/test_baseline*.R` as well. Manifest/provenance assertions use
descriptive fields; scientific fitting parity checks remain. The archive tests
use hand-specified aggregate software records, not biological estimates.
See bundled validation receipts for exact outcomes and unrun checks.

Production validation still requires the private dictionary/counts/crosswalk,
resolved Telescope environment, whole-library alignment/strand review, measured
resource pilot, and explicit cluster submission. Local software tests do not
certify scientific identifiability, causal assumptions, or production readiness.
