# Promoter-window H1/H5/H6/H7 analysis on Seadragon

## Scope and prerequisites

This runbook starts from completed all-SNP baseline layers
(`docs/BASELINE_NB_WORKFLOW.md`) and immutable private count batches. Each
promoter-window radius is an independent analysis: radius 3000 bp is primary,
radius 0 and 1000 bp are sensitivities. Window-dependent analyses are computed
per radius; upstream baseline results and external results are read, never
modified. Commands below are operator actions; the repository performs no
remote execution or automatic submission.

Roots (set again in every new shell or allocation):

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export SCRIPTS_DIR="$PROJECT/scripts"
export GXE=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/qtl_results/gxe
export CONFIG_DIR="$PROJECT/data/private/promoter_configs_ppBMI_v1"
export OPERATOR="$PROJECT/data/private/promoter_operator_settings.ppBMI_v1.yaml"
export POE_CMDSTAN_ENV="$(<"$SCRIPTS_DIR/config/cmdstan_env.path")"
unset CMDSTAN
cd "$SCRIPTS_DIR"
```

The six cohort IDs are GUSTO, mouse_Babr_CAST, mouse_FVB_CAST_E12,
maize_B73_Ki11, maize_B73_W22 and fish_Pprolifica_AB. Each baseline layer's
full input configuration is `data/<COHORT>/private/<COHORT>.baseline_nb.yaml`;
the helpers verify the configured locations and compatibility.

Use the cluster `bin/Rscript_sif` and the configured R library. Every
standalone R entry point initializes `R/library_paths.R` before loading
packages or modules. `POE_R_LIBRARY` overrides `config/r_library.path`. An
invalid selected library stops immediately; initialization does not install
packages.

Login nodes are limited to light work (two cores, 12 GB). Run `te-inputs`,
H1 fitting pilots and CmdStan setup on an allocated compute node, for example
an LSF interactive allocation from the site's interactive queue:

```bash
bsub -Is -q "$INTERACTIVE_QUEUE" -n 4 -M 16 -R "rusage[mem=16]" -W 8:00 /bin/bash
```

## 0. Software environment

### 0A. R library and CmdStan tools

Check the selected library:

```bash
bin/Rscript_sif -e 'cat(R.version.string, "\n"); print(.libPaths()); cat("yaml: ", find.package("yaml"), "\n")'
bin/Rscript_sif -e 'for (p in c("brms", "cmdstanr", "posterior")) {stopifnot(requireNamespace(p, quietly=TRUE)); cat(p, as.character(packageVersion(p)), find.package(p), "\n")}'
```

Missing packages are installed into this same shared library through the
site package-management procedure. Do not introduce Conda R or a second R
package library.

CmdStan and its build tools come from a **tool-only** Conda environment
(`poe_cmdstan`). `config/cmdstan_env.path` holds its prefix; `POE_CMDSTAN_ENV`
overrides it. The solve should provide CmdStan and build tools, not `r-base`
or R packages. Keep the prefix on a filesystem bound into the R container
(`/rsrch5` is bound).

```bash
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda create --prefix "$POE_CMDSTAN_ENV" --override-channels -c conda-forge \
  --strict-channel-priority cmdstan cxx-compiler make
mkdir -p "$PROJECT/data/private/runtime_checks"
chmod 700 "$PROJECT/data/private/runtime_checks"
conda list --prefix "$POE_CMDSTAN_ENV" --explicit \
  > "$PROJECT/data/private/runtime_checks/cmdstan_conda_explicit.txt"
unset CMDSTAN
```

No Conda activation is needed for fits. The installation root is
`<prefix>/bin/cmdstan`. [1] The initializer prepends `<prefix>/bin`, sets make
and the Conda C++ compiler, and calls `cmdstanr::set_cmdstan_path()`. It keeps
the R library and does not import host `LD_LIBRARY_PATH` or Python settings.
An explicit `CMDSTAN` pointing elsewhere stops startup. The prefix must be an
absolute path without whitespace or make/shell special characters.

On an allocated compute node, with no fits running, configure the compiler
once, then check and compile a smoke model. `configure` backs up `make/local`,
keeps existing TBB/linker settings, adds a final compiler setting and rebuilds
CmdStan with one core. [2]

```bash
set -o pipefail
bin/Rscript_sif bin/check_cmdstan_runtime.R configure 2>&1 \
  | tee "$PROJECT/data/private/runtime_checks/cmdstan_configure.log"
bin/Rscript_sif bin/check_cmdstan_runtime.R check 2>&1 \
  | tee "$PROJECT/data/private/runtime_checks/cmdstan_check.log"
bin/Rscript_sif bin/check_cmdstan_runtime.R compile 2>&1 \
  | tee "$PROJECT/data/private/runtime_checks/cmdstan_compile.log"
```

Run each command only after the preceding one succeeds. `check` validates
paths, executables and make's selected compiler without modifying anything.
`compile` also runs CmdStanR's toolchain check, compiles a minimal smoke model
and runs its help command; it does not sample. Fitting sessions perform
read-only runtime checks only. If compiler/sysroot or dynamic-library checks
fail, inspect the captured output with site support; do not switch R runtimes
or import a host linker path. Repeat setup after moving, recreating or
updating the tool environment, and keep its explicit package list with the logs.

### 0B. H1 backend check and RStan stack rebuild

brms 2.23.0 samples through CmdStanR and then converts output into RStan
classes, so RStan must load cleanly. [3] Check the backend without fitting:

```bash
bin/Rscript_sif bin/repair_rstan_stack.R check
```

The fitting preflight performs the same namespace loading and empty RStan class
construction before sampling, outside the per-task error handler, and never
installs packages. A passing check does not establish sampler readiness.

If the check fails, the library owner rebuilds RcppParallel 5.1.11-2 (R
displays 5.1.11.2), StanHeaders 2.32.10 and RStan 2.32.7 together under
controlled settings, rather than adding a global Conda linker override. [4, 5]
brms 2.23.0 and CmdStanR 0.99.0 are not replaced. Existing prerequisites must
satisfy the package requirements; a missing dependency or compiler stops the
procedure without an automatic upgrade.

**Login node: source downloads.** Interactive nodes have no internet access.
On an internet-connected login node:

```bash
export RSTAN_SOURCES="$PROJECT/data/private/runtime_checks/rstan_sources"
bin/Rscript_sif bin/repair_rstan_stack.R fetch "$RSTAN_SOURCES"
```

`fetch` downloads the three pinned archives from CRAN and validates archive
layout and package/version metadata. It does not compile, load RStan or modify
the shared packages. Validated archives are reused; a failed transfer is never
published as a completed archive. Stop on a fetch error; a network timeout
does not establish that the pinned version is absent from CRAN.

**Allocated node: build from the local sources.** In a compute allocation,
set the roots again and stage into a new, unused work directory:

```bash
export RSTAN_SOURCES="$PROJECT/data/private/runtime_checks/rstan_sources"
export RSTAN_BUILD="$PROJECT/data/private/runtime_checks/rstan_build_$(date +%Y%m%d)"
bin/Rscript_sif bin/repair_rstan_stack.R stage "$RSTAN_BUILD" --sources "$RSTAN_SOURCES"
```

With `--sources`, the helper copies and validates all three local archives
before the first build and never downloads. Prerequisites are checked in
dependency order. Compilation uses the container gcc/g++ with one compile
process, private Makevars and cleared external include/link/TBB overrides.
Global Makevars, R configuration, permanent library selection and CmdStan
`make/local` are untouched. Builds use a temporary staging library; fresh-process
loading, RStan object construction and linkage inspection must pass before
promotion. The helper reports disk-space estimates and stops when the space is
unavailable. Logs and staging files are kept; inspect `build_<PACKAGE>.log`
after a failure, and use a new work directory for another staging attempt.

Before `apply`, coordinate a maintenance window and stop every session using
the affected packages, including other users' sessions. The flag confirms that
coordination; the helper cannot inspect other users' processes.

```bash
bin/Rscript_sif bin/repair_rstan_stack.R apply "$RSTAN_BUILD" --maintenance-confirmed
bin/Rscript_sif bin/repair_rstan_stack.R check
```

`apply` revalidates staging, verifies the incoming package copies (contents
and permission bits, including nested and empty directories), keeps the
original package directories in a private backup under the shared library,
promotes with same-filesystem directory renames and validates again from a
fresh process using the normal library. A promotion or validation error
attempts automatic rollback. Backups, rejected candidates, logs and
`state.rds` are kept. If `apply` stops at `Package copy verification failed`
before promotion, keep the same work directory and, within the maintenance
window, run `apply` again: it prepares a new private incoming copy and leaves
the earlier copy for inspection. Do not rebuild, delete partial copies, remove
locks or change shared permissions by hand. Copy verification does not
establish shared-library access for other accounts or site ACL behavior.

To restore the prior package directories, stop all affected sessions and run:

```bash
bin/Rscript_sif bin/repair_rstan_stack.R rollback "$RSTAN_BUILD" --maintenance-confirmed
```

Rollback also handles an interrupted promotion with a missing active package
and a matching kept lock. It refuses to overwrite later external package
changes or another transaction's lock. If automatic rollback reports an error,
keep all state and backups and inspect them before manual changes.

Do not run fits until `check` succeeds.

## 1. Configure the promoter-window runs

All radii use inclusive interval intersections clipped to native contigs.
`tss_te_label` and its TE IDs keep exact-overlap definitions;
`promoter_te_label` and `promoter_*` fields record the active window. The
event −500/+100 annotation stays separate.

### 1A. Check the upstream configurations

The Seadragon project file tree lists these paths relative to `$PROJECT`;
readability and scientific validity are checked below and by later preflights. [6]

| Path | Use |
|---|---|
| `data/<COHORT>/private/<COHORT>.baseline_nb.yaml` (six cohorts) | Baseline input configuration |
| `data/GUSTO/private/qtl_context_v1/phenotype_map.tsv` | Phenotype crosswalk; verify mapping coverage |
| `data/GUSTO/private/qtl_context_v1/qtl_inputs.tsv` | Baseline QTL manifest |

Stop and resolve any `MISSING` entry before generation:

```bash
(
  missing=0
  for COHORT in GUSTO mouse_Babr_CAST mouse_FVB_CAST_E12 maize_B73_Ki11 maize_B73_W22 fish_Pprolifica_AB; do
    P="$PROJECT/data/$COHORT/private/$COHORT.baseline_nb.yaml"
    if [[ -f "$P" && -r "$P" ]]; then
      printf 'OK %s\n' "$P"
    else
      printf 'MISSING %s\n' "$P" >&2
      missing=1
    fi
  done
  for NAME in phenotype_map.tsv qtl_inputs.tsv; do
    P="$PROJECT/data/GUSTO/private/qtl_context_v1/$NAME"
    if [[ -f "$P" && -r "$P" ]]; then
      printf 'OK %s\n' "$P"
    else
      printf 'MISSING %s\n' "$P" >&2
      missing=1
    fi
  done
  [[ "$missing" == 0 ]]
)
```

The external `$GXE` directory is outside the tree's scope. Do not use baseline
QTL lead tables or `qtl_inputs.tsv` as GxE inputs.

### 1B. Configure the external ppBMI results

```bash
cd "$SCRIPTS_DIR" && bash bin/Rscript_sif bin/setup_external_ppbmi.R "$PROJECT"
```

The helper reads `config/external_ppbmi.yaml`, which contains the producer
path and GxE root for this deployment. It reads only the small producer
config, runtime exposure summary, aggregate summary and aggregate headers. It
stops on missing files, incompatible headers or disagreement with the supplied
metadata; it does not scan aggregate rows, refit anything or submit jobs.

The live records declare `ppBMI`, `transform=none`, inclusive limits 14–70,
`enabled=1`, runtime `binary=False`, 417 nonmissing and zero missing
observations, mean 24.939616306954434 and SD 5.520665327307657. The producer
uses a pooled continuous z-score across retained ancestries with population SD
(`ddof=0`). [7]

The helper writes four private outputs:

- `data/GUSTO/private/qtl_context_v1/gxe_tables.ppBMI_v1.tsv`
- `data/GUSTO/private/qtl_context_v1/ppBMI_exposure.ppBMI_v1.tsv`
- `data/private/promoter_operator_settings.ppBMI_v1.yaml`
- `data/GUSTO/private/qtl_context_v1/ppBMI_setup.ppBMI_v1.yaml`

Local `enabled=TRUE` means **include the external ppBMI results**; it does not
indicate a ppBMI analysis run within this project. Type, transformation, range
and scaling are populated automatically. The receipt records the metadata,
source paths, small-file sizes, headers and outstanding operator declarations.
It does not certify aggregate rows, mapping coverage or an immutable copy of
the production config.

The helper starts from `data/private/promoter_operator_settings.yaml` if
present, otherwise the template; an operator at another location can be
supplied as the third argument. H1 settings, optional exposure settings and
QTL/BED overrides are kept; only GUSTO's ppBMI manifests and phenotype-map
paths are set. The map is read from the GUSTO baseline YAML. Unmentioned
QTL/BED inputs inherit from that YAML. Identical repeated runs reuse identical
outputs; any differing output stops the helper before it writes. Edit the
written operator directly afterwards. Do not run simultaneous setups for one
project.

### 1C. Review the operator settings

Review `h1_bayes` in `$OPERATOR`. The sampling settings apply to every
reciprocal cohort and radius: 4 chains, 1000 warmup and 1000 post-warmup
iterations, `adapt_delta: 0.99`, `max_treedepth: 12`, `seed: 1147`,
`parallel_chains: 4`. `adapt_delta` must be at least 0.95 and below 1. All
cohorts and radii of one analysis use the same settings; a different setting
is a separate analysis in fresh output directories.

`h1_bayes.count_provenance` and `h1_bayes.producing_pipeline` are optional
descriptive annotations, not model parameters. Missing values, empty values
and the shipped placeholders mean not recorded; missing annotations stay
missing and are never converted to claims of unique assignment or duplicate
correction. Invalid explicit assignment labels are rejected. Required count
support, verified biological clusters, dosage, model, priors and diagnostic
thresholds are unaffected. Optional `exposure_by_cohort` requires actual
H6/H7 inputs; omission keeps the upstream `exposure_nb`. Verify phenotype-map
contents and identifier coverage before the context stage, and review
independently verified clusters per reciprocal direction. Missing fish cluster
support yields explicit ineligibility; no support threshold is relaxed.

### 1D. Generate the configurations

```bash
bash bin/Rscript_sif bin/prepare_promoter_configs.R "$PROJECT" "$OPERATOR" "$CONFIG_DIR"
```

The generator checks the operator mapping, numerical sampling controls,
required placeholders and manifest-path syntax, reads all six upstream YAMLs
and validates merged H1 settings before writing. It requires absolute upstream
paths and an empty destination. It writes 18 configurations
(`<COHORT>.r<RADIUS>.yaml`) and `CONFIGS.tsv`; each sets
`baseline_nb.output_dir` to `data/<COHORT>/results/promoter_r<RADIUS>_run1`.

Review `CONFIGS.tsv`, all 18 YAMLs, repeat-reference scope, dosage, builds,
parental assignments, `h1_bayes` and H6/H7 private input paths. Generation
validates configuration mechanics, not private input provenance, manifest
contents or biological readiness. Preparation never changes source count
batches.

### 1E. Optional unfilled operator scaffold

A separate unfilled operator scaffold is not needed for section 1B. The block
below creates one from the template only when no operator file exists at that
path; it never replaces an existing file.

<!-- operator-create -->
```bash
OPERATOR="$PROJECT/data/private/promoter_operator_settings.yaml"
CONFIG_DIR="$PROJECT/data/private/promoter_configs_run1"
(
  set -e
  umask 077
  mkdir -p "$(dirname "$OPERATOR")"
  if [[ -e "$OPERATOR" || -L "$OPERATOR" ]]; then
    printf 'Keeping existing operator file: %s\n' "$OPERATOR"
  else
    cp -n config/promoter_operator_settings.example.yaml "$OPERATOR"
  fi
  if [[ ! -f "$OPERATOR" || ! -r "$OPERATOR" ]]; then
    printf 'Operator file is not a readable regular file: %s\n' "$OPERATOR" >&2
    exit 1
  fi
)
```

After using this block, set `OPERATOR` and `CONFIG_DIR` again from the roots
at the top before later sections.

## 2. Repeat resources

`te-inputs` reads the repeat track prepared by `regulatory-resources` (and,
for fish, `fish-repeats`) under each configuration's `regulatory.output_dir`.
Check the resource receipts per cohort before window preparation:

```bash
bin/Rscript_sif poe_regulatory.R "$CONFIG_DIR/<COHORT>.r3000.yaml" regulatory-status
```

Resources are shared by the three radii of a cohort. Preparation and fish
repeat discovery are described in `docs/REGULATORY_RESOURCES.md`. Missing
repeats leave TE comparisons unavailable; a partial repeat reference cannot
certify nonTE.

## 3. Prepare one cohort and radius

Select the cohort and radius (set these again in each new shell):

```bash
COHORT=mouse_Babr_CAST
RADIUS=3000
CONFIG="$CONFIG_DIR/$COHORT.r$RADIUS.yaml"
POE_INPUT="$PROJECT/data/$COHORT/private/$COHORT.baseline_nb.yaml"
export NB_LOG_DIR="$PROJECT/data/$COHORT/logs/promoter_r$RADIUS"
OUT=$(awk '/^baseline_nb:/{f=1;next} /^[^ #]/{f=0} f && $1=="output_dir:"{gsub(/["'\'']/,"",$2); print $2; exit}' "$CONFIG"); echo "$OUT"
```

On an allocated compute node:

```bash
bin/Rscript_sif poe_baseline.R "$CONFIG" prepare-window "$POE_INPUT"
bin/Rscript_sif poe_baseline.R "$CONFIG" te-inputs
```

`prepare-window` refuses nonempty destinations. It requires completed upstream
prepare/combine receipts, compatible biological/reference/model settings and
available absolute private batch pointers. It copies allowlisted small POE
tables, inventories, manifests and preparation metadata; the destination
receipts describe this copy. It excludes context, TE comparisons, indexed
jobs, fits and bootstrap chunks. No upstream POE refitting is required.

Review `window_preparation.tsv`, `TE_availability.tsv` and
`transcript_TE_annotation.tsv` in `$OUT` before fitting.

## 4. Reciprocal-cohort Bayesian H1

Reciprocal cohorts are mouse_Babr_CAST, mouse_FVB_CAST_E12, maize_B73_Ki11,
maize_B73_W22 and fish_Pprolifica_AB. The model is specified in
`docs/H1_BAYES_MODELS.md`.

### 4A. Prepare H1 inputs

Submit preparation as one single-host job; the helper derives worker and
memory controls from the requested CPUs and memory and refuses arrays:

```bash
unset NB_ARRAY_END H1_PREP_WORKERS H1_PREP_MEMORY_GB
NB_CPUS=4 NB_MEMORY_GB=16 bin/submit_baseline_nb.sh "$CONFIG" h1-inputs
```

Alternatively run it on an existing four-core, 16-GB allocation:

```bash
H1_PREP_WORKERS=4 H1_PREP_MEMORY_GB=16 bin/Rscript_sif poe_baseline.R "$CONFIG" h1-inputs
```

Use one route, never both, and never run two preparations in one output
directory. Each required count batch is loaded once; workers compute
per-transcript coverage and same-gene consistency checks across batches,
including informative rows on transcripts excluded by the count gates. The
coordinator writes three tables and writes completion last. An incomplete
preparation restarts from the batch indexes; a completed one is reused without
count loading. Requested workers are capped by batch count, scheduler CPU
controls and estimated memory, so the effective count may be lower. The
estimate allows 1 GiB plus four annotation sizes, two feature-index sizes and
one batch-index size for the coordinator, and four times the largest required
batch size plus 256 MiB per worker (1024-based units). This is a planning
guard, not RSS enforcement. Missing batch sizes or a missing budget select
batch-wise serial execution with a message; a one-worker estimate above the
budget stops before loading. No fitting or compilation occurs here.

Review `h1_transcript_coverage.tsv`, `h1_gene_coverage.tsv` and
`h1_bayes_jobs.tsv`. Count-eligible transcripts can include upstream GLMM
failures. Job indices run from 1 to the number of jobs, in gene-identifier
order. Count jobs and gene sizes:

```bash
N=$(( $(wc -l < "$OUT/h1_bayes_jobs.tsv") - 1 )); echo "jobs: $N"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next}{gsub(/"/,""); n++; t=$c["n_count_eligible"]; s+=t; if(t>m)m=t}END{printf "jobs=%d mean_tx=%.1f max_tx=%d\n",n,(n?s/n:0),m}' "$OUT/h1_bayes_jobs.tsv"
```

If there are zero jobs, skip the fit array and run `te-combine` (section 4D)
to write complete unavailable coverage.

### 4B. Pilot and fit array

Every gene job fits five prior profiles on the full data and the primary
profile once per support-preserving cluster omission. Private task objects
allow an interrupted gene to continue from its last completed task with the
same settings and seeds. Each fit uses four simultaneous one-thread chains on
four CPUs; set `parallel_chains` explicitly to match a smaller allocation.

Submit a pilot of the first genes, then inspect diagnostics (section 4C) and
the LSF logs before the full array. As a reference, mouse_Babr_CAST genes
with up to 10 eligible transcripts used at most 3 GB and about 11 minutes per
gene at `adapt_delta` 0.95; 0.99 samples more slowly. These are starting
allocations, not guarantees for other cohorts.

```bash
NB_CPUS=4 NB_MEMORY_GB=8 NB_WALL=4:00 NB_MAX_CONCURRENT=4 NB_ARRAY_END=$(( N < 20 ? N : 20 )) \
  bin/submit_baseline_nb.sh "$CONFIG" h1-fit ARRAY_INDEX
grep -L "Successfully completed" "$NB_LOG_DIR"/h1-fit.*.out
grep -h "Max Memory" "$NB_LOG_DIR"/h1-fit.*.out | sort | uniq -c
grep -h "Run time" "$NB_LOG_DIR"/h1-fit.*.out | awk '{print $4}' | sort -n | tail -1
```

`NB_LOG_DIR` must be set; `ARRAY_INDEX` is passed literally and becomes the
LSF array index. The helper always submits indices `1..NB_ARRAY_END`; genes
already completed exit after printing that they are complete. Size memory and
wall time from the measured pilot, then submit the full array:

```bash
NB_CPUS=4 NB_MEMORY_GB=8 NB_WALL=4:00 NB_MAX_CONCURRENT=16 NB_ARRAY_END=$N \
  bin/submit_baseline_nb.sh "$CONFIG" h1-fit ARRAY_INDEX
bjobs -w | grep nb_h1-fit
```

If a completed gene has recorded `fit_failure` tasks, normal fitting stops for
that gene with an explicit retry instruction. After a passing backend check
(section 0B), run on an allocated node:

```bash
bin/Rscript_sif poe_baseline.R "$CONFIG" h1-retry-failed JOB_INDEX
```

The retry moves the gene's tables and receipt to the private `h1_recovery`
folder, reuses compatible saved tasks and reruns missing or failed tasks. A
completed `diagnostic_failure` task is reused, not refitted. Sampling output
that failed before a task object was saved is not guaranteed recoverable.

### 4C. Completeness and diagnostics

Gene receipts missing (empty output means all present):

```bash
for j in $(seq 1 "$N"); do printf -v f '%s/receipts/h1_%06d.rds' "$OUT" "$j"; [[ -f $f ]] || echo "$j"; done | paste -sd, -
```

Task status by profile (cluster omissions pooled as `LOCO`). Text fields in
the result tables are quoted, so the commands strip quotes first:

```bash
awk -F'\t' 'FNR==1{delete c; for(i=1;i<=NF;i++)c[$i]=i;next}{gsub(/"/,""); o=$c["omitted_cluster"]; k=(o==""||o=="NA")?$c["profile"]:"LOCO"; n[k" "$c["status"]]++}END{for(k in n)print n[k],k}' "$OUT"/h1_bayes/gene_*_diagnostics.tsv | sort -k2,2 -k1,1rn
```

Primary full-data fit per gene:

```bash
awk -F'\t' 'FNR==1{delete c; for(i=1;i<=NF;i++)c[$i]=i;next}{gsub(/"/,"")} $c["profile"]=="primary" && ($c["omitted_cluster"]==""||$c["omitted_cluster"]=="NA"){g=FILENAME; sub(/.*gene_0*/,"",g); sub(/_diag.*/,"",g); printf "%s\t%s\tdiv=%s\ttd=%s\trhat=%s\tess=%s/%s\tebfmi=%s\n",g,$c["status"],$c["divergences"],$c["treedepth_hits"],$c["max_rhat"],$c["min_bulk_ess"],$c["min_tail_ess"],$c["min_ebfmi"]}' "$OUT"/h1_bayes/gene_*_diagnostics.tsv | sort -n
```

Genes with `fit_failure` tasks:

```bash
grep -l fit_failure "$OUT"/h1_bayes/gene_*_diagnostics.tsv | sed 's/.*gene_0*//; s/_diag.*//' | paste -sd, -
```

Ready means Rhat <1.01, bulk/tail ESS ≥400, no divergences or tree-depth hits,
and per-chain E-BFMI ≥0.3 for model parameters and derived contrasts. One
divergence fails readiness even when a percentage display rounds to 0.0%.
`fit_failure`, `diagnostic_failure` and `unavailable_LOCO_support` remain
visible; thresholds are fixed. Different sampler settings form a separate
analysis in a fresh output directory.

### 4D. Combine

Only after every gene receipt exists:

```bash
bin/Rscript_sif poe_baseline.R "$CONFIG" te-combine
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next}{gsub(/"/,""); n[$c["fit_status"]]++}END{for(k in n)print n[k],k}' "$OUT/h1_coverage.tsv" | sort -rn
column -t -s$'\t' "$OUT/h1_cohort_descriptive.tsv" | sed 's/"//g'
```

`h1_coverage.tsv` gives each gene's primary fit status (`ready`,
`diagnostic_failure`, `fit_failure`, `missing_job` or `ineligible`).
`h1_cohort_descriptive.tsv` reports equal-gene averages of full-data posterior
means and medians per profile with the number of ready genes; the averages
include genes that failed diagnostics and carry no combined uncertainty. See
`docs/H1_BAYES_MODELS.md` for priors, estimands, cancellation diagnostics and
interpretation limits.

## 5. GUSTO descriptive H1

Use the GUSTO configuration for the active radius:

```bash
COHORT=GUSTO
CONFIG="$CONFIG_DIR/$COHORT.r$RADIUS.yaml"
export NB_LOG_DIR="$PROJECT/data/$COHORT/logs/promoter_r$RADIUS"
OUT=$(awk '/^baseline_nb:/{f=1;next} /^[^ #]/{f=0} f && $1=="output_dir:"{gsub(/["'\'']/,"",$2); print $2; exit}' "$CONFIG"); echo "$OUT"
bin/Rscript_sif poe_baseline.R "$CONFIG" prepare-window "$PROJECT/data/$COHORT/private/$COHORT.baseline_nb.yaml"
bin/Rscript_sif poe_baseline.R "$CONFIG" te-inputs
bin/Rscript_sif poe_baseline.R "$CONFIG" te-pilot 5
N=$(( $(wc -l < "$OUT/within_gene_jobs.tsv") - 1 )); echo "mixed genes: $N"
```

Size the array from `TE_pilot_resources.tsv` (projected seconds per gene),
then:

```bash
NB_CPUS=1 NB_ARRAY_END=$N bin/submit_baseline_nb.sh "$CONFIG" te-fit ARRAY_INDEX
bin/Rscript_sif poe_baseline.R "$CONFIG" te-combine
```

If there are no mixed genes, skip `te-pilot` and `te-fit` and still run
`te-combine`. A gene stopped at the wall-time limit continues from its last
finished bootstrap chunk when resubmitted alone (`te-fit INDEX`).

The class contrasts use the supplied transcript POE. Bootstrap diagnostics
refit counts from the private batches for this radius; fields named
`*_bootstrap_centered_tail_prob` and intervals require ≥90% valid refits. The
primary magnitude field is `magnitude_sq_corrected` (class mean of
`b^2 - SE^2`), which removes the noise inflation of raw `|b|` when TE and
nonTE transcripts differ in precision; raw `magnitude` is secondary. Review
`TE_mean_se` and `nonTE_mean_se` in the class rows. No H1 discovery threshold
is assigned.

## 6. GUSTO baseline and ppBMI H5

Use the private GxE and ppBMI manifests written in section 1B. Participant-level
exposure matrices, genotype files and phenotype BEDs are not needed for
summary-result integration.

The aggregate summary has 18 merged ancestry/modality groups labelled exactly
`ppBMI` and 396 chromosome-labelled entries; the producer's aggregator glob can
include both. [8] The adapter keeps **only `exposure == "ppBMI"`**. Do not
strip `.chr` suffixes or combine shard rows with merged rows. Header checks and
this summary do not certify every aggregate row; the context stage keeps strict
validation. The adapter reads only these aggregates:

```bash
test -s "$GXE/aggregate/gxe_tier1_results.tsv.gz"
test -s "$GXE/aggregate/gxe_tier1_meta.tsv.gz"
test -s "$GXE/aggregate/gxe_tier2_results.tsv.gz"
bin/Rscript_sif poe_baseline.R "$CONFIG" context
```

Remove an unavailable optional tier from the private manifest and document
the missing source; never point it to a significant-only table or variant
shard. Tier 1 keeps all tested features, b_gi/b_gi_se and nominal versus
feature-calibrated statistics. Primary source support is qval ≤0.05 within
ancestry × modality × exposure. Main b_g/b_e effects remain distinct. Tier 1
ACAT synthesis keeps ancestry-only and discordant-lead features; an IVW
variant effect requires concordant leads. Tier 2 q_meta is an exploratory
source flag because the producer applies BH after selecting the minimum cis
meta-p without within-cis calibration.

Review `gxe_ppBMI_context.tsv.gz`, `gxe_ppBMI_links.tsv.gz`, mapping/resource
coverage and exposure provenance. Baseline and GxE outputs, ancestries, tiers,
modalities, exact-transcript/event/gene contexts and exact/proximal
same-insertion evidence remain separate. Same-insertion requires a variant
inside that insertion. Untested and unavailable evidence stays missing. GUSTO
cohort overlap with the atlas requires independent sample-provenance review;
H5 provides contextual association evidence. See `docs/H5_ppBMI_CONTEXT.md`.

## 7. Between-gene comparison and H6/H7

For each cohort and radius, after its H1 combine:

```bash
bin/Rscript_sif poe_baseline.R "$CONFIG" between-genes
# Every cohort; GUSTO already ran it in section 6. Nonhuman cohorts record
# context as not applicable, which report-ready requires.
bin/Rscript_sif poe_baseline.R "$CONFIG" context
bin/Rscript_sif poe_baseline.R "$CONFIG" report-ready
bin/Rscript_sif poe_baseline.R "$CONFIG" status
```

The between-gene method is described in `docs/REGULATORY_RESOURCES.md`.
`status` lists each baseline stage as `complete_verified` or
`missing_or_stale`.

H6/H7 `baseline_dir` points to this radius-specific POE/annotation layer.
Native family RNA quantification is reused where its catalogue is unchanged.
Feature-family targets, target exclusions, predictor scaling, jobs, fits and
multiplicity corrections are computed in fresh H6/H7 output folders.

```bash
bin/Rscript_sif poe_exposure.R "$CONFIG" preflight
bin/Rscript_sif poe_exposure.R "$CONFIG" prepare
bin/Rscript_sif poe_exposure.R "$CONFIG" pilot 9
```

Review `target_family_links.tsv` for insertion IDs, distances, radius and
exact/proximal-only flags. Same-family insertions produce one transcript-family
target; different families produce separate targets. RNA predictors are
family-wide. The target gene span ±1 kb exclusion is fixed; an intronic
insertion can identify the family while its own reads are excluded.

Use the measured H6/H7 array commands in `EXPOSURE_NB_WORKFLOW.md` for `fit`,
GUSTO `mediation`, `combine`, selected bootstrap diagnostics and `enrichment`.
H7 is maternal PFAS → target-excluded family RNA → signed POE. ppBMI stays in
H5 context. H5 follow-up keeps gene-cluster uncertainty and separate
support/multiplicity families. No H1/H5/H6 significance gate selects H7 jobs.
See `docs/H6_H7_MODELS.md`.

## 8. Completion and export

Verify all planned rows, failures, missing resources and sensitivity
unavailability before interpreting estimates. Run all three radii without
selecting a radius by significance. Export only approved aggregate tables and
review all free text for human-data privacy.

```bash
bin/Rscript_sif poe_exposure.R "$CONFIG" export
```

Use `bin/export_results.R --baseline-nb` for POE/H1/between-gene/H5 aggregates
and `--exposure-nb` for completed H6/H7 aggregate snapshots, following their
documented required files. Raw counts, posterior model objects, batch manifests
with private pointers, participant identifiers and exposure values remain private.

## References

1. CmdStan Installation. [cited 2026 Nov 16]. Available from: https://mc-stan.org/docs/cmdstan-guide/installation.html
2. Install CmdStan or clean and rebuild an existing installation — install_cmdstan • cmdstanr. [cited 2026 Nov 16]. Available from: https://mc-stan.org/cmdstanr/reference/install_cmdstan.html
3. https://raw.githubusercontent.com/cran/brms/2.23.0/R/backends.R. [cited 2026 Nov 16]. Available from: https://raw.githubusercontent.com/cran/brms/2.23.0/R/backends.R
4. https://raw.githubusercontent.com/cran/rstan/2.32.7/src/Makevars. [cited 2026 Nov 16]. Available from: https://raw.githubusercontent.com/cran/rstan/2.32.7/src/Makevars
5. https://raw.githubusercontent.com/cran/RcppParallel/5.1.11-2/R/tbb.R. [cited 2026 Nov 16]. Available from: https://raw.githubusercontent.com/cran/RcppParallel/5.1.11-2/R/tbb.R
6. seadragon_Placenta_EvoDevo_filetree.txt. Seadragon project file tree.
7. https://raw.githubusercontent.com/sbresnahan/Placenta_xQTL_atlas/main/07_gxe_mapping/50_build_gxe_inputs.py. [cited 2026 Nov 16]. Available from: https://raw.githubusercontent.com/sbresnahan/Placenta_xQTL_atlas/main/07_gxe_mapping/50_build_gxe_inputs.py
8. https://raw.githubusercontent.com/sbresnahan/Placenta_xQTL_atlas/main/07_gxe_mapping/57_aggregate_gxe.py. [cited 2026 Nov 16]. Available from: https://raw.githubusercontent.com/sbresnahan/Placenta_xQTL_atlas/main/07_gxe_mapping/57_aggregate_gxe.py
