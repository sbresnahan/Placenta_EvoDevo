# Native promoter-window H1/H5/H6/H7 execution

## Scope and prerequisites

Use completed upstream gene/transcript POE models, preparation metadata and
immutable original private count batches. Recompute window-dependent analyses.
Keep all existing external results intact. Commands below are operator actions;
this source package performs no remote execution or automatic submission.

**Resuming after successful ppBMI setup:** apply the H1 metadata-gate repair,
keep the existing private operator/manifests/receipt and unchanged count files,
and continue directly at section 1D. Do not rerun ppBMI setup or fill optional
H1 annotations merely to pass generation. Earlier setup receipts may still
list those fields as blockers; retain that historical receipt unchanged.

**Resuming an interrupted slow H1 preparation:** after the metadata-gate repair,
apply the H1 parallel-preparation repair and use section 3A. Keep the existing
promoter configuration, ppBMI setup and completed upstream stages. Do not
repeat `prepare-window` or delete outputs to restart this incomplete stage.

**Resuming after the RStan/TBB loading error:** use section 0A to validate and
repair the existing shared R library, then section 3B for gene-1 restart or
explicit failed-task retry. Keep completed H1 preparation and all ppBMI files.

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
Every standalone R entry point initializes `R/library_paths.R` before loading
packages or modules. `POE_R_LIBRARY` overrides `config/r_library.path`.
An invalid selected library stops immediately; initialization does not install
packages. Keep the existing R container and shared R package library.

## 0. R library and CmdStan tools

Check the actual selected library before retrying configuration generation:

```bash
bin/Rscript_sif -e 'cat(R.version.string, "\n"); print(.libPaths()); cat("yaml: ", find.package("yaml"), "\n")'
```

If YAML is absent, install it into that same R-compatible shared library through
the existing package-management procedure. Do not introduce Conda R or a second
R package library. Bayesian H1 also needs brms, cmdstanr and posterior in the
selected library:

```bash
bin/Rscript_sif -e 'for (p in c("brms", "cmdstanr", "posterior")) {stopifnot(requireNamespace(p, quietly=TRUE)); cat(p, as.character(packageVersion(p)), find.package(p), "\n")}'
```

Install a **tool-only** environment once. The proposed environment is named
`poe_cmdstan`; it is not assumed to exist already. Review the Conda solve:
it should provide CmdStan and build tools, not `r-base` or R packages.
`config/cmdstan_env.path` is the default prefix. Set `POE_CMDSTAN_ENV` to
override it, or edit that one-line file so fresh batch sessions use the same
selection. Keep the prefix on a filesystem bound into the R container; the
default under `/rsrch5` is already bound.

```bash
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
export POE_CMDSTAN_ENV="$(<config/cmdstan_env.path)"
conda create --prefix "$POE_CMDSTAN_ENV" --override-channels -c conda-forge \
  --strict-channel-priority cmdstan cxx-compiler make
mkdir -p "$PROJECT/data/private/runtime_checks"
chmod 700 "$PROJECT/data/private/runtime_checks"
conda list --prefix "$POE_CMDSTAN_ENV" --explicit \
  > "$PROJECT/data/private/runtime_checks/cmdstan_conda_explicit.txt"
unset CMDSTAN
```

R points to the tools within its own session, as with CMake. No Conda activation
is needed for fits. Conda's installation root is `<prefix>/bin/cmdstan`. [1]
The initializer prepends `<prefix>/bin`, sets make and the Conda C++ compiler,
and calls `cmdstanr::set_cmdstan_path()`. It preserves the R library and does
not import host `LD_LIBRARY_PATH` or Python settings. If an explicit `CMDSTAN`
points elsewhere, startup fails; unset it or select the matching prefix.
The prefix must be an absolute path without whitespace or make/shell special
characters. A missing tool environment does not block non-CmdStan stages.

On an allocated compute node, configure the compiler **once**, with no fitting
jobs running. This backs up `make/local`, preserves existing TBB/linker settings,
adds a final compiler setting, and explicitly rebuilds CmdStan using one core.
The setup modifies the tool environment, not the R library. CmdStanR documents
`make/local` compiler settings and recommends rebuilding after changes. [2]

```bash
set -o pipefail
bin/Rscript_sif bin/check_cmdstan_runtime.R configure 2>&1 \
  | tee "$PROJECT/data/private/runtime_checks/cmdstan_configure.log"
bin/Rscript_sif bin/check_cmdstan_runtime.R check 2>&1 \
  | tee "$PROJECT/data/private/runtime_checks/cmdstan_check.log"
bin/Rscript_sif bin/check_cmdstan_runtime.R compile 2>&1 \
  | tee "$PROJECT/data/private/runtime_checks/cmdstan_compile.log"
```

Run each command only after the preceding command succeeds. `check` validates
paths, executable availability and make's selected compiler without modifying
the installation or rebuilding. `compile` also runs CmdStanR's toolchain check,
compiles a minimal software smoke model and invokes its help command to check
executable loading; it does not sample. Temporary smoke artifacts are removed.
Fitting sessions perform read-only runtime checks, never shared compiler setup.

Do not launch H1 fits until setup and compile-only checks succeed in the real
container. If compiler/sysroot or dynamic-library checks fail, stop and inspect
the captured output with site support; do not switch R runtimes or import a
host linker path as an automatic workaround. Local source-delivery mocks cannot
establish Seadragon compiler/ABI compatibility. Repeat setup/preflight after
moving, recreating or updating the tool environment, and keep its explicit
package list with the logs. After this preflight, run the representative
actual-model pilot in section 3 before sizing production jobs. Keep compiler
output and model objects in private storage.

### 0A. RStan/TBB loading repair and H1 backend preflight

brms 2.23.0 uses CmdStanR for sampling and then converts output into RStan
classes. Successful chains alone do not show that this conversion succeeded.
The reported failure also occurred in a fresh `library(rstan)` session.
The dependency listing showed RcppParallel libtbb.so.2 files and LeafCutter
C++ runtime paths; it did not establish every build-time cause. [3]
The pinned RStan build obtains flags from RcppParallel, whose TBB paths can
be configured externally. The repair builds these packages together under
controlled settings instead of adding a global Conda linker override. [4, 5]

Stop the old fit first. Put the incremental backend repair patch in the scripts
directory; it applies after the H1 parallel-preparation repair:

```bash
git apply --check Placenta_EvoDevo_H1_backend_recovery_fix.patch && git apply Placenta_EvoDevo_H1_backend_recovery_fix.patch
```

For a checkout that already has the backend repair, apply the offline-source
follow-up once. It applies after the backend repair above:

```bash
git apply --check Placenta_EvoDevo_H1_offline_rstan_fix.patch && git apply Placenta_EvoDevo_H1_offline_rstan_fix.patch
```

Keep the repository's `R/`, `bin/` and `poe_baseline.R` together at the checkout
root. On Seadragon that checkout is deployed as `$PROJECT/scripts`; `scripts`
is not an additional directory to create inside the repository.

Check the backend without fitting:

```bash
bin/Rscript_sif bin/repair_rstan_stack.R check
```

This initial check is expected to fail while RStan is broken. The new fitting
preflight performs namespace loading and empty RStan class construction before
new sampling, outside the per-task error handler. It never installs packages.
Passing checks do not establish sampler readiness.

The library owner may stage a clean, same-version repair on an allocated node.
Building requires R 4.3.1 and an allocation identity passed through the
container. It preserves RcppParallel 5.1.11-2 (R displays 5.1.11.2), StanHeaders
2.32.10 and RStan 2.32.7; brms 2.23.0 and CmdStanR 0.99.0 are not replaced.
Existing prerequisites must satisfy the package requirements. A missing
dependency or compiler stops the procedure without an automatic upgrade.

**Login-node phase: source downloads only.** Interactive nodes have no internet
access. If currently inside the failed interactive session, leave it with:

```bash
exit
```

On the internet-connected login node, restore the paths and fetch the exact
archives. Use the existing container wrapper; do not recreate an environment:

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export SCRIPTS_DIR="$PROJECT/scripts"
cd "$SCRIPTS_DIR"
export RSTAN_SOURCES="$PROJECT/data/private/runtime_checks/rstan_sources_v1"
bin/Rscript_sif bin/repair_rstan_stack.R fetch "$RSTAN_SOURCES"
```

Fetch downloads from CRAN and validates archive layout and package/version
metadata. It does not compile, require an allocation, load RStan, or modify the
shared packages. A failed transfer is not published as a completed archive;
already validated archives can be reused. The cache retains `sources.tsv`.
Stop on a fetch error; do not start another allocation until all three sources
have passed. A network timeout does not establish that the pinned version is
absent from CRAN.

**Allocated-node phase: build from local sources.** Request another interactive
allocation using the same successful scheduler command as before; no scheduler
queue or account is guessed here. Repeat any site module setup needed to make
the existing container executable available. Within the new allocation, restore
these variables explicitly:

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export SCRIPTS_DIR="$PROJECT/scripts"
cd "$SCRIPTS_DIR"
export CONFIG="$PROJECT/data/private/promoter_configs_ppBMI_v1/mouse_Babr_CAST.r3000.yaml"
export POE_CMDSTAN_ENV="$(<config/cmdstan_env.path)"
unset CMDSTAN
export RSTAN_SOURCES="$PROJECT/data/private/runtime_checks/rstan_sources_v1"
export RSTAN_REPAIR="$PROJECT/data/private/runtime_checks/rstan_repair_v2"
bin/Rscript_sif bin/repair_rstan_stack.R stage "$RSTAN_REPAIR" --sources "$RSTAN_SOURCES"
```

Keep the failed `rstan_repair_v1` work directory unchanged. Use v2 only if it is
unused; choose another new name for later attempts. No Conda activation,
environment installation, CmdStan reconfiguration, configuration generation,
ppBMI setup or H1 preparation rerun is needed.

With `--sources`, the helper copies and validates all three local source
archives before starting the first build and never falls back to a download.
Missing, empty, invalid or wrong-version sources stop staging. The original
`stage WORK` form still downloads and is unsuitable for disconnected nodes.
Prerequisites are checked in dependency order before each build.
Compilation uses the container gcc/g++ with one compile process, private
Makevars and cleared external include/link/TBB
overrides. Global Makevars, R configuration, permanent library selection and
CmdStan make/local are untouched. Builds use a temporary staging library;
fresh-process loading, RStan object construction and linkage inspection must
pass before promotion. The helper reports disk-space estimates and stops when
the estimated space is unavailable; memory use and build time are not guaranteed.
Logs and staging files are retained. If a build fails, inspect its
`build_<PACKAGE>.log`; do not apply or silently change package versions.
Use a new work location for another staging attempt so earlier logs survive.

Before apply, coordinate a maintenance window and stop all sessions using the
affected packages, including other users' sessions. Stop any earlier repair
process too. The flag below explicitly confirms that coordination; the helper
cannot inspect every other user's process.

```bash
bin/Rscript_sif bin/repair_rstan_stack.R apply "$RSTAN_REPAIR" --maintenance-confirmed
```

Apply revalidates staging, verifies incoming package copies, preserves original
directories in a private backup under the shared library, and promotes them
with same-filesystem directory renames. It validates again from a fresh process
using the normal library. A promotion/validation error attempts automatic
rollback. Backups, rejected candidates, logs and `state.rds` are retained.
Do not run fits until normal-library validation succeeds:

```bash
bin/Rscript_sif bin/repair_rstan_stack.R check
```

To restore the prior package directories, first stop all affected sessions and
any repair process, then explicitly request rollback:

```bash
bin/Rscript_sif bin/repair_rstan_stack.R rollback "$RSTAN_REPAIR" --maintenance-confirmed
```

Rollback also handles an interrupted promotion with a missing active package
and a matching retained repair lock. It refuses to overwrite subsequent external
package changes and refuses another transaction's lock. It restores the old
library state, which was known to be broken; rollback does not make H1 runnable.
If automatic rollback reports an error, keep all state/backups and inspect it
before manual changes. Do not delete package locks or directories blindly.
No actual R 4.3.1 package rebuild was executed as part of local source validation.

## 1. Configure independent window runs

Primary radius is 3000 bp; sensitivities are 0 and 1000 bp. All use inclusive
interval intersections clipped to native contigs. `tss_te_label` and its TE IDs
retain exact-overlap definitions. `promoter_te_label` and `promoter_*` fields
record the active window. The event −500/+100 annotation stays separate.

### 1A. Reuse the recorded upstream paths

The supplied Seadragon file tree lists these paths relative to `$PROJECT`.
It records locations, not current readability or scientific validation. [6]

| Path | Action |
|---|---|
| `data/GUSTO/private/GUSTO.baseline_nb.yaml` | Reuse |
| `data/mouse_Babr_CAST/private/mouse_Babr_CAST.baseline_nb.yaml` | Reuse |
| `data/mouse_FVB_CAST_E12/private/mouse_FVB_CAST_E12.baseline_nb.yaml` | Reuse |
| `data/maize_B73_Ki11/private/maize_B73_Ki11.baseline_nb.yaml` | Reuse |
| `data/maize_B73_W22/private/maize_B73_W22.baseline_nb.yaml` | Reuse |
| `data/fish_Pprolifica_AB/private/fish_Pprolifica_AB.baseline_nb.yaml` | Reuse |
| `data/GUSTO/private/qtl_context_v1/phenotype_map.tsv` | Reuse location; verify mapping coverage |
| `data/GUSTO/private/qtl_context_v1/qtl_inputs.tsv` | Preserve existing baseline QTL configuration |

Check the current files without changing them. Stop and resolve any `MISSING`
entry before generation:

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
QTL lead tables or `qtl_inputs.tsv` as GxE inputs. The versioned files below are
new setup outputs, not paths claimed to exist in the snapshot.

### 1B. Automatically configure the existing external ppBMI results

After installing the ppBMI setup repair, run this single line from any directory.
Copy only the command text, without Markdown fences:

```bash
cd /rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/scripts && bash bin/Rscript_sif bin/setup_external_ppbmi.R /rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
```

No previous ppBMI analysis in Placenta_EvoDevo is required. The helper reads
`config/external_ppbmi.yaml`, which already contains the actual producer path
and GxE root supplied for this deployment. It reads only the small producer
config, runtime exposure summary, aggregate summary and aggregate headers.
It stops on missing files, incompatible headers or disagreement with the supplied
metadata; it does not scan aggregate rows, refit anything or submit jobs.

The supplied live records declare `ppBMI`, `transform=none`, inclusive limits
14–70, `enabled=1`, runtime `binary=False`, 417 nonmissing and zero missing
observations, mean 24.939616306954434 and SD 5.520665327307657. The inspected
producer uses a pooled continuous z-score across retained ancestries with
population SD (`ddof=0`). [7]

The helper writes four private versioned outputs:

- `data/GUSTO/private/qtl_context_v1/gxe_tables.ppBMI_v1.tsv`
- `data/GUSTO/private/qtl_context_v1/ppBMI_exposure.ppBMI_v1.tsv`
- `data/private/promoter_operator_settings.ppBMI_v1.yaml`
- `data/GUSTO/private/qtl_context_v1/ppBMI_setup.ppBMI_v1.yaml`

Local `enabled=TRUE` means **include the existing external ppBMI results**;
it does not assert a completed local ppBMI run. Type, transformation, range and
scaling are populated automatically. The receipt records the metadata, source
paths, small-file sizes, headers and outstanding operator declarations.
It does not certify aggregate rows, mapping coverage or an immutable copy of
the original production config.

The helper uses `data/private/promoter_operator_settings.yaml` if present,
otherwise the supplied template. An existing operator at a different location
can be supplied as the third helper argument, after the project and source-path
YAML. H1 declarations, optional exposure settings and QTL/BED overrides are
preserved; only GUSTO's ppBMI manifests and phenotype-map paths are updated.
The map is read from the live upstream GUSTO YAML. Unmentioned QTL/BED inputs
continue to inherit from that YAML.

Rerunning with identical inputs reuses identical outputs. Any differing output
stops the helper before it writes any other output. After editing the versioned
operator, proceed to the generator; do not rerun setup expecting it to overwrite
your edits. Keep existing outputs if a source changes and review a new version.
Do not launch simultaneous setup invocations for the same project.

### 1C. Reuse unchanged counts without a new annotation audit

`h1_bayes.count_provenance` and `h1_bayes.producing_pipeline` are optional
descriptive annotations, not H1 model parameters. Missing values, empty values
and the exact shipped annotation placeholders mean not recorded. There is no
need to check pipeline versions or supply these labels to integrate ppBMI.
The generator reads an existing versioned operator with those placeholders
directly, without modifying it or rerunning setup.

Existing valid annotations are retained. Per-cohort upstream H1 settings are
inherited, and supplied operator fields are merged without letting an unfilled
annotation erase a recorded value. Missing annotations remain missing in H1
settings and outputs; they are not converted to claims of unique assignment or
duplicate correction. Invalid explicit assignment labels are still rejected.
No count batch, count value or upstream preparation file is changed.

This removes a metadata gate, not scientific input checks. The required count
support, verified biological clusters, dosage, model, priors and diagnostic
thresholds remain unchanged. Optional `exposure_by_cohort` still requires actual
H6/H7 inputs; omission preserves upstream `exposure_nb`.
Verify phenotype-map contents and identifier coverage before the context stage.
Review independently verified clusters per reciprocal direction. Missing fish
cluster support yields explicit ineligibility; no support threshold is relaxed.

### 1D. Generate the independent configurations

With completed ppBMI setup, run this single line. No manual edit of the two
optional H1 annotations is required:

```bash
cd /rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/scripts && bash bin/Rscript_sif bin/prepare_promoter_configs.R /rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo /rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/data/private/promoter_operator_settings.ppBMI_v1.yaml /rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/data/private/promoter_configs_ppBMI_v1
```

The generator checks the operator mapping, numerical sampling controls, genuinely
required placeholders and supplied manifest-path syntax before writing configs.
It reads all six upstream YAMLs and validates merged H1 settings before writing.
It requires absolute upstream paths and an empty config destination; do not erase
a nonempty destination to retry.

Review `CONFIGS.tsv`, all 18 YAMLs, repeat-reference scope, dosage, builds,
parental assignments and H6/H7 private input paths. Successful generation
validates configuration mechanics, not private input provenance, manifest
contents or biological readiness. Later stage-specific preflights still apply.
Preparation never changes source count batches.

Set these shell variables for the later sections, again as a single line:

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo; export SCRIPTS_DIR="$PROJECT/scripts"; export GXE=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/qtl_results/gxe; OPERATOR="$PROJECT/data/private/promoter_operator_settings.ppBMI_v1.yaml"; CONFIG_DIR="$PROJECT/data/private/promoter_configs_ppBMI_v1"; cd "$SCRIPTS_DIR"
```

### 1E. Optional original operator scaffold

This block is retained for operators who explicitly want a separate unfilled
scaffold; it is **not required** by the automatic setup above. It never replaces
an existing original operator file. Creating it does not complete declarations,
and the versioned ppBMI operator is still the input for section 1D.

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

If you use this optional block, restore the versioned `OPERATOR` and `CONFIG_DIR`
variables with the single-line command in section 1D before later stages.

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

### 3A. Batch-wise preparation and restart

Stop an old foreground `h1-inputs` with Ctrl+C before applying this repair.
Wait for the shell prompt; never run two preparations in one output directory.
Place `Placenta_EvoDevo_H1_parallel_prepare_fix.patch` in the scripts directory.
This incremental repair applies after the H1 metadata-gate repair, not to an
unmodified repository checkout. From the scripts directory:

```bash
git apply --check Placenta_EvoDevo_H1_parallel_prepare_fix.patch && git apply Placenta_EvoDevo_H1_parallel_prepare_fix.patch
```

Do not regenerate configurations or rerun ppBMI setup. On the existing
four-core, 16-GB interactive allocation, with `CONFIG` still selecting
`mouse_Babr_CAST.r3000.yaml`, restart in one line:

```bash
H1_PREP_WORKERS=4 H1_PREP_MEMORY_GB=16 bin/Rscript_sif poe_baseline.R "$CONFIG" h1-inputs
```

Each required count batch is loaded once. R workers compute unchanged
per-transcript coverage and retain same-gene consistency checks across batches,
including informative rows on transcripts excluded by the final count gates.
The coordinator writes the same three tables and writes completion last.
An incomplete preparation restarts from the batch indexes; there is no
per-gene checkpoint. A completed preparation is reused without count loading.

Requested workers are capped by batch count, scheduler CPU controls when
available and estimated memory. Progress gives requested/effective workers,
batch start/completion, elapsed seconds and final eligible mixed-class jobs.
The effective count may therefore be below four. No fitting or compilation
occurs in `h1-inputs`.

`H1_PREP_WORKERS` defaults to one. Parallel loading needs an explicit total
`H1_PREP_MEMORY_GB` budget and valid batch `object_bytes`. Missing sizes select
batch-wise serial execution with a clear message; a missing budget also selects
serial execution. This does not restore per-gene batch rereading. The estimate
allows 1 GiB plus four annotation sizes, two feature-index sizes and one
batch-index size for the coordinator; each worker allows four times the largest
required batch size plus 256 MiB. GB controls use 1024-based units. This is a
planning guard, not RSS enforcement; memory or storage measurements may require
fewer workers. A known one-worker estimate above the budget stops before loading.
The wrapper forwards controls through the clean container and restricts native
threads to one per preparation process.

For a new single-host batch allocation, the helper derives controls from
requested CPUs and memory and refuses preparation arrays or excess controls:

```bash
unset NB_ARRAY_END H1_PREP_WORKERS H1_PREP_MEMORY_GB
NB_CPUS=4 NB_MEMORY_GB=16 bin/submit_baseline_nb.sh "$CONFIG" h1-inputs
```

Use either the interactive call or submission, not both. Submission is optional:
login-node preparation remains possible within the site's two-core/12-GB
limits; use controls that fit the actual allocation and stay below the guard.
Do not run actual fitting on a login node. Without explicit controls the
preparation call remains batch-wise serial:

```bash
bin/Rscript_sif poe_baseline.R "$CONFIG" h1-inputs
```

Review `h1_transcript_coverage.tsv`, `h1_gene_coverage.tsv` and
`h1_bayes_jobs.tsv`. Count-eligible transcripts can include upstream GLMM
failures. The job index is the `job_index` column.

### 3B. H1 fitting and diagnostics

After a backend repair, confirm section 0A's read-only check succeeds.
Keep the existing `CONFIG` and preparation outputs.

Run one representative gene on an allocated compute node, measuring elapsed
time, peak RSS and task storage. The normal job runs five prior profiles and
all support-preserving leave-one-cluster-out tasks. Its private task objects
permit restart after a completed task. Include compiler startup in the pilot.
Use four CPUs for four simultaneous one-thread chains. Reduce
`parallel_chains` explicitly to match a smaller allocation.

```bash
bin/Rscript_sif poe_baseline.R "$CONFIG" h1-fit 1
```

An interrupted gene without a completion receipt resumes with that normal
command. Compatible saved tasks are reused; missing tasks run with the same
settings/seeds. Saved table compatibility compares typed columns and row order,
not process-local data.table metadata. If a completed gene has recorded
`fit_failure` tasks, normal fitting now stops with an explicit retry instruction:

```bash
bin/Rscript_sif poe_baseline.R "$CONFIG" h1-retry-failed 1
```

Use the retry action only after backend validation. It preserves the prior
gene outputs/receipt in the private `h1_recovery` folder, reuses compatible saved
tasks and reruns missing/failed tasks. A completed `diagnostic_failure` task is
reused, not automatically refitted; one divergence still fails zero-divergence
readiness even when a percentage display rounds to 0.0%. Sampling output that
failed before a task object was saved is not guaranteed recoverable.
Task progress and errors are printed. A completion receipt means tasks were
recorded, not that every task passed diagnostics.

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

Use the versioned private GxE and ppBMI manifests generated in section 1B; do not
replace them with generic templates. No prior local ppBMI analysis is needed.
Participant-level exposure matrices, genotype files and phenotype BEDs are
unnecessary for summary-result integration.

The supplied aggregate summary has 18 merged ancestry/modality groups labelled
exactly `ppBMI` and 396 chromosome-labelled entries. The inspected aggregator's
glob/parser can include both merged and chromosome-labelled files. [8]
The adapter keeps **only `exposure == "ppBMI"`**. Do not strip `.chr` suffixes or
combine shard rows with merged rows. Header checks and this summary do not
certify every aggregate row; retain the context stage's strict validation.

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

## References

1. CmdStan Installation. [cited 2026 Nov 16]. Available from: https://mc-stan.org/docs/cmdstan-guide/installation.html
2. Install CmdStan or clean and rebuild an existing installation — install_cmdstan • cmdstanr. [cited 2026 Nov 16]. Available from: https://mc-stan.org/cmdstanr/reference/install_cmdstan.html
3. https://raw.githubusercontent.com/cran/brms/2.23.0/R/backends.R. [cited 2026 Nov 16]. Available from: https://raw.githubusercontent.com/cran/brms/2.23.0/R/backends.R
4. https://raw.githubusercontent.com/cran/rstan/2.32.7/src/Makevars. [cited 2026 Nov 16]. Available from: https://raw.githubusercontent.com/cran/rstan/2.32.7/src/Makevars
5. https://raw.githubusercontent.com/cran/RcppParallel/5.1.11-2/R/tbb.R. [cited 2026 Nov 16]. Available from: https://raw.githubusercontent.com/cran/RcppParallel/5.1.11-2/R/tbb.R
6. seadragon_Placenta_EvoDevo_filetree.txt. File provided by the user: seadragon_Placenta_EvoDevo_filetree.txt.
7. https://raw.githubusercontent.com/sbresnahan/Placenta_xQTL_atlas/main/07_gxe_mapping/50_build_gxe_inputs.py. [cited 2026 Nov 16]. Available from: https://raw.githubusercontent.com/sbresnahan/Placenta_xQTL_atlas/main/07_gxe_mapping/50_build_gxe_inputs.py
8. https://raw.githubusercontent.com/sbresnahan/Placenta_xQTL_atlas/main/07_gxe_mapping/57_aggregate_gxe.py. [cited 2026 Nov 16]. Available from: https://raw.githubusercontent.com/sbresnahan/Placenta_xQTL_atlas/main/07_gxe_mapping/57_aggregate_gxe.py
