# Tests

Run the two site checks on Seadragon from the scripts directory with
`SCRIPTS_DIR` (and optionally `POE_SITE_FILE`) set as in the
[walkthrough](../docs/Placenta_EvoDevo_step_by_step.md#part-1-shared-setup).

## Setup test (once per environment)

```bash
tests/test_setup.sh
```

Confirms that every tool the pipeline calls is reachable: host commands
(`bash`, `gzip`, `curl`, `bsub`, `singularity`, `/usr/bin/time`), the native
tools activated exactly as the workers activate them (`ornaments`,
`create_personalized_transcriptome`, `samtools`, `bcftools`, `python`), the
`sratoolkit` module, and the R image, R library and R packages. Rerun it after
changing `config/seadragon.sh`, the R image/library, or the native environments.

Expected last line: `Setup test passed: all required tools are available.`

## Input test (per cohort, after retrieval)

```bash
tests/test_inputs.sh "$CONFIG"
```

Reads the cohort YAML and its manifests and confirms that everything the
reference and later stages read is present:

| Profile | Checked |
|---|---|
| Public crosses (maize, mouse, fish) | Manifest linkage; every FASTQ mate at its ENA byte size or SRA receipt size; every required resource (and its decompressed copy); the retrieval receipt; marker sources (and the AGPv3 liftover chain for maize); genome, annotation |
| GUSTO | Manifest/trio linkage; every private local FASTQ; precomputed VCF, origin matrix and FASTQ ID map; trio and covariate tables; genome, annotation, curated transcript FASTA |

Repeat and regulatory resources are retrieved and checked by
`regulatory-resources` itself (see `docs/REGULATORY_RESOURCES.md`).

Expected last line: `<cohort>: all N input checks passed; inputs are ready.`
On failure, every missing or mis-sized input is listed and the script exits 1.

## Local model and archive software checks

Every standalone R script sources `R/library_paths.R` before loading packages
or package-dependent modules. The initializer selects `POE_R_LIBRARY` when set;
direct Rscript calls otherwise read `config/r_library.path`. An invalid selected
path stops startup. Keep site-specific library settings in this shared
configuration.

```bash
bin/Rscript_sif tests/test_r_library_startup.R "$SCRIPTS_DIR"
```

This check validates startup ordering for all standalone R entry points,
explicit and configured-default libraries, invalid-path failure, and generation
of all 18 promoter configurations with `yaml` hidden from the initial search
path. Its container check uses an executable mock without running a container
or contacting Seadragon. The sampler check stops at initialization without
reading a fixture or fitting.

`test_cmdstan_runtime.R SOURCE_ROOT [RECEIPT_DIR]` checks the tool-only Conda
prefix, compiler/make paths, conflicting installations, setup-only rebuild,
preserved TBB/linker settings, direct fitting startup and saved-task resumption.
It runs the real wrapper against an executable container double that starts a
clean R session. CmdStanR calls and compiler/model executables are test doubles:
this test performs no real compilation or sampling and cannot validate site
ABI compatibility. Run `bin/check_cmdstan_runtime.R configure`, then `check`
and `compile` through `bin/Rscript_sif` on an allocated Seadragon node as
documented in the runbook. Only `configure` edits/rebuilds the shared tool
environment; do not run it concurrently with fits.

Run `tests/test_baseline*.R`, `tests/test_exposure_modules.R`,
`tests/test_exposure_workflow.R`, `tests/test_exposure_edges.R`,
`tests/test_exposure_archive.R` and `tests/test_repository_contract.R` using the
same R library. Each accepts the repository path as its first argument; exposure
and repository-contract tests accept a second directory for TSV receipts.
`bash tests/test_exposure_deployment.sh "$PWD"` exercises mocked scheduling only.

The independent archive and repository-contract tests use explicit software
records and local files, not private data, network transfers or production fits.
Numerical model parity tests reuse observed package data with artificial labels;
they do not establish any biological finding. Site setup/input tests still need
the deployment environment and real cohort inputs.

## Promoter-window, H1 and ppBMI contracts

`test_h1_backend_recovery.R SOURCE_ROOT [RECEIPT_DIR]` tests read-only backend
failure handling, stop-before-sampling behavior, saved-table compatibility,
strict settings/mapping checks, interrupted-gene resumption and explicit
failed-task retry with preserved task files and prior failure outputs. Fits
are software test doubles; no sampler, compiler or real backend is invoked.

`test_rstan_repair.R SOURCE_ROOT [RECEIPT_DIR]` tests pinned source archives,
dependency checks, clean build controls, stage/apply separation, maintenance
confirmation, verified promotion, automatic rollback, hard-interruption recovery,
foreign locks and refusal to overwrite later package changes. It also checks
login-node fetch without an allocation or backend loading, validated-cache
reuse, failed and partial downloads, real fetch CLI parsing, offline staging
with a downloader that errors if called, and rejection of missing, empty,
invalid or wrong-version archives before the first build. It uses fake
package directories and build/load callbacks. Private-umask fixtures check
exact file and directory permissions, empty directories, package-root access
bits after promotion, strict rejection of changed content or modes, retained
partial incoming copies and retry from the same staged build without compiling
again. The private process umask and backup parent remain unchanged.
Actual subprocesses check
environment isolation and the fetch CLI; none install or load native R packages. These
tests do not establish that the Seadragon shared library is usable.

`test_h1_parallel_prepare.R SOURCE_ROOT [RECEIPT_DIR]` compares a frozen
per-gene reference loop with batch-wise one-worker and two-worker execution.
All three TSVs must match byte-for-byte, with additional in-memory checks of
coverage types and missing values. Fixed software fixtures cover cross-batch
genes, missing transcripts, support exclusions, rounding, cluster conflicts
on count-ineligible transcripts, invalid indexes/counts, worker errors,
memory/CPU controls and completed-stage reuse. Per-read event files count
actual count-file opens, including worker processes. The small instrumented
benchmark measures fixture overhead and I/O reuse, not Seadragon throughput.

`test_h1_prepare_deployment.R SOURCE_ROOT [RECEIPT_DIR]` checks forwarding
through the clean container, native thread limits, resource-derived scheduler
defaults, single-host placement and rejection of preparation arrays or
controls above allocated CPUs/memory. It uses executable argument mocks only.

`test_external_ppbmi_setup.R SOURCE_ROOT [RECEIPT_DIR]` exercises metadata-only
setup using local software fixtures and the observed small summary records
supplied for this deployment. It checks automatic type/range/scale declarations,
header compatibility, source and override preservation, exact merged labels,
optional H1 annotations, private outputs and identical/conflicting reruns.
Aggregate bodies deliberately contain an unparsed software sentinel: the helper
must not read result rows. This test does not validate actual aggregate rows,
phenotype-map contents or count provenance. `test_gxe_contract.R` additionally
checks that the real adapter excludes chromosome-labelled duplicate rows before
row validation and rejects shard-only input.

`test_h1_optional_annotations.R SOURCE_ROOT [RECEIPT_DIR]` uses a
two-placeholder operator with absent upstream H1 blocks. It generates all 18
configs and initializes the 15 reciprocal H1 settings without fitting or compiler
access. It verifies unrecorded optional metadata, exact valid-setting identity,
cohort-specific inheritance, retained numerical checks and unchanged count
sentinels, upstream YAMLs, operator, ppBMI manifests and setup receipt.
No scientific count data are created or modified by these setup fixtures.

`test_promoter_operator_setup.R SOURCE_ROOT [RECEIPT_DIR]` executes the exact
operator-file creation block from the runbook. It checks private file creation,
existing-file preservation, field-specific preflight failures before writes,
all-six-source validation, generation of 18 configurations and inheritance of
baseline QTL/BED and per-cohort exposure settings. All paths and declarations
used for generation in this test are explicit local software fixtures, not
evidence about the real private configurations. The template's recorded
phenotype-map location is checked as a location only, not as validated content.

Run `test_promoter_contract.R`, `test_promoter_deployment.R`,
`test_h1_bayes_contract.R` and `test_gxe_contract.R` with the repository path.
The H1 draw-contrast checks additionally use posterior when installed.
`test_h1_workflow.R` requires the Bayesian packages and the configured CmdStan
tool environment; it deliberately fails model invocations to verify complete
failure coverage.

`test_h1_sampler.R SOURCE_ROOT PRIVATE_TEST_DIR OBSERVED_COUNT_FIXTURE.rds`
compiles and samples the actual joint model. The fixture uses observed
`glmmTMB::Salamanders` counts with artificial software labels and must contain
`data` and `annotation`. Set `H1_SMOKE_FULL=true` for four chains with
1000 warmup and 1000 retained iterations each. The default short check uses
250 each and can fail computational readiness while passing backend execution.
Keep model objects and fixture posterior values in the private test directory.
Only validation receipts and diagnostics belong in a source delivery.
