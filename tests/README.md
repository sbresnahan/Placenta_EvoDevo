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

Configured `additional_analysis.te_rmsk` is also checked. Regulatory-module
resources are retrieved and checked by `regulatory-resources` itself.

Expected last line: `<cohort>: all N input checks passed; inputs are ready.`
On failure, every missing or mis-sized input is listed and the script exits 1.

## Local model and archive software checks

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

Run `test_promoter_contract.R`, `test_promoter_deployment.R`,
`test_h1_bayes_contract.R` and `test_gxe_contract.R` with the repository path.
The H1 draw-contrast checks additionally use posterior when installed.
`test_h1_workflow.R` requires the Bayesian packages and a configured CmdStan
path; it deliberately fails model invocations to verify complete failure coverage.

`test_h1_sampler.R SOURCE_ROOT PRIVATE_TEST_DIR OBSERVED_COUNT_FIXTURE.rds`
compiles and samples the actual joint model. The fixture uses observed
`glmmTMB::Salamanders` counts with artificial software labels and must contain
`data` and `annotation`. Set `H1_SMOKE_FULL=true` for four chains with
1000 warmup and 1000 retained iterations each. The default short check uses
250 each and can fail computational readiness while passing backend execution.
Keep model objects and fixture posterior values in the private test directory.
Only validation receipts and diagnostics belong in a source delivery.
