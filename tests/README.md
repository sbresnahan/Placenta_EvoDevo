# Tests

Two checks, both run on Seadragon from the scripts directory with
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
