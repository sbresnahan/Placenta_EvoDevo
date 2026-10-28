# Setup

Use the [walkthrough](Placenta_EvoDevo_step_by_step.md#part-1-shared-setup)
to create a cohort's working directories and configuration.

## Tool check

```bash
00_setup/preflight.sh "$CONFIG"
```

Expected output: `Required tools and R packages are accessible.`
`--execute` is a compatible alias.

| Component | Setting |
|---|---|
| Native environment | `/rsrch5/home/epi/bhattacharya_lab/software/ornaments/bin/activate` |
| samtools / bcftools | `samtools-1.16.1` / `bcftools-1.16` |
| R image | `/risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif` |
| R library | `/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1` |
| R binds | `/rsrch5`, `/rsrch9`, `/tmp` |
| R packages | data.table, yaml, curl, GenomicRanges, IRanges, S4Vectors, Rsamtools, Biostrings |

Site settings are in `config/seadragon.sh`; cohort paths/settings are in
`config/<cohort>.private.yaml`. Tool checks work before cohort input tables
are populated.

## Settings

Select actual file paths, strand modes, biological cluster IDs, and model
covariates in the cohort's configuration and tables. Existing approval fields
and extra legacy sign-off columns are ignored.

The [schema](../config/SCHEMA.md) describes each input table.
Reference and analysis stages validate the inputs they use.
