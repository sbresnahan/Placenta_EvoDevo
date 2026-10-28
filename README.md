# Placenta_EvoDevo

Parent-of-origin expression analysis of maize endosperm, mouse extraembryonic
tissues, and human GUSTO placenta using Ornaments.

## Start here

Follow the [step-by-step walkthrough](docs/Placenta_EvoDevo_step_by_step.md):
set project paths, select a cohort, check tools, and submit the first stage.

Scripts: `/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/scripts`

Data: `/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/data/<cohort>`

| Cohort | Samples |
|---|---|
| maize_B73_Ki11 | Four PE endosperm pools, 20 DAP |
| maize_B73_W22 | Six PE endosperm pools, 14 DAP |
| mouse_Babr_CAST | Five SE extraembryonic ectoderm units, E7.5 |
| mouse_FVB_CAST_E12 | Four PE female placenta units, E12.5 |
| GUSTO | Local human placenta RNA and linked trios |

## Project stages

| Stage | Action |
|---|---|
| [Setup](00_setup/README.md) | Check native tools and R packages |
| [Retrieval](01_retrieval/README.md) | Download public reads/resources |
| [Reference](02_reference/README.md) | Build transcript and promoter catalogs |
| [Genotypes](03_genotypes/README.md) | Prepare parental markers or child inputs |
| [Ornaments](04_ornaments/README.md) | Build indices, run pilots, quantify RNA |
| [Assignment](05_parental_expression/README.md) | Orient maternal/paternal expression |
| [Models](06_promoter_models/README.md) | Fit eligible gene/promoter comparisons |
| [Export](07_exports/README.md) | Write optional human aggregates |

Stage commands use `bin/submit_pipeline.sh "$CONFIG" STAGE --submit`.
Omitting `--submit` prints a plan. Results include completion status, QC tables,
pilot metrics, and model estimates.

Settings are described in the [schema](config/SCHEMA.md).
See [operations](docs/OPERATIONS.md) for logs/status,
[validation](docs/VALIDATION.md) for software-check coverage, and
[attribution](ATTRIBUTION.md) for software/data credits.
