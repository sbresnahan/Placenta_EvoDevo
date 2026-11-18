# Project files

## Entry points

| Entry | Purpose |
|---|---|
| `tests/test_setup.sh` | Tool and R-package availability |
| `tests/test_inputs.sh` | Per-cohort input readiness after retrieval |
| `bin/submit_pipeline.sh` | Stage submission |
| `poe.R` | Stage, status and aggregate actions |
| `run_pipeline.R` | Job graph planning |
| `runtime.R` | Worker phase dispatch |
| `bin/Rscript_sif` | R container and library initialization |
| `bin/run_worker.sh` | Host/native and R worker phases |
| `poe_regulatory.R` | Repeat/regulatory resource stages and status |
| `poe_baseline.R`, `bin/submit_baseline_nb.sh` | Baseline POE layer, promoter windows, H1, between-gene, H5 |
| `poe_exposure.R`, `bin/submit_exposure_nb.sh` | H6/H7 exposure models |
| `bin/export_results.R` | Promoter-model, baseline and H6/H7 archives |

Numbered module directories contain worker scripts and short READMEs.
See the [module map](README.md#project-stages).

## R implementation

| Files | Purpose |
|---|---|
| `config.R` | Paths, settings and table schemas |
| `acquisition.R` | Downloads and FASTQ integrity |
| `reference.R`, `annotation.R` | Reference coordinates and promoter catalog |
| `parental.R`, `public.R` | Parental orientation and public models |
| `models.R` | Human promoter inference |
| `pipeline.R`, `stages.R` | Quantification, assignment and completion |
| `seadragon.R`, `lsf.R`, `scheduler.R` | Host/R interface and scheduler |
| `common.R`, `workflow.R`, `exports.R` | Shared, completion and aggregate utilities |
| `regulatory_*.R` | Repeat resources, fish repeats, context and between-gene helpers |
| `baseline_*.R`, `window_prepare.R` | Baseline POE layer, TE annotation, within-gene contrasts, context |
| `h1_*.R`, `cmdstan_runtime.R`, `rstan_repair.R` | Bayesian H1 inputs, fits, backend and RStan stack |
| `gxe_ppbmi.R`, `external_ppbmi_setup.R` | ppBMI GxE context |
| `exposure_*.R`, `te_expression.R` | H6/H7 inputs, models, mediation and archive |

## Supporting files

- `config/`: cohort examples, site settings and [schema](config/SCHEMA.md).
- `manifests/`: public run/library/unit tables, resource tables, human schemas.
- `docs/`: [walkthrough](docs/Placenta_EvoDevo_step_by_step.md), [analysis runbook](docs/Seadragon_analysis_runbook.md), model specifications and operations.
- `tests/`: setup and per-cohort input tests ([README](tests/README.md)).
- `evidence/`: dataset catalog and [scientific resource review](evidence/review.md).
- `LICENSE`, [attribution](ATTRIBUTION.md): software/data credits and terms.
