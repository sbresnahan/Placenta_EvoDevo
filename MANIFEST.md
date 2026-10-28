# Project files

## Entry points

| Entry | Purpose |
|---|---|
| `00_setup/preflight.sh` | Tool and R-package availability |
| `bin/submit_pipeline.sh` | Stage submission |
| `poe.R` | Stage, status and aggregate actions |
| `run_pipeline.R` | Job graph planning |
| `runtime.R` | Worker phase dispatch |
| `bin/Rscript_sif` | R container and library initialization |
| `bin/run_worker.sh` | Host/native and R worker phases |

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

## Supporting files

- `config/`: cohort examples, site settings and [schema](config/SCHEMA.md).
- `manifests/`: public run/library/unit tables, resource tables, human schemas.
- `docs/`: [walkthrough](docs/Placenta_EvoDevo_step_by_step.md), setup, operations,
  software-check coverage and inventories.
- `tests/`: software fixtures and assertions.
- `evidence/`: dataset catalog and [scientific resource review](evidence/review.md).
- `LICENSE`, [attribution](ATTRIBUTION.md): software/data credits and terms.
