# Project inventory

## Entry points and modules

- `bin/submit_pipeline.sh`: host submission interface for retrieve, reference,
  pilot, explicitly accepted scale and models.
- `00_run_pipeline.sh`: optional short LSF driver.
- `run_pipeline.R`: configuration-requiring graph planner; does not call host
  `bsub` from R.
- `poe.R`: metadata/stage/status and local aggregate CLI.
- `runtime.R`: R worker-phase dispatcher.
- `00_setup/` through `07_exports/`: numbered entry points, checked-in LSF
  scripts and module READMEs. See [the module map](README.md#project-stages).
- `bin/`: R Singularity bridge, site initialization, host inventory, native-tool
  wrappers, worker runner and scheduler wrappers.

## Scientific implementation

| Files under `R/` | Responsibility |
|---|---|
| `config.R` | Configuration/path/table validation |
| `acquisition.R` | Bounded verified acquisition and FASTQ checks |
| `reference.R`, `annotation.R` | Reference/transcript/SNP coordinates and promoter catalog |
| `parental.R` | Typed-trio orientation and same-locus allocations |
| `public.R` | Public-cross orientation, coverage-selected sentinels and reciprocal models |
| `models.R` | Human paired-promoter conditional inference |
| `pipeline.R`, `stages.R` | Quantification, assignment, stage contracts and authorization |
| `seadragon.R`, `lsf.R`, `scheduler.R` | Host/R phases, library job graphs and submission interfaces |
| `common.R`, `workflow.R`, `exports.R` | Integrity, receipt and aggregate-privacy utilities |

## Settings, tables and documentation

- `config/`: five example cohort YAMLs, trusted site settings and
  [schema](config/SCHEMA.md).
- `manifests/`: public run/library/unit selections and resource contracts,
  plus empty private schemas. Copy operational tables into the data root.
- `tests/`: deterministic software tests and tiny real-samtools coordinate
  integration when available; not biological cohort results.
- `docs/`: [step-by-step workflow](docs/Placenta_EvoDevo_step_by_step.md), setup, operations and validation.
  Auxiliary TSV inventories/environment metadata are not execution instructions.
- `evidence/`: resource catalog, human recommendations and
  [scientific review](evidence/review.md).

The package contains no participant payloads, raw reads, genomes, variants,
custom human annotation, installed binaries, local/private settings, logs or
completion receipts. Code terms are in `LICENSE`; dependency/data terms are
described in [attribution](ATTRIBUTION.md).
