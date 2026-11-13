# Placenta_EvoDevo

Parent-of-origin expression analysis of maize endosperm, mouse extraembryonic
tissues, placental fish, and human GUSTO placenta using Ornaments.

## Start here

Follow the [step-by-step walkthrough](docs/Placenta_EvoDevo_step_by_step.md):
set project paths, select a cohort, check tools, and submit the first stage.

For the additive, non-destructive reanalysis of existing maize results, GUSTO
precomputed parent-of-origin assignment, and cross-cohort comparisons, use the
[revised runbook](docs/REVISED_RUNBOOK.md).

Scripts: `/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/scripts`

Data: `/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/data/<cohort>`

| Cohort | Samples |
|---|---|
| maize_B73_Ki11 | Four PE endosperm pools, 20 DAP |
| maize_B73_W22 | Six PE endosperm pools, 14 DAP |
| mouse_Babr_CAST | Five SE extraembryonic ectoderm units, E7.5 |
| mouse_FVB_CAST_E12 | Four PE female placenta units, E12.5 |
| fish_Pprolifica_AB | *Poeciliopsis prolifica* Line A × Line B; reciprocal placental RNA, sample curation pending |
| GUSTO | Local placenta RNA, phased transcriptomic genotypes and precomputed chromosome parental-origin labels |

The fish profile includes the GEO/SRA study, Dryad fixed-SNP resources, RefSeq
genome and independent Line A PacBio HiFi run. Its initial run tables are
intentionally empty pending an accession-to-cross/tissue audit and a prepared
Line A/Line B SNP table. See the walkthrough before attempting submission.

## Project stages

| Stage | Action |
|---|---|
| [Setup](tests/README.md) | `tests/test_setup.sh`: check native tools and R packages |
| [Retrieval](01_retrieval/README.md) | Download public reads/resources, then `tests/test_inputs.sh "$CONFIG"` |
| [Reference](02_reference/README.md) | Build transcript and promoter catalogs |
| [Genotypes](03_genotypes/README.md) | Prepare parental markers or child inputs |
| [Ornaments](04_ornaments/README.md) | Build indices, run pilots, quantify RNA |
| [Assignment](05_parental_expression/README.md) | Orient maternal/paternal expression |
| [Models](06_promoter_models/README.md) | Fit eligible gene/promoter comparisons |
| [Additional analyses](docs/REVISED_RUNBOOK.md) | Separate TSS/donor questions, relaxed reciprocal callability, optional TE tests and GUSTO gene POE; frozen primary results unchanged |
| [Export](07_exports/README.md) | Write optional human aggregates |

Stage commands use `bin/submit_pipeline.sh "$CONFIG" STAGE --submit`.
Omitting `--submit` prints a plan. Results include completion status, QC tables,
pilot metrics, and model estimates.

For GUSTO recovery after a completed shared reference, use
`reference --resume --sample J1001 --submit` to validate one child, then
`reference --resume --submit` to submit only incomplete child/index jobs.
See the [GUSTO recovery instructions](docs/REVISED_RUNBOOK.md#6-human-gusto--precomputed-phased-transcript-variants-and-chromosome-parental-origin-private).

`STAGE=additional` requires a completed primary `models` receipt and creates
`data/<cohort>/additional_v1/` without changing primary results or receipts.
Use `bin/Rscript_sif poe.R "$CONFIG" additional-status` to verify completion.

Settings are described in the [schema](config/SCHEMA.md).
For opt-in hypothesis #1 (same-first-donor TE/non-TE TSS POE magnitude) and #5
(human cis-QTL/chromatin context), see the [regulatory modules runbook](docs/REGULATORY_MODULES.md).
These use separate `additional_v2` outputs, species-specific resource retrieval,
and separately scheduled fish repeat discovery/masking; existing results are unchanged.
See [operations](docs/OPERATIONS.md) for logs/status and
[attribution](ATTRIBUTION.md) for software/data credits.
