# Placenta_EvoDevo

Parent-of-origin expression analysis of maize endosperm, mouse extraembryonic
tissues, placental fish, and human GUSTO placenta using Ornaments.

## Start here

Follow the [step-by-step walkthrough](docs/Placenta_EvoDevo_step_by_step.md):
set project paths, select a cohort, check tools, and submit the first stage.

For the all-SNP baseline analysis (gene and transcript POE, TE-TSS within- and
between-gene comparisons, human regulatory context) and its results archive,
use the [baseline runbook](docs/rerun_baseline_poe.md); outputs are described in
the [baseline schema](docs/BASELINE_NB_SCHEMA.md).

For completed upstream POE layers, follow
`docs/Seadragon_analysis_runbook.md`: create independent radius-3000, radius-0
and radius-1000 runs with `prepare-window`, fit Bayesian H1 for reciprocal
cohorts, integrate baseline and ppBMI H5 context, then prepare H6/H7.
The model specification is `docs/H1_BAYES_MODELS.md`.

Scripts: `/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/scripts`

Data: `/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/data/<cohort>`

| Cohort | Samples |
|---|---|
| maize_B73_Ki11 | Four PE endosperm pools, 20 DAP |
| maize_B73_W22 | Six PE endosperm pools, 14 DAP |
| mouse_Babr_CAST | Five SE extraembryonic ectoderm units, E7.5 |
| mouse_FVB_CAST_E12 | Four PE female placenta units, E12.5 |
| fish_Pprolifica_AB | *Poeciliopsis prolifica* Line A × Line B; reciprocal placental RNA |
| GUSTO | Local placenta RNA, phased transcriptomic genotypes and precomputed chromosome parental-origin labels |

The fish profile includes the GEO/SRA study, Dryad fixed-SNP resources, RefSeq
genome and independent Line A PacBio HiFi run. The accession-to-cross/tissue
audit and the Line A/Line B SNP table are documented in `docs/fish_*.md`; see the
walkthrough before attempting submission.

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
| Additional analyses (`STAGE=additional`) | Separate TSS/donor questions, relaxed reciprocal callability, optional TE tests and GUSTO gene POE in `additional_v1/` |
| [Baseline NB](docs/rerun_baseline_poe.md) | All-SNP gene/transcript POE, TE-TSS comparisons, human context, results archive |
| [Export](07_exports/README.md) | Write optional human aggregates |

Stage commands use `bin/submit_pipeline.sh "$CONFIG" STAGE --submit`.
Omitting `--submit` prints a plan. Results include completion status, QC tables,
pilot metrics, and model estimates.

For GUSTO recovery after a completed shared reference, use
`reference --resume --sample J1001 --submit` to validate one child, then
`reference --resume --submit` to submit only incomplete child/index jobs.
See [GUSTO in the walkthrough](docs/Placenta_EvoDevo_step_by_step.md#h4-build-the-catalog-and-child-indices).

`STAGE=additional` requires a completed primary `models` receipt and creates
`data/<cohort>/additional_v1/` without changing primary results or receipts.
Use `bin/Rscript_sif poe.R "$CONFIG" additional-status` to verify completion.

Settings are described in the [schema](config/SCHEMA.md).
For opt-in hypothesis #1 (same-first-donor TE/non-TE TSS POE magnitude) and #5
(human cis-QTL/chromatin context), see the [regulatory modules runbook](docs/REGULATORY_MODULES.md).
These use separate `additional_v2` outputs, species-specific resource retrieval,
and separately scheduled fish repeat discovery/masking.
## Results archive

To hand results to a reporting or research workflow, bundle every cohort's
primary models, QC, `additional_v1`, hypothesis #1 (`te-tss`, `te-genes`) and hypothesis #5
(`human-context`) outputs with the method docs:

```bash
mkdir -p "$DATA_ROOT/exports"
bin/Rscript_sif bin/export_results.R --out "$DATA_ROOT/exports/results_$(date +%Y%m%d).tar.gz" \
  --force config/*.private.yaml
```

Without `--force` it stops if any applicable part is incomplete; with `--force`
it skips those parts (e.g. fish `te-tss` while repeat discovery runs) and lists
them in the archive's `README.md` and `AVAILABILITY.tsv`. Human tables are
aggregate-only: GUSTO tables with sample/family/parent identifiers are withheld.

See [operations](docs/OPERATIONS.md) for logs/status and
[attribution](ATTRIBUTION.md) for software/data credits.
