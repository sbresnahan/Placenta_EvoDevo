# Placenta_EvoDevo

Parent-of-origin expression analysis of maize endosperm, mouse extraembryonic
tissues, placental fish, and human GUSTO placenta using Ornaments.

## Start here

Follow the [step-by-step walkthrough](docs/Placenta_EvoDevo_step_by_step.md):
set project paths, select a cohort, check tools, and submit the first stage.

The all-SNP baseline layer (gene and transcript POE from the original count
batches) is built with the [baseline workflow](docs/BASELINE_NB_WORKFLOW.md);
outputs are described in the [baseline schema](docs/BASELINE_NB_SCHEMA.md).

From a completed baseline layer, the [analysis runbook](docs/Seadragon_analysis_runbook.md)
creates independent radius-3000, radius-0 and radius-1000 promoter-window runs
with `prepare-window`, fits Bayesian H1 for reciprocal cohorts and descriptive
H1 for GUSTO, integrates baseline and ppBMI H5 context, and prepares H6/H7.
Model specifications: [H1](docs/H1_BAYES_MODELS.md), [H5](docs/H5_ppBMI_CONTEXT.md),
[H6/H7](docs/H6_H7_MODELS.md).

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
| [Repeat resources](docs/REGULATORY_RESOURCES.md) | Repeat tracks, fish repeat discovery/masking and human context inputs |
| [Baseline NB](docs/BASELINE_NB_WORKFLOW.md) | All-SNP gene/transcript POE layer and its results archive |
| [Promoter-window analyses](docs/Seadragon_analysis_runbook.md) | H1, between-gene comparison, H5 context and H6/H7 per window radius |
| [Export](07_exports/README.md) | Write optional human aggregates |

Stage commands use `bin/submit_pipeline.sh "$CONFIG" STAGE --submit`.
Omitting `--submit` prints a plan. Results include completion status, QC tables,
pilot metrics, and model estimates.

For GUSTO child/index jobs after a completed shared reference, use
`reference --resume --sample J1001 --submit` to validate one child, then
`reference --resume --submit` to submit only incomplete child/index jobs.
See [GUSTO in the walkthrough](docs/Placenta_EvoDevo_step_by_step.md#h4-build-the-catalog-and-child-indices).

Settings are described in the [schema](config/SCHEMA.md).

## Results archives

The promoter-window and baseline analyses have their own archives:
`bin/export_results.R --baseline-nb` (POE layer, H1, between-gene, H5) and
`--exposure-nb` (H6/H7); see the runbook. To bundle every cohort's primary
promoter models, QC and manifests with the method docs:

```bash
mkdir -p "$DATA_ROOT/exports"
bin/Rscript_sif bin/export_results.R --out "$DATA_ROOT/exports/results_$(date +%Y%m%d).tar.gz" \
  --force config/*.private.yaml
```

Without `--force` it stops if any applicable part is incomplete; with `--force`
it skips those parts and lists them in the archive's `README.md` and
`AVAILABILITY.tsv`. Human tables are
aggregate-only: GUSTO tables with sample/family/parent identifiers are withheld.

See [operations](docs/OPERATIONS.md) for logs/status and
[attribution](ATTRIBUTION.md) for software/data credits.
