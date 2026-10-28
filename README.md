# Placenta_promoter_POE

Raw-read parent-of-origin expression analysis for maize endosperm, mouse
extraembryonic tissues and human GUSTO placenta. The workflow prepares
references and parental markers, quantifies RNA with Ornaments, assigns
maternal/paternal expression and fits eligible gene/promoter models.

## Start here

**Follow [the Placenta_EvoDevo step-by-step workflow](docs/Placenta_EvoDevo_step_by_step.md).**
This is the primary, self-contained operating guide. Complete shared
setup, choose one species and one cohort, then follow that section through
results. Commands run from:

```text
/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/scripts
```

Cohort inputs, scratch, logs and outputs go under:

```text
/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/data/<cohort>
```

Installed software and the configured shared human reference inputs remain
at their source locations. The package does not move existing data.

| Species | Profile | Scope |
|---|---|---|
| Maize | `maize_B73_Ki11` | 20-DAP endosperm; four biological pools, eight technical runs |
| Maize | `maize_B73_W22` | 14-DAP endosperm; six ear-pool libraries |
| Mouse | `mouse_Babr_CAST` | E7.5 extraembryonic ectoderm; five SE libraries; explicit B6J marker proxy |
| Mouse | `mouse_FVB_CAST_E12` | E12.5 placenta; four PE libraries |
| Human | `GUSTO` | Private placenta RNA and independently verified typed trios |

Supplied: scripts, example configurations, public run tables and empty human
schemas. Required locally: reads, references, variants, the human annotation,
verified input contracts and applicable approvals.

## Project stages

| Module | Task | Recommended action |
|---|---|---|
| [00_setup](00_setup/README.md) | Settings, R packages and tool checks | `preflight.sh` |
| [01_retrieval](01_retrieval/README.md) | Verified public acquisition | `retrieve` |
| [02_reference](02_reference/README.md) | Transcript/promoter catalog | `reference` |
| [03_genotypes](03_genotypes/README.md) | Parental markers or typed child inputs | Included in `reference` |
| [04_ornaments](04_ornaments/README.md) | Indices, pilots and quantification | `reference`, `pilot`, `scale` |
| [05_parental_expression](05_parental_expression/README.md) | Parental assignment | Included in `scale` |
| [06_promoter_models](06_promoter_models/README.md) | Callability and inference | `models` |
| [07_exports](07_exports/README.md) | Private human aggregate review | Separate explicit actions |

`bin/submit_pipeline.sh` is the host submission interface; `bin/Rscript_sif`
runs R in Singularity. Native tools and `bsub` run on the host. Commands without
`--submit` are dry-runs. A completed pilot never starts scale automatically.
Do not submit individual modules in addition to the same driver graph.

## Reference pages

- [Setup checklist](docs/SETUP.md) and [configuration schema](config/SCHEMA.md)
- [Monitoring, manual submission and troubleshooting](docs/OPERATIONS.md)
- [Validation and live-site checks](docs/VALIDATION.md)
- [Code inventory](MANIFEST.md) and [attribution](ATTRIBUTION.md)
- [Scientific resource review](evidence/review.md)

Local tests do not establish actual Seadragon/cohort performance. Expected
allocations and promoter coordinates alone do not prove imprinting or an
intragenomic-conflict mechanism.
