# Configuration and table schema

Each command uses one YAML with `schema_version: 2`.
The five example configurations provide the paths and settings.
Relative paths resolve against the YAML's directory.

## Configuration

| Section | Contents |
|---|---|
| `profile`, `species`, `cohort_id` | Public-cross or human-trio analysis and cohort identity |
| `paths`, `templates` | Inputs, outputs, operational tables, repeated library/run paths |
| `tools` | Host wrappers, R bridge, utility paths, R library |
| `execution` | Backend, R image, binds, site file |
| `ornaments` | `semantics: ref_alt`, native threads |
| `reference` | Build, genome/annotation, contig aliases |
| `parents`, `dosage` | Focal/alternate parents and maternal/paternal dosage |
| `markers` | Adapter, source/sample identities, quality filters, proxy description |
| `pilot` | Strand/throughput fragment counts and diagnostic thresholds |
| `retrieval` | `method: https` (default) or `method: sra` for public Seadragon cohorts; `sra_threads` (default 1), retries and connection settings |
| `lsf` | Queue, numeric GB units, cores, memory/reservation and wall time per stage |
| `models` | Depth/family/orientation/bootstrap settings and covariate names |
| `fragments` | SE fragment mean/SD and measured/assumed origin |
| `exports` | `enabled` and positive integer `min_cell` |

Manual approval fields and pilot runtime/storage caps are retired.
Old fields in existing private YAMLs are ignored.

The tool-only setup reader resolves tool, site, image, and inventory paths.
Stage readers also validate configuration values and relevant cohort tables.

## Paths

The data root is
`/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/data/<cohort>/`.
Working YAMLs belong in `scripts/config/`.

Allowed template tokens are `{sample}`, `{run}`, `{mate}`, `{mode}`, `{strand}`,
`{stage}`, and `{job}`. Run/sample/unit keys use letters, numbers, underscores,
periods, and hyphens. `${POE_SITE_ROOT}` is the optional environment substitution.

FASTA indices use adjacent `.fai` paths.
Ornaments writes `allele_counts.txt`, `thetas.txt`, and `tpms.txt` within
each configured `quant_dir`.

GUSTO uses `reference.genome_source` for the shared FASTA and `paths.annotation`
for the shared HPLRv2 GTF. Derived references stay in the cohort directory.

## Tables

| Table | Required content |
|---|---|
| Sources | Unique run ID, linked library ID, layout, mate URL/MD5/byte fields, read count |
| Libraries | Unique sample/unit IDs, layout and explicit strand mode; human tables add `child_id` |
| Units | Unique unit ID, biological `cluster_id`, maternal/paternal identities, stage, tissue |
| Resources | ID, HTTPS URL, destination/partial/response paths, bytes/MD5, required flag, build; optional decompression destination |
| Marker sources | Unique source ID, path, build |
| Trios | Child, mother, father, family IDs |
| Covariates | Sample ID and configured covariate columns |

Human source rows add actual `local_1`/`local_2` FASTQ paths.
SRA mode uses existing SRR/ERR/DRR `run_id` accessions without renaming them.
It creates FASTQs from NCBI SRA archives, so the manifest's ENA gzip MD5 and
byte counts **do not apply** to those generated files. The host writes a
per-run `.sra.tsv` checksum receipt and the R finish gate verifies generated
checksums, mate consistency and manifest fragment counts. Existing ENA files
are accepted only when their original manifest MD5/size match.
PE read counts include both mates. Multiple runs link to one primary library.
Missing cluster IDs retain descriptive public results.

Prepared marker tables have `chrom`, `genomic_pos`, `ref`, `alt`,
`focal_index`, and `alternate_index`: one-based positions, single-base alleles,
and opposite homozygous indices (0=REF, 1=ALT).

Legacy manual sign-off columns in existing tables are ignored.
