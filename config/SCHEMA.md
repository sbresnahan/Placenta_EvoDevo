# Configuration and table schema

Each command uses one YAML with `schema_version: 2`.
The six example configurations provide the paths and settings.
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
| `retrieval` | `method: https` (default) or `method: sra` for public Seadragon cohorts; `sra_threads` (default 1), `max_parallel_jobs` (default 2), retries and connection settings |
| `lsf` | Queue, numeric GB units, cores, memory/reservation and wall time per stage |
| `models` | Public depth and minimum callable unit thresholds; human family/orientation/bootstrap settings and covariate names |
| `fragments` | SE fragment mean/SD and measured/assumed origin |
| `exports` | `enabled` and positive integer `min_cell` |

Manual approval fields and pilot runtime/storage caps are retired.
Old fields in existing private YAMLs are ignored.

Stage readers validate configuration values and relevant cohort tables.

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

GUSTO's `gusto_precomputed.enabled: true` adapter uses the existing phased
transcriptomic VCF (`vcf`), chromosome Left/Right parent-of-origin matrix
(`origin_matrix`), and infant/RNA ID map (`fastq_id_map`). Set
`phase_origin.enabled: false`. It copies the original ESPRESSO transcript FASTA
from `reference.transcript_source` and parses the matching ESPRESSO GTF at
`paths.annotation`. The underscore in VCF column `010-04002_B010-04002`
separates an upstream label from the infant identity `B010-04002`;
the FASTQ ID table links that child to `J1010` (for example). Genomic hg38
coordinates are parsed from VCF `ID`, while `CHROM` and `POS` refer to
transcript coordinates. The matrix holds one origin call per chromosome and
Left/Right haplotype, not one row per SNP. Biallelic SNPs and indels enter
the personalized transcriptome/index; only SNPs are mapped to genomic loci
for parental count orientation. Only defined origins at phased heterozygous
SNP genotypes yield parental counts. Missing allele origin is `NA`.
The standard human `typed_vcf`/`personalization_vcf` fields remain required
in schema 2 for legacy compatibility but are not accessed by this adapter.
Unspecified human FASTQ read counts are allowed, with pilot runtime projection
reported as `NA`; all non-human pilot read totals remain mandatory.

## Tables

| Table | Required content |
|---|---|
| Sources | Unique run ID, linked library ID, layout, mate URLs and byte sizes, read count |
| Libraries | Unique sample/unit IDs, layout and explicit strand mode; human tables add `child_id` |
| Units | Unique unit ID, biological `cluster_id`, maternal/paternal identities, stage, tissue |
| Resources | ID, HTTPS URL, destination/partial/response paths, byte size (except CyVerse VCFs), required flag, build; optional decompression destination |
| Marker sources | Unique source ID, path, build |
| Trios | Legacy: child, mother, father, family IDs. Precomputed GUSTO: child and family IDs only |
| Covariates | Sample ID and configured covariate columns |

Human source rows add actual `local_1`/`local_2` FASTQ paths.
SRA mode uses existing SRR/ERR/DRR `run_id` accessions without renaming them.
It creates FASTQs from NCBI SRA archives, so the manifest's ENA gzip byte
counts **do not apply** to those generated files. The host writes a
per-run `.sra.tsv` size receipt and the R finish gate verifies recorded sizes,
mate consistency and manifest fragment counts. Existing ENA files are
accepted when their original manifest byte sizes match. Legacy MD5 columns
may remain in manifests but are ignored.
PE read counts include both mates. Multiple runs link to one primary library.
Missing or repeated cluster IDs do not gate public-cross model fitting.
Public models require `models.min_callable_units` ≥3 (default 3), a single
preselected marker with depth ≥`models.min_depth` (default 10 in examples)
in those units, and observations from both reciprocal cross directions.
The gene event or two promoter events need not be callable in every unit.
Marker selection prioritizes callable unit count, then minimum depth, before
looking at parental allelic imbalance. P-values are exploratory unit-level
calculations: repeated embryos from one litter, seeds from one ear, or pooled
samples may violate independence and make standard errors anti-conservative.
Almost-zero residual root-mean-square (≤ square-root machine epsilon times
the larger of 1 and the largest absolute response) is flagged as
`degenerate_residual_variance`: the effect remains descriptive, while SE,
confidence limits, p-value, and q-value are `NA`. This prevents extreme
p-values driven by tied allele ratios in small samples. Four callable units
are recommended for the primary sensitivity run; keep three-unit results
exploratory even when the variance test passes. For the sensitivity run, set
`models.min_callable_units: 4`.
The cluster provenance field remains in the unit manifest for interpretation.
The human trio gates remain unchanged (≥20 depth, ≥10 families, etc.).

Prepared marker tables have `chrom`, `genomic_pos`, `ref`, `alt`,
`focal_index`, and `alternate_index`: one-based positions, single-base alleles,
and opposite homozygous indices (0=REF, 1=ALT).

Legacy manual sign-off columns in existing tables are ignored.
