# Configuration and manifest schema

Each invocation requires one explicit YAML using `schema_version: 2`.
The five example files define the complete structure. Paths are not merged
from another configuration. Values are never evaluated as shell/R code.

## Paths and settings

Relative YAML paths **and path columns in tables** resolve against the YAML's
directory. Keep private copies in `scripts/config/` when using the supplied
relative wrappers/site file. Only `${POE_SITE_ROOT}` is allowlisted for optional
environment substitution; it must actually be set if used.

Project data paths point to
`/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/data/<cohort>/`.
Operational tables belong in its `manifests/` directory; copy supplied templates
there using the [step-by-step workflow](../docs/Placenta_EvoDevo_step_by_step.md) before preflight.
Human shared input references and installed software are explicit exceptions.

| Section | Required meaning |
|---|---|
| `schema_version`, `profile`, `species`, `cohort_id` | Schema 2; `public_cross` maize/mouse or `human_trio` human; safe unique cohort identifier |
| `paths` | Single-file/root inputs, outputs, operational manifests, receipts and snapshot |
| `templates` | Repeated per-run/library/mode/strand/job paths |
| `tools` | Project host wrappers, R bridge, utilities and installed R-library path |
| `execution` | `backend: seadragon`, image/binds, trusted site file and POSIX approval |
| `ornaments` | Required REF/ALT output semantics and native thread count |
| `reference` | Build, compatible genome/annotation files, annotation approval and explicit contig spelling map |
| `parents`, `dosage` | Actual focal/alternate identities and reviewed tissue/ploidy dosage |
| `markers` | `panzea`, `mgp` or `prepared`, exact samples, quality/source and proxy checks |
| `pilot`, `retrieval` | Bounded prefixes/diagnostics, budgets, retry/timeouts and internet-host approval |
| `lsf` | Verified memory/environment semantics and per-stage queue/CPU/limit/reservation/time |
| `models`, `human` | Depth/family/orientation/bootstrap floors, covariates and genotype/identity approvals |
| `fragments` | One YAML mean/SD/origin setting applied to SE libraries only |
| `exports` | Disabled initially; local policy, cell threshold, review token and suppression review |

Keep all base path/template/tool keys even for unused branches. Required
Seadragon keys include `paths.native_inventory`, `catalog_receipt`,
`marker_receipt`, `job_registry` and templates `child_receipt`,
`library_pilot_report`, `library_pilot_receipt`, `library_metrics`.
`provenance_hints` supplies context only; it does not establish linkage.
Settings are recorded for the operator; they are not compared with earlier settings to decide whether results are current.

## Software and reference identity

Infrastructure defaults live in `seadragon.sh` and `bin/Rscript_sif`.
`execution.site_file` must match `POE_SITE_FILE`; image overrides must match
`execution.image`. The installed commands must be available in the configured native environments. Inventories record command paths, available versions, site/image paths, and Python package information. They are not compared with an earlier installation.

Keep `ornaments.semantics: ref_alt`; the count reader still requires the eight-column REF/ALT output contract.

For GUSTO, `reference.genome_source` is read-only. The wrapper copies it to
`paths.genome` before indexing. `paths.annotation` names the exact shared
HPLRv2 GTF; normalized GTF/transcripts/indices stay in the cohort data directory.

`paths.genome_fai` and `paths.transcripts_fai` must be the adjacent `.fai`
files required by samtools/Rsamtools. Quantifier outputs must declare
`allele_counts.txt`, `thetas.txt` and `tpms.txt` directly inside `quant_dir`.
Do not rename those upstream output contracts.

Allowed template tokens are `{sample}`, `{run}`, `{mate}`, `{mode}`, `{strand}`,
`{stage}` and `{job}`. Substituted values contain only letters, digits, dots,
underscores and hyphens; unknown/unresolved tokens fail validation.

## Manifest schemas

Supplied `manifests/*.tsv` provide headers and public selections. Edit the
operational copies, not source templates. All tables must keep exact unique
keys and linkage; a sample label is not automatically a biological replicate.

- **Sources:** unique `run_id`, `sample_id`, layout, HTTPS URL/MD5/byte columns
  per mate and `read_count`/`base_count` for projection. One row is one physical
  run. PE `read_count` is twice fragment count. Human sources retain common
  headers and add actual `local_1/local_2`; use real read counts.
- **Libraries:** unique `sample_id` and `unit_id`, one primary library per
  unit, SE/PE, `fr`/`rf`/`unstranded` and `strand_verified`. SE parameters come
  solely from YAML. Human libraries require `child_id`, `identity_verified`.
- **Units:** unique `unit_id`, `cluster_id`, `independence_verified`, actual
  maternal/paternal identities, tissue and stage. Unresolved independence
  blocks inference; one public configuration cannot mix stages/tissues.
- **Resources:** unique `resource_id`, HTTPS URL, destination/partial/response,
  expected bytes/MD5, required flag, build and optional `decompress_to`.
  Unresolved required contracts cannot be downloaded.
- **Marker sources:** unique `source_id`, path, and build. Prepared
  data have `chrom`, `genomic_pos`, `ref`, `alt`, `focal_index`,
  `alternate_index`; indices are opposite homozygous 0/1 states.
- **Trios:** exact unique child, mother, father and family linkage; no conflicting
  pedigree identities. One offspring per family is selected independently of
  imbalance for human inference.
- **Covariates:** unique `sample_id` with the approved model covariates and
  exact linked rows; no missing covariates in inference.

Empty private schemas contain no participants and intentionally fail until
populated with real approved inputs. Other cohorts in the evidence catalog are
not automatically included in the four configured public selections.
