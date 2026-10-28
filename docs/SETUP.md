# Setup reference

Use [the walkthrough](WALKTHROUGH.md) for the execution order. This page explains
the settings and checks needed to make those commands runnable.

## Project directories

| Purpose | Location |
|---|---|
| Script root | `/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/scripts` |
| Data root | `/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/data` |
| Cohort files | `data/<cohort>/` with separate manifests, raw, reference, private, work, scratch, logs and receipts |
| Private settings | `scripts/config/<profile>.private.yaml` |

The operational data paths are explicit in each YAML. Changing `DATA_ROOT` in
your shell does not rewrite them. If you edit a root, update its related
`paths`, `templates` and resource/marker-source table paths consistently.
The shipped `manifests/` folder contains templates, not the operational copies.
Do not copy templates over edited cohort tables.

Use a POSIX filesystem for work, temporary files, locks and atomic receipts.
Verify write permissions and reserve raw/reference/index/scratch headroom.
Creating directories or accepting the example budget does not establish
sufficient disk or memory.

## Software environments

Native host tools are activated by repository wrappers. R is isolated in
Singularity; native tools and `bsub` must not be assumed to work inside it.

| Setting | Default |
|---|---|
| Module initialization | `/etc/profile.d/modules.sh` |
| Conda executable | `/risapps/rhel8/miniforge3/24.5.0-0/bin/conda` |
| samtools environment | `samtools-1.16.1` |
| bcftools environment | `bcftools-1.16` |
| Ornaments activation | `/rsrch5/home/epi/bhattacharya_lab/software/ornaments/bin/activate` |
| Installed commands after activation | `ornaments`, `create_personalized_transcriptome` |
| R image | `/risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif` |
| R library | `/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1` |
| R binds | `/rsrch5`, `/rsrch9`, `/tmp` |

Inspect trusted `config/seadragon.sh`; do not replace it with arbitrary
participant-supplied shell code. `POE_SITE_FILE` and `execution.site_file` must
identify the same file. An R-image override must match `execution.image`.
The bridge's actual binds must cover all inputs/outputs.

Required R packages: **data.table, digest, yaml, curl, GenomicRanges, IRanges,
S4Vectors, Rsamtools and Biostrings**. R must be at least 4.3.1. The preflight
checks packages through the R bridge. Missing/incompatible packages need a
reviewed site fix; the project does not modify the shared library automatically.
Host utilities include Bash, GNU time, gzip, curl, git and the site scheduler/
Singularity executables.

## Configuration checklist

Copy the selected example beside itself in `config/`. All five examples use
the same required schema; retain keys even when their branch does not use
them. See [the schema](../config/SCHEMA.md) for table headers and path rules.

| Exact key or table field | What to enter/check | Timing |
|---|---|---|
| `paths.*`, `templates.*` | Correct absolute cohort destinations, no cross-cohort collisions | Before any stage |
| `execution.site_file`, `execution.image` | Same trusted site file/image used by shell wrappers | Before preflight |
| `ornaments.binary_sha256` | Reviewed digest of the installed binary, not of its wrapper | Before execution |
| `ornaments.installation_provenance`, `ornaments.installation_provenance_sha256` | Reviewed identity record and its digest if source identity needs it | Before execution |
| `execution.posix_work_approved` | True only after checking work/scratch/receipt storage | Before submission |
| `lsf.memory_units_verified`, `lsf.memory_per_job_verified` | Actual site interpretation of limit and reservation | Before submission |
| `lsf.execution_environment_verified` | Worker access to image, library, native tools and input/output storage | Before submission |
| `lsf.stages.*` | Queue, `cpus`, `memory_value`, `reservation_value`, `wall_time` | Before reference/pilot; review against measurements |
| `lsf.resource_limits_approved` | Reviewed total resource budget, not just one job's RAM | Before pilot/scale acceptance |
| `retrieval.internet_host_approved` | Approved internet/transfer execution location | Public retrieval only |
| Resource rows | Correct HTTPS source, expected byte size/MD5, build and destinations | Before acquiring that resource |
| `reference.genome_sha256`, `reference.annotation_sha256` | Verified exact uncompressed input bytes | After files available, before reference |
| `reference.annotation_approved` | Correct tissue/build/annotation and scientific review | Before reference |
| Marker-source rows, `markers.provenance_approved` | Correct source path/build, verified SHA256 and actual parent/proxy provenance | Before public markers |
| `markers.proxy*` | Explicit applicable proxy, approval and reason; no silent strain substitution | Before public markers |
| Library `strand_mode`, `strand_verified` | Actual protocol review, not guessed from labels | Before pilot |
| Unit `cluster_id`, `independence_verified` | True independent unit/cluster definition | Before inference; review before pilot |
| `dosage.*` | Maize normal endosperm 2:1; mouse/human 1:1 | Before processing |
| `fragments.*` | Single explicitly labelled SE assumption; PE does not use it | Before reference/pilot |
| `human.genotype_qc_approved`, library `identity_verified` | Independent real RNA/DNA/pedigree QC and exact trio linkage | Before human preflight |
| `models.covariates`, `human.covariates_approved` | Approved covariate list, exact linked rows and centering | Before human models; configure before pilot |
| `exports.*` | Disabled during analysis; explicit local disclosure review separately | After completed human analysis |

Approval fields are attestations, not instructions to fabricate evidence.
Empty/unresolved input contracts must be completed; setting a flag does not
create a reference, genotype table, independent replicate or identity link.

Fingerprinting is conservative: even a later approval/resource-setting edit
can invalidate reference/pilot/scale receipts. Complete settings before reference
preparation. If pilot results require changes, rerun affected stages rather
than editing receipts or reusing an accepted report.

### Calculate file checksums

For a real file that you have independently verified against its intended
source/build, obtain its exact SHA256:

```bash
read -r -p "Absolute path to verified file: " VERIFIED_FILE
[[ -f "$VERIFIED_FILE" && "$VERIFIED_FILE" == /* ]] || { printf 'Provide an existing absolute path\n' >&2; exit 1; }
sha256sum -- "$VERIFIED_FILE"
```

Use the **uncompressed** genome and annotation for the YAML reference hashes.
Use the exact file named by a marker-source row for its SHA256 (including
compression if the source path ends in `.gz`). MD5/byte-size acquisition
contracts and SHA256 analysis identities serve different checks.

## Installed-tool identity

With `SCRIPTS_DIR`, `POE_SITE_FILE` and absolute `CONFIG` set:

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
```

The first command checks metadata. The second inventories actual host binaries
and checks R dependencies separately. It does not install tools, download
data or submit jobs. The pending inventory path is printed before acceptance.
Review that file privately if a gate fails.

The supported Ornaments source identity is:

```text
2c4b1a2c18616f9de1cb307749688daaf3313d17
```

The inventory records real executable paths, versions and SHA256 values, source
identity when available, site/image identity and a Python dependency digest.
`tools.ornaments_source` is an installation prefix, not proof of source identity.
Set `POE_ORNAMENTS_SOURCE` only to a real clean source checkout that genuinely
establishes installation provenance. A matching checkout alone does not prove
that an arbitrary binary was built from it; review installation evidence.

If source identity is unavailable, prepare a reviewed private installation
record at `ornaments.installation_provenance`. Its required YAML structure is:

```yaml
reviewed: true
commit: 2c4b1a2c18616f9de1cb307749688daaf3313d17
binary_sha256: REVIEWED_INSTALLED_ORNAMENTS_SHA256
personalizer_sha256: REVIEWED_INSTALLED_CONSOLE_ENTRYPOINT_SHA256
evidence: LOCAL_BUILD_OR_INSTALL_EVIDENCE_LINKING_THESE_FILES_TO_THE_SOURCE
```

Uppercase values above are **placeholders**, not usable digests/evidence.
Complete them only from actual reviewed files and build/install records.
Enter the installed Ornaments digest in `ornaments.binary_sha256` and the
installation-record file's exact digest in
`ornaments.installation_provenance_sha256`. Re-run preflight.

**Success:** an accepted inventory at `paths.native_inventory` and valid metadata.
**Failure:** unresolved/mismatched identity, missing library/image/package or
bad manifests must be corrected locally. Never waive the pin or automatically
rebuild the installed software.

Metadata preflight does not certify that all reference/read inputs exist,
that LSF policy is correct or that models will be eligible. GUSTO's empty
schemas cannot pass until real private linkage/QC is supplied.

## Resource review

The YAML's queue/CPU/RAM/time values are initial templates, not measurements.
Confirm the site's memory limit and reservation units; do not translate values
from another cluster by guesswork. Native thread count cannot exceed allocated
CPUs. Review public input totals from manifests and actual local GUSTO inventory.

Use the bounded pilot's measured RSS/time and actual index sizes to choose
allocations and review total runtime/storage. Projection covers quantification,
not all downloading, indexing, startup or queue delays. Endpoint work size
is not a disk high-water mark. The R image is hashed during host inventory;
its I/O overhead also needs site observation.

## Analysis rules

These are interpretation/eligibility contracts, not optional setup shortcuts:

- Promoter groups use a fixed TSS anchor, a primary 50-bp window and compatible
  first donor; no transitive chaining. The 25/100-bp catalogs are structural
  sensitivity views, not evidence of tested expression-model robustness.
- Expected allocations are combined at the same event/genomic locus, not
  summed across SNPs as independent molecules. Missing parental evidence
  remains NA.
- One sentinel per entity is selected by coverage across every configured
  primary unit with depth at least 20, not by imbalance. Promoter-pair selection
  is likewise coverage-based.
- Public gene outcomes are `log2((M + 0.5)/(P + 0.5))` minus
  `log2(maternal_dose/paternal_dose)`. Promoter contrasts are within-unit
  differences, so this dosage offset cancels.
- Public inference uses an unweighted intercept plus reciprocal-cross
  orientation. Residual degrees of freedom are independent-unit n minus
  model rank. One direction, repeated/unresolved clusters, missing units,
  insufficient coverage or nonpositive residual df blocks POE inference.
- BH adjustment is separate by cohort, stage, species and endpoint. Do not
  combine results as extra replicates or silently change the testing family.
- Human inference uses the paired natural-log promoter contrast, cis-orientation
  terms and approved centered covariates. Only exact equal/opposite cis
  redundancy is removed; otherwise full rank is required.
- Human eligibility requires depth at least 20/event, ten paired families,
  three maternal-REF and three maternal-ALT orientations at both sentinels and
  at least five residual df. HC3 t inference and 2,000 family bootstraps
  use seed 1147 and require at least 90% valid bootstrap fits.
- Human `q <= 0.05` and absolute contrast at least `log(2)` flag conditional
  discordance, not definitive imprinting. Uncertainty is conditional on
  quantification and marker assumptions.

Maternal seed-coat/decidua/blood admixture, cis effects, genotype proxies and
developmental composition remain alternative explanations. Human whole-gene
analysis/masking and TE comparisons are outside the implemented workflow.
