# Project walkthrough

Complete **shared setup once**, then follow **one species and one profile at
a time**. Do not paste the whole guide into an unattended shell: wait and check
completion at each stop.

- [Shared setup](#shared-setup)
- [Maize](#maize)
- [Mouse](#mouse)
- [Human GUSTO](#human-gusto)
- [Reading the results](#reading-the-results)

## Shared setup

### 1. Extract the scripts into the project

The ZIP contains a top-level `scripts/` folder. Put the ZIP in an approved
location, then enter its actual absolute location when prompted:

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
mkdir -p "$PROJECT"
read -r -p "Absolute path to Placenta_promoter_POE.zip: " PACKAGE
[[ -f "$PACKAGE" && "$PACKAGE" == /* ]] || { printf 'Provide an existing absolute ZIP path\n' >&2; exit 1; }
[[ ! -e "$PROJECT/scripts" ]] || { printf 'scripts already exists; inspect it before extracting\n' >&2; exit 1; }
unzip "$PACKAGE" -d "$PROJECT"
cd "$PROJECT/scripts"
export SCRIPTS_DIR="$PWD"
export POE_SITE_FILE="$PWD/config/seadragon.sh"
export DATA_ROOT="$PROJECT/data"
umask 077
mkdir -p "$DATA_ROOT"
chmod +x bin/* bin/native/* 00_run_pipeline.sh 00_setup/*.sh */*.lsf 04_ornaments/*.sh
```

Standard `unzip` preserves executable permissions; the `chmod` command also
covers extractors that do not. Do not overwrite an active scripts directory.
The `scripts/` directory itself is the project root—there is no extra
`Placenta_promoter_POE/` level.

In a new shell, restore the environment:

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
cd "$PROJECT/scripts"
export SCRIPTS_DIR="$PWD"
export POE_SITE_FILE="$PWD/config/seadragon.sh"
export DATA_ROOT="$PROJECT/data"
umask 077
```

All remaining commands run from this directory.

### 2. Check the software and storage

Inspect `config/seadragon.sh`:

| Component | Configured location or environment |
|---|---|
| Conda | `/risapps/rhel8/miniforge3/24.5.0-0/bin/conda` |
| samtools | `samtools-1.16.1` |
| bcftools | `bcftools-1.16` |
| Ornaments activation | `/rsrch5/home/epi/bhattacharya_lab/software/ornaments/bin/activate` |
| R image | `/risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif` |
| R library | `/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1` |

The wrappers activate native tools separately. R uses `bin/Rscript_sif`; native
tools and `bsub` run on the host. Do not install another Ornaments environment
or activate all environments together. The R image needs `/rsrch5`, `/rsrch9`
and `/tmp` visible.

All cohort files use this layout:

```text
Placenta_EvoDevo/
  scripts/                 scripts, configuration examples, table templates
  data/
    <cohort>/
      manifests/           operational run/library/unit/resource tables
      raw/                 raw FASTQs
      reference/           downloaded/prepared references and shared index
      private/             private inputs, child indices, plans and snapshots
      work/                quantification and parental allocations
      scratch/             temporary/prefix/native execution files
      logs/                tool and LSF logs
      receipts/            validated completion records
      model_*.tsv          final configured model tables
```

Human shared reference inputs are exceptions: keep their exact source paths,
specified in the GUSTO section. Software stays at its installed locations.
Check permissions, POSIX behavior, storage capacity and queue access; the guide
does not establish that the cluster paths are available.

### 3. Prepare settings without bypassing checks

Each species section creates a private YAML beside its example in `config/`
and copies table templates into `data/<cohort>/manifests/`. Copies do not
overwrite existing edited files. Do not keep participant payloads in `scripts/`.

| Edit/review | What it controls |
|---|---|
| Private YAML | References/hashes, installed-tool identity, resources, approvals and result paths |
| Sources table | Physical runs, FASTQ paths, expected checksums/sizes |
| Libraries table | Primary libraries, SE/PE layout and independently reviewed strand |
| Units table | Actual parents, tissue/stage and biological independence |
| Resources table | Approved public asset URLs, sizes/checksums and destinations |
| Marker-source table | Parental source paths, compatible build and exact SHA256 |
| GUSTO trio/covariate tables | Private linkage and approved model covariates |

Use the exact-key [setup checklist](SETUP.md#configuration-checklist).
Paths in table columns resolve against the YAML directory, not the table
directory; supplied operational destinations are absolute.

Before retrieval, complete tool, storage, internet-host and LSF checks. After
retrieval, enter verified uncompressed genome/annotation hashes and marker
source hashes before reference preparation. **Do not turn on every approval
flag as a shortcut.**

### 4. Recognize submission versus completion

- Without `--submit`, stage commands print a plan and may write planning/status
  files, but do not submit jobs.
- With `--submit`, the driver returns after graph submission, not completion.
- Inspect returned numeric IDs and logs. `bjobs` lists active jobs; `bjobs -a`
  includes recent jobs. For a specific job, run `bjobs -l` followed by its
  actual numeric ID. A job disappearing from the default list is not proof
  of success.
- `bin/Rscript_sif poe.R "$CONFIG" status` validates files/receipts, not live
  scheduler state. Downstream stages can correctly show blocked before they
  are run.
- Wait for every job in a stage, including its collector, before continuing.
  Do not submit modules manually as well as through the same driver graph.

Once reference work starts, keep settings fixed. Changes to settings, resources,
manifests, inputs or tools can invalidate affected work. If a pilot reveals
that changes are necessary, correct them and repeat affected prerequisites and
the pilot before accepting a report. See [operations](OPERATIONS.md).

## Maize

### 1. Choose one maize profile

Run **one** block:

**B73/Ki11: 20-DAP endosperm, four pools/eight technical runs.**

```bash
export PROFILE=maize_B73_Ki11
```

**B73/W22: 14-DAP endosperm, six ear-pool libraries.**

```bash
export PROFILE=maize_B73_W22
```

Prepare the selected cohort:

```bash
export COHORT_DATA="$DATA_ROOT/$PROFILE"
mkdir -p "$COHORT_DATA"/{manifests,raw,reference,private,work,scratch,logs,receipts}
cp -n manifests/"${PROFILE}"_*.tsv "$COHORT_DATA/manifests/"
cp -n "config/${PROFILE}.example.yaml" "config/${PROFILE}.private.yaml"
export CONFIG="$PWD/config/${PROFILE}.private.yaml"
```

`cp -n` preserves any existing edited copy; inspect it before reuse. All maize
commands below apply to this `CONFIG`. Finish one profile, then repeat the
section separately for the other.

### 2. Complete the maize inputs and preflight

Edit the private YAML and operational manifests:

| Profile | Requirements before reference preparation |
|---|---|
| B73/Ki11 | B73 RefGen_v4/AGPv4 genome; exact `F1maize.FINAL.gff`; Panzea HapMap3.2.1 raw/nonimputed AGPv4 calls; verified source SHA256 values |
| B73/W22 | Compatible B73 RefGen_v4 genome/annotation; prepared opposite-homozygote markers; complete resource/marker-source contracts |

B73/W22 needs a prepared table with columns:

```text
chrom genomic_pos ref alt focal_index alternate_index
```

The allele indices must be opposite homozygous 0/1 states. Missing marker or
annotation assets must be supplied; retrieval does not discover them for you.
For B73/Ki11, review Panzea parents, duplicate lines and quality filters.
External inbred calls are not genotypes of the original experimental plants.

Review `strand_mode`/`strand_verified` from the real protocol. B73/Ki11's
configured `fr` setting needs independent review and pilot diagnostics.
Review ear/pool independence; technical runs do not increase biological n.
Use dosage **2:1 only for normal triploid endosperm**, not embryo/whole kernel.

Complete tool/LSF/storage review using [setup](SETUP.md), then:

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
```

**Continue when:** configuration/manifests pass and an accepted native inventory
is written. A pending inventory needs [installed-tool review](SETUP.md#installed-tool-identity).
Preflight does not certify reference availability or biological callability.

### 3. Retrieve the public maize inputs

On the approved internet-enabled host:

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

Wait for successful retrieval and inspect logs and `paths.resource_receipt`.
Required asset rows without verified byte sizes/checksums block acquisition.
Conditional assets must have their contracts completed before use.

**Continue when:** raw reads and required resources pass checks. Enter verified
**uncompressed** genome/annotation SHA256 in `reference.genome_sha256` and
`reference.annotation_sha256`. Complete operational marker-source SHA256 values
and annotation/source approvals. Finish these edits before reference work.

### 4. Prepare the reference, markers and index

```bash
bin/submit_pipeline.sh "$CONFIG" reference
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

The driver submits reference catalog, parental markers/transcript SNPs and one
shared Ornaments index in order. Modules 03 and 04 are already included.
Wait for all three jobs, then:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
```

**Continue when:** reference status is complete and the transcript FASTA,
catalog, markers and shared index validate.

### 5. Run and review the pilot

```bash
bin/submit_pipeline.sh "$CONFIG" pilot
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

Wait for every library pilot and the collector. Read `paths.pilot_report`.
Check strand agreement, pseudoalignment, processed fragments, elapsed time,
`peak_rss_kb`, index bytes and projected time for every library. Review total
raw/reference/index/work/scratch capacity. Endpoint work bytes are a lower
bound, not a continuously monitored disk maximum.

**Continue when:** pilot status is complete, all libraries pass and measured
resources fit the reviewed budget. Record the report SHA256 printed by the
collector. Changed settings require affected stages/pilot to be repeated.

### 6. Accept that report and quantify all reads

```bash
read -r -p "Accepted maize pilot-report SHA256: " PILOT_SHA
[[ "$PILOT_SHA" =~ ^[a-f0-9]{64}$ ]] || { printf 'Invalid SHA256\n' >&2; exit 1; }
bin/submit_pipeline.sh "$CONFIG" scale --accept "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --submit --accept "$PILOT_SHA"
```

The first command is still a dry-run, not authorization. Submission runs each
complete library, dependent parental assignment and a completion collector.
Wait for all jobs, then verify quantification completion:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
```

### 7. Fit models and inspect results

```bash
bin/submit_pipeline.sh "$CONFIG" models
bin/submit_pipeline.sh "$CONFIG" models --submit
```

After the model job completes:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
```

Read `paths.model_gene`, `paths.model_promoter` and `paths.model_pairs` under
`$COHORT_DATA`. Inspect eligibility/status before estimates or q-values.
Insufficient depth, unresolved independence, missing cross direction or invalid
rank/degrees of freedom yields a limitation, not a forced POE result.

## Mouse

### 1. Choose one mouse profile

Run **one** block:

**Babr/CAST: E7.5 extraembryonic ectoderm, five SE libraries.**

```bash
export PROFILE=mouse_Babr_CAST
```

**FVB/CAST: E12.5 placenta, four PE libraries.**

```bash
export PROFILE=mouse_FVB_CAST_E12
```

Prepare the selected cohort:

```bash
export COHORT_DATA="$DATA_ROOT/$PROFILE"
mkdir -p "$COHORT_DATA"/{manifests,raw,reference,private,work,scratch,logs,receipts}
cp -n manifests/"${PROFILE}"_*.tsv "$COHORT_DATA/manifests/"
cp -n "config/${PROFILE}.example.yaml" "config/${PROFILE}.private.yaml"
export CONFIG="$PWD/config/${PROFILE}.private.yaml"
```

Inspect existing copies before reuse. Run the other profile separately.
Do not add E16.5 libraries to E12.5 or mix tissues/crosses.

### 2. Complete mouse inputs and preflight

Both profiles need compatible **GRCm38** genome/annotation and joint MGP v5
calls with complete resource/marker-source contracts and SHA256 values.
Require called opposite homozygotes and parent-specific `FI=1`. Missing
positions in a strain-specific variant-only VCF are not reference calls.

| Profile | Specific review |
|---|---|
| Babr/CAST | Actual parent `C57BL_6J_Babr`; explicit `REFERENCE_C57BL_6J` proxy; CAST_EiJ calls; proxy approval/reason |
| FVB/CAST | Joint-source FVB_NJ and CAST_EiJ calls; PE protocol; animal/litter independence |

A compatible Babr genome-wide source is unresolved. The B6J-reference proxy is
an assumption, not a measured Babr genotype; do not substitute B6N.
For Babr/CAST retain the single assumed SE setting:

```yaml
fragments:
  fragment_mean_nt: 250
  fragment_sd_nt: 50
  fragment_parameter_origin: assumed_user_defined
```

This is not a measured distribution or tested sensitivity grid. The PE
FVB/CAST profile does not use SE arguments. Use dosage **1:1**. Review strand
and animal/litter independence: multiple placentas are not necessarily
independent litters. Supply missing reference/marker assets before use.

Complete [setup](SETUP.md), then:

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
```

**Continue when:** manifests and native inventory pass. Complete any pending
build review. Independence may block inference even if technical work proceeds.

### 3. Retrieve public mouse inputs

On the approved internet-enabled host:

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

Wait for successful retrieval and inspect `paths.resource_receipt`. Retrieval
does not find missing promoter annotations or establish Babr genotypes.

**Continue when:** reads/resources pass verification. Enter uncompressed
genome/annotation SHA256 values, marker-source SHA256 values and reviewed
annotation/provenance approvals before reference submission.

### 4. Prepare the reference and index

```bash
bin/submit_pipeline.sh "$CONFIG" reference
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

Reference, markers and shared index jobs are automatically linked in order.
Wait for all jobs, then:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
```

**Continue when:** reference completes and markers/index validate. Check actual
strain labels, contig/REF agreement and proxy qualifications.

### 5. Run and review the mouse pilot

```bash
bin/submit_pipeline.sh "$CONFIG" pilot
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

Wait for all libraries and the collector. Read `paths.pilot_report`; inspect
strand, pseudoalignment, fragments, SE/PE arguments, proxy flags, peak RSS and
projected time. Include raw/reference/index/scratch needs in the storage budget.
Endpoint work size is only a lower bound.

**Continue when:** pilot status is complete, all libraries pass and resource
review is complete. Do not claim robustness to untested fragment assumptions.

### 6. Accept the mouse report and quantify all reads

```bash
read -r -p "Accepted mouse pilot-report SHA256: " PILOT_SHA
[[ "$PILOT_SHA" =~ ^[a-f0-9]{64}$ ]] || { printf 'Invalid SHA256\n' >&2; exit 1; }
bin/submit_pipeline.sh "$CONFIG" scale --accept "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --submit --accept "$PILOT_SHA"
```

Wait for complete-library quantification, parental assignment and the collector:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
```

### 7. Fit models and inspect mouse results

```bash
bin/submit_pipeline.sh "$CONFIG" models
bin/submit_pipeline.sh "$CONFIG" models --submit
```

After the model job completes:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
```

Read `paths.model_gene`, `paths.model_promoter` and `paths.model_pairs` under
`$COHORT_DATA`. Keep Babr proxy/SE qualifications with its results. Unresolved
litter independence, insufficient coverage or invalid design blocks inference;
do not relax gates to obtain p-values.

## Human GUSTO

### 1. Prepare the private configuration and tables

```bash
export PROFILE=GUSTO
export COHORT_DATA="$DATA_ROOT/$PROFILE"
mkdir -p "$COHORT_DATA"/{manifests,raw,reference,private,work,scratch,logs,receipts}
cp -n manifests/GUSTO_*.tsv "$COHORT_DATA/manifests/"
cp -n config/GUSTO.example.yaml config/GUSTO.private.yaml
export CONFIG="$PWD/config/GUSTO.private.yaml"
```

Inspect existing copies before reuse. GUSTO uses local private inputs.
**Skip public retrieval.** Supplied GUSTO tables have headers only; populate
the private operational copies with independently checked real data.

### 2. Prepare verified linkage and local files

| YAML entry | Required input |
|---|---|
| `paths.sources` | Physical runs, real `local_1/local_2` FASTQs, layout and read counts |
| `paths.libraries` | Primary libraries, child IDs, strand/protocol and independent RNA/DNA identity checks |
| `paths.units` | Family/cluster independence, actual parents, tissue and stage |
| `paths.trios` | Exact child, mother, father and family linkage |
| `paths.covariates` | Exact sample rows and approved centered covariates named in `models.covariates` |
| `paths.typed_vcf` | Independently QC'd typed child/mother/father genotypes |
| `paths.personalization_vcf` | Approved child variants for personalization |

Arrange approved raw inputs under `$COHORT_DATA/raw` and private genotype/
covariate files under the cohort directory. This is a required local preparation
step: no private files have been moved, fetched or linked by this package.
Fill the exact real FASTQ paths rather than inventing a naming pattern.

Keep required schema columns, including common source/resource headers. Choose
one offspring per family independently of imbalance. Imputed child variants
can personalize the reference but do not become typed parental anchors.
Verify RNA/DNA identity and pedigree; filenames alone are not linkage evidence.
Use dosage **1:1**. Complete local genotype, covariate, privacy and storage review.

### 3. Verify the shared human references and tools

Keep these exact source references in the YAML:

```text
Genome (reference.genome_source):
/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon.fasta

Annotation (paths.annotation):
/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/reference/HPLRv2.0.annotated.PANTRY.gtf
```

These are read-only input exceptions to the project data root. Verify actual
availability and supplied SHA256 values. Do not guess a different GTF path.
Reference preparation copies the genome into `$COHORT_DATA/reference` before
indexing; shared source FASTA/index files are not rewritten.

Complete [setup](SETUP.md), keep aggregate release disabled, then:

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
```

**Continue when:** real linkage/manifests and tool identity pass and all private
inputs exist with approved access. Empty schemas or unresolved identity cannot
be bypassed. Example memory requests are templates, not GUSTO sizing results.

### 4. Build references and child-specific indices

```bash
bin/submit_pipeline.sh "$CONFIG" reference
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

The graph includes the shared catalog and typed/personalization preparation
plus an index for each configured library/child. Do not submit these again
manually. Wait for every job:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
```

**Continue when:** shared reference completes **and every child/index job
succeeds**. The shared-reference status alone does not certify all child
indices; pilot jobs also check their own child/index receipts.

### 5. Run and review all-library pilots

```bash
bin/submit_pipeline.sh "$CONFIG" pilot
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

Wait for all libraries and the collector. Inspect `paths.pilot_report`, strand
and typed-anchor callability, actual linked family count and measured
memory/time. Include all child indices, reads, references and scratch in the
budget. Raw participant count is not the callable promoter-pair denominator.

**Continue when:** all configured libraries pass, pilot status is complete and
measured resources fit the approved budget. Full processing is not automatic.

### 6. Accept the GUSTO report and quantify all reads

```bash
read -r -p "Accepted GUSTO pilot-report SHA256: " PILOT_SHA
[[ "$PILOT_SHA" =~ ^[a-f0-9]{64}$ ]] || { printf 'Invalid SHA256\n' >&2; exit 1; }
bin/submit_pipeline.sh "$CONFIG" scale --accept "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --submit --accept "$PILOT_SHA"
```

This is fresh raw-read quantification, not reuse of expression count tables.
Wait for every quantification/assignment job and completion collector:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
```

### 7. Fit and inspect human promoter models

```bash
bin/submit_pipeline.sh "$CONFIG" models
bin/submit_pipeline.sh "$CONFIG" models --submit
```

After the model job completes:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
```

Read `paths.model_promoter` and `paths.model_pairs` under `$COHORT_DATA`.
Human gene analysis is deferred: `paths.model_gene` records that status.

Inference needs depth at least 20/event, at least ten paired families, at least
three maternal-REF and three maternal-ALT orientations per sentinel, approved
centered covariates, full rank and at least five residual degrees of freedom.
HC3 t inference and 2,000 family bootstraps use seed 1147 with at least 90%
valid fits. Failed gates produce statuses, not forced estimates or p-values.

### 8. Optionally review and release local aggregates

Keep the completed analysis configuration fixed. For disclosure review, copy it:

```bash
cp -n "$CONFIG" config/GUSTO.release.private.yaml
export RELEASE_CONFIG="$PWD/config/GUSTO.release.private.yaml"
```

Inspect an existing release copy before reuse. Keep the same result paths.
First set the locally approved `exports.policy_approved` and integer
`exports.min_cell` in the release copy only; leave `exports.enabled` off.
No universal cell-size rule is supplied:

```bash
bin/Rscript_sif 07_exports/prepare_candidate.R "$RELEASE_CONFIG"
```

Inspect the candidate and exact printed SHA256. Complete manual and
complementary-suppression review. Only after approval, set the release copy's
`exports.approved_candidate_sha256`, `exports.review_token`,
`exports.complementary_suppression_reviewed` and `exports.enabled`:

```bash
bin/Rscript_sif 07_exports/release_aggregates.R "$RELEASE_CONFIG"
```

**Result:** the exact reviewed aggregate table at `paths.release`; no transfer
or publication occurs. Export actions check candidate schema/policy/hash, not
the full analysis receipts. Confirm completed analysis and result integrity
before review. Do not restart the analysis using the release copy: altered
settings have a different fingerprint.

All participant tables, local paths, logs, command plans, counts, indices,
receipts and candidate files remain on approved private storage.

## Reading the results

| Output key | Meaning |
|---|---|
| `paths.status` | File/receipt completion and blocker reasons; not live scheduler state |
| `paths.pilot_report` | All-library diagnostics and bounded resource projection |
| `paths.metrics` | Full quantification metrics after collection |
| `paths.annotation_audit` | Reference/annotation checks |
| `paths.model_gene` | Public gene results; human deferred status |
| `paths.model_promoter` | Promoter results and callability statuses |
| `paths.model_pairs` | Selected promoter-pair/sentinel inputs |
| `templates.assignment_qc` | Per-library mapping/orientation checks |
| `paths.deferred` | Endpoints outside implemented scope |

A successful technical job can still yield an ineligible biological comparison.
Missing evidence is NA, not zero. SNPs are not extra independent molecules and
technical runs are not extra biological replicates. Promoter structure and
expected allocations do not prove jointly parent/event-informative molecules,
imprinting or intragenomic conflict. Consider maternal seed-coat/decidua/blood
contamination, cis effects, proxy genotypes and biological independence.

For failures, use [operations](OPERATIONS.md). For local-test versus live-site
limitations, see [validation](VALIDATION.md).
