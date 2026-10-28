# Placenta_EvoDevo walkthrough

Complete the shared setup, then follow the recipe for your cohort.

# Part 1. Shared setup

## 1. Set the project paths

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export DATA_ROOT="$PROJECT/data"
export SCRIPTS_DIR="$PROJECT/scripts"
export POE_SITE_FILE="$SCRIPTS_DIR/config/seadragon.sh"
umask 077
cd "$SCRIPTS_DIR"
```

## 2. Prepare a cohort

Define this terminal helper to select a cohort and create its working copies:

```bash
use_cohort() {
  export COHORT="${1:?Cohort name required}"
  export CONFIG="$SCRIPTS_DIR/config/$COHORT.private.yaml"
  export COHORT_DIR="$DATA_ROOT/$COHORT"
  mkdir -p "$COHORT_DIR"/{manifests,raw,reference,private,work,scratch,logs,receipts}
  for SOURCE in "$SCRIPTS_DIR/manifests/${COHORT}_"*.tsv; do
    DEST="$COHORT_DIR/manifests/$(basename "$SOURCE")"
    [[ -e "$DEST" ]] || cp -- "$SOURCE" "$DEST"
  done
  [[ -e "$CONFIG" ]] || cp -- "$SCRIPTS_DIR/config/$COHORT.example.yaml" "$CONFIG"
}
use_cohort maize_B73_Ki11
```

Edit sample tables in `$COHORT_DIR/manifests/` and settings in `$CONFIG`.
The helper preserves existing working copies.

## 3. Check the tools

```bash
00_setup/preflight.sh "$CONFIG"
```

Expected output: `Required tools and R packages are accessible.`

Setup checks the native tools, scheduler command, R image, library, and packages.
The existing `--execute` option performs the same check.

| Tool setting | Existing location |
|---|---|
| Ornaments environment | `/rsrch5/home/epi/bhattacharya_lab/software/ornaments` |
| samtools environment | `samtools-1.16.1` |
| bcftools environment | `bcftools-1.16` |
| Conda | `/risapps/rhel8/miniforge3/24.5.0-0/bin/conda` |
| R image | `/risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif` |
| R library | `/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1` |

Native tools run on the host. `bin/Rscript_sif` runs R and initializes its library.
Site settings are in `config/seadragon.sh`.

## 4. Job settings

The configuration stores queue, cores, memory, and wall time under `lsf.stages`.
Memory values are numeric GB, with matching `memory_value` and
`reservation_value`. [562] Public retrieval uses `transfer`. [561]

| Work | Queue | Typical cores / GB | Wall time |
|---|---|---|---|
| Retrieval | transfer | 1 / 4 | `02:00` |
| Reference, pilot, models | short | 2 / 8 | Set `02:59` |
| Quantification, assignment | medium | 2 / 8 | `12:00` |

For existing private configurations, set reference, pilot, and model
`wall_time` values to `'02:59'`.

Your site limits are: short under 3 hours, medium 4–24 hours, long 25–120 hours,
20 cores and 300 GB RAM per job. LSF enforces job limits. Pilot reports provide
usage and runtime projections for choosing resources.

## 5. Follow submissions

Each recipe lists the submission and the outputs to inspect when its jobs finish.

```bash
bjobs -a
ls -lt "$COHORT_DIR/logs/lsf/"
ls -lt "$COHORT_DIR/private/"jobs_*.submitted.tsv
bin/Rscript_sif poe.R "$CONFIG" status
```

The terminal prints each submitted job ID and the job-record filename.
View an individual job or log with:

```bash
read -r -p "Job ID: " JOB_ID
bjobs -l "$JOB_ID"
read -r -p "Log filename: " JOB_LOG
less -S "$JOB_LOG"
```

`status` reports completed stages and their required outputs. After an input or
setting changes, rerun the affected stages and their dependents. A command
without `--submit` prints a plan.

# Part 2. Maize B73/Ki11 — 20 DAP endosperm

## K1. Select the cohort

```bash
use_cohort maize_B73_Ki11
```

Four paired-end endosperm pools use eight sequencing runs:

| Maternal parent | Library | Runs |
|---|---|---|
| B73 | Sample_15 | ERR2889292, ERR2889293 |
| B73 | Sample_16 | ERR2889297, ERR2889298 |
| Ki11 | Sample_13 | ERR2889302, ERR2889303 |
| Ki11 | Sample_14 | ERR2889307, ERR2889308 |

The supplied strand is `fr`, dosage is maternal:paternal 2:1, and compressed RNA
totals 52,116,556,804 bytes. The tables already group runs into the four pools.

## K2. Retrieve the inputs

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

One job downloads RNA, the B73 RefGen_v4/AGPv4 genome, `F1maize.FINAL.gff`
from Zenodo record 2611319, and ten raw AGPv4 chromosome VCFs.
When it finishes:

```bash
less -S "$COHORT_DIR/resource_receipt.tsv"
ls -lh "$COHORT_DIR/raw/"
ls -lh "$COHORT_DIR/reference/"
```

Retrieval checks source MD5/bytes, gzip integrity, and FASTQ structure/mates.

## K3. Build the reference

```bash
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

Three jobs build the transcript/promoter catalog, B73/Ki11 markers, and index.
Marker preparation uses depth ≥5, allele fraction ≥0.9, and REF agreement ≥99%.
When the jobs finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/annotation_audit.tsv"
less -S "$COHORT_DIR/markers.tsv"
ls -lh "$COHORT_DIR/reference/personal.index"
```

Expected status: `reference complete`.

## K4. Run the pilot

```bash
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

Five jobs run four library pilots and collect their metrics. When they finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/pilot_report.tsv"
```

Expected status: `pilot complete`. Diagnostics use 100,000 fragments per strand
and a one-million-fragment throughput prefix per library.
The report includes pseudoalignment, time, RSS, index size, endpoint work bytes,
and projected quantification time. Projections are informative.

## K5. Quantify all four libraries

```bash
bin/submit_pipeline.sh "$CONFIG" scale --submit
```

Nine jobs perform quantification, parental assignment, and collection.
When they finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/work/Sample_15/complete/fr/assignment_qc.tsv"
```

Expected status: `quantification complete`. Each library has an
`assignment_qc.tsv` under its configured strand directory.

## K6. Fit the models

```bash
bin/submit_pipeline.sh "$CONFIG" models --submit
```

When the model job finishes:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/model_gene.tsv"
less -S "$COHORT_DIR/model_promoter.tsv"
less -S "$COHORT_DIR/model_pairs.tsv"
```

Expected status: `models complete`. Public effects use the 2:1 dosage offset.

# Part 3. Maize B73/W22 — 14 DAP endosperm

## W1. Select the cohort and retrieve RNA/genome

```bash
use_cohort maize_B73_W22
00_setup/preflight.sh "$CONFIG"
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

This cohort has six paired-end pools from PRJNA623806:
B73 maternal SRR11509920/SRR11509919/SRR11509910 and W22 maternal
SRR11509916/SRR11509915/SRR11509914. Compressed RNA totals 18,265,330,625 bytes.
The supplied strand is `rf` and dosage is 2:1.

One retrieval job downloads the configured RNA and B73v4 genome.
When it finishes:

```bash
less -S "$COHORT_DIR/resource_receipt.tsv"
ls -lh "$COHORT_DIR/reference/genome.fa"
```

## W2. Supply the annotation and markers

Provide a B73 RefGen_v4 annotation with transcript/exon relationships:

```bash
read -r -p "B73 RefGen_v4 annotation filename: " ANNOTATION_INPUT
cp -n -- "$ANNOTATION_INPUT" "$COHORT_DIR/reference/source_annotation.gtf"
read -r -p "B73/W22 prepared marker table filename: " MARKER_FILE
nano "$COHORT_DIR/manifests/maize_B73_W22_marker_sources.tsv"
```

The marker table uses one-based positions and columns:

```text
chrom	genomic_pos	ref	alt	focal_index	alternate_index
```

REF/ALT are single-base A/C/G/T alleles. Focal B73 and alternate W22 have
opposite homozygous indices, with 0=REF and 1=ALT.
In `marker_sources.tsv`, enter the source ID, absolute marker-table path,
and build `B73_RefGen_v4`. The configured adapter is `prepared`.

Sample settings are in:

```bash
nano "$COHORT_DIR/manifests/maize_B73_W22_libraries.tsv"
nano "$COHORT_DIR/manifests/maize_B73_W22_units.tsv"
```

`cluster_id` records the biological sampling cluster.

## W3. Build the reference

```bash
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

When the three reference/marker/index jobs finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/annotation_audit.tsv"
ls -lh "$COHORT_DIR/reference/personal.index"
```

Expected status: `reference complete`.

## W4. Run the pilot

```bash
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

When the seven pilot/collector jobs finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/pilot_report.tsv"
```

Expected status: `pilot complete`.

## W5. Quantify the six libraries

```bash
bin/submit_pipeline.sh "$CONFIG" scale --submit
```

When the 13 quantification/assignment/collector jobs finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/work/BW_endo1/complete/rf/assignment_qc.tsv"
```

Expected status: `quantification complete`.

## W6. Fit the models

```bash
bin/submit_pipeline.sh "$CONFIG" models --submit
```

When the model job finishes:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/model_gene.tsv"
less -S "$COHORT_DIR/model_promoter.tsv"
less -S "$COHORT_DIR/model_pairs.tsv"
```

Expected status: `models complete`.

# Part 4. Mouse Babr/CAST — E7.5 extraembryonic ectoderm

## B1. Select the cohort and retrieve RNA

```bash
use_cohort mouse_Babr_CAST
00_setup/preflight.sh "$CONFIG"
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

Five single-end wild-type units use Babr maternal SRR8360867/SRR8360868 and
CAST maternal SRR8360857/SRR8360873/SRR8360872.
Compressed RNA totals 2,162,452,856 bytes.
The supplied strand is `unstranded`, dosage is 1:1, and SE fragment
mean/SD are assumed 250/50 nt.

When the retrieval job finishes:

```bash
less -S "$COHORT_DIR/resource_receipt.tsv"
ls -lh "$COHORT_DIR/raw/"
```

## B2. Supply the GRCm38 reference and joint MGP source

```bash
read -r -p "GRCm38 FASTA filename: " GENOME_INPUT
read -r -p "GRCm38 annotation filename: " ANNOTATION_INPUT
read -r -p "Joint MGP v5 GRCm38 VCF.gz filename: " MGP_INPUT
cp -n -- "$GENOME_INPUT" "$COHORT_DIR/reference/genome.fa"
cp -n -- "$ANNOTATION_INPUT" "$COHORT_DIR/reference/source_annotation.gtf"
cp -n -- "$MGP_INPUT" "$COHORT_DIR/reference/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
nano "$COHORT_DIR/manifests/mouse_Babr_CAST_marker_sources.tsv"
nano "$COHORT_DIR/manifests/mouse_Babr_CAST_units.tsv"
```

The joint source is
`https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz`.
The marker-source row uses the local VCF path and build `GRCm38`.
Alternatively, enter the publisher's actual bytes/MD5 and `required=TRUE`
in the resources table to download this source during retrieval.

The actual focal parent is `C57BL_6J_Babr`; markers use the
`REFERENCE_C57BL_6J` proxy and `CAST_EiJ`. This proxy is recorded in the output.
MGP preparation uses called parental `FI=1` evidence.
Fill `cluster_id` from the biological/litter sampling records when available.

## B3. Build the reference

```bash
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

When the three catalog/marker/index jobs finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/annotation_audit.tsv"
ls -lh "$COHORT_DIR/reference/personal.index"
```

Expected status: `reference complete`.

## B4. Run the pilot

```bash
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

When the six pilot/collector jobs finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/pilot_report.tsv"
```

Expected status: `pilot complete`. Reports include the SE fragment parameters
and marker-proxy metadata.

## B5. Quantify the five libraries

```bash
bin/submit_pipeline.sh "$CONFIG" scale --submit
```

When the 11 quantification/assignment/collector jobs finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/work/SRR8360867/complete/unstranded/assignment_qc.tsv"
```

Expected status: `quantification complete`.

## B6. Fit the models

```bash
bin/submit_pipeline.sh "$CONFIG" models --submit
```

When the model job finishes:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/model_gene.tsv"
less -S "$COHORT_DIR/model_promoter.tsv"
less -S "$COHORT_DIR/model_pairs.tsv"
```

Expected status: `models complete`. Missing or shared cluster IDs produce
`independence_unresolved` descriptive results.

# Part 5. Mouse FVB/CAST — E12.5 placenta

## F1. Select the cohort and retrieve RNA

```bash
use_cohort mouse_FVB_CAST_E12
00_setup/preflight.sh "$CONFIG"
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

Four female paired-end placenta units come from GSE75957 / PRJNA305855 / SRP067322:
CAST maternal SRR3085950/SRR3085951 (`Pl_E12_5_CF_1/2`) and
FVB maternal SRR3085952/SRR3085953 (`Pl_E12_5_FC_1/2`).
Compressed RNA totals 22,013,620,572 bytes.
The supplied strand is `rf` and dosage is 1:1.

When the retrieval job finishes:

```bash
less -S "$COHORT_DIR/resource_receipt.tsv"
ls -lh "$COHORT_DIR/raw/"
```

## F2. Supply the GRCm38 reference and joint MGP source

```bash
read -r -p "GRCm38 FASTA filename: " GENOME_INPUT
read -r -p "GRCm38 annotation filename: " ANNOTATION_INPUT
read -r -p "Joint MGP v5 GRCm38 VCF.gz filename: " MGP_INPUT
cp -n -- "$GENOME_INPUT" "$COHORT_DIR/reference/genome.fa"
cp -n -- "$ANNOTATION_INPUT" "$COHORT_DIR/reference/source_annotation.gtf"
cp -n -- "$MGP_INPUT" "$COHORT_DIR/reference/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
nano "$COHORT_DIR/manifests/mouse_FVB_CAST_E12_marker_sources.tsv"
nano "$COHORT_DIR/manifests/mouse_FVB_CAST_E12_units.tsv"
```

The joint source is
`https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz`.
Use build `GRCm38` and the local VCF path in the marker-source table.
For public acquisition of this VCF, fill its actual publisher bytes/MD5
and set `required=TRUE` in the resources table.

The configured MGP samples are `FVB_NJ` and `CAST_EiJ`, using called `FI=1`
parental evidence. Enter actual litter/biological clusters in `cluster_id`
when available.

## F3. Build the reference

```bash
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

When the three reference/marker/index jobs finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/annotation_audit.tsv"
ls -lh "$COHORT_DIR/reference/personal.index"
```

Expected status: `reference complete`.

## F4. Run the pilot

```bash
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

When the five pilot/collector jobs finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/pilot_report.tsv"
```

Expected status: `pilot complete`.

## F5. Quantify all four libraries

```bash
bin/submit_pipeline.sh "$CONFIG" scale --submit
```

When the nine quantification/assignment/collector jobs finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/work/Pl_E12_5_CF_1/complete/rf/assignment_qc.tsv"
```

Expected status: `quantification complete`.

## F6. Fit the models

```bash
bin/submit_pipeline.sh "$CONFIG" models --submit
```

When the model job finishes:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/model_gene.tsv"
less -S "$COHORT_DIR/model_promoter.tsv"
less -S "$COHORT_DIR/model_pairs.tsv"
```

Expected status: `models complete`. Missing or shared cluster IDs produce
`independence_unresolved` descriptive results.

# Part 6. Human GUSTO — private trio-linked placenta analysis

## H1. Select the cohort and check tools

```bash
use_cohort GUSTO
00_setup/preflight.sh "$CONFIG"
```

The supplied human tables contain headers. Populate them with the actual
local RNA, genotype, sample, trio, and covariate records.

## H2. Populate the input tables

```bash
nano "$COHORT_DIR/manifests/GUSTO_sources.tsv"
nano "$COHORT_DIR/manifests/GUSTO_libraries.tsv"
nano "$COHORT_DIR/manifests/GUSTO_units.tsv"
nano "$COHORT_DIR/manifests/GUSTO_trios.tsv"
nano "$COHORT_DIR/manifests/GUSTO_covariates.tsv"
nano "$CONFIG"
```

| Table | Contents |
|---|---|
| Sources | One row per run/lane group; `run_id`, `sample_id`, layout, local FASTQ mate paths, bytes/MD5, and full-file read count |
| Libraries | `sample_id`, `unit_id`, `child_id`, `layout`, `strand_mode` |
| Units | `unit_id`, `cluster_id`, maternal/paternal IDs, stage, tissue |
| Trios | `child_id`, `mother_id`, `father_id`, `family_id` |
| Covariates | `sample_id` plus the columns named in `models.covariates` |

Use one primary library per child and one selected child per family.
PE `read_count` counts both mates; the projection converts it to fragments.
Source URLs can be empty for local human files.
The strand setting is `fr`, `rf`, or `unstranded`. SE libraries use the
configured fragment mean/SD and origin.

Set `paths.typed_vcf` to the typed trio VCF and `paths.personalization_vcf`
to the separate personalization VCF. The default destinations are:

```text
data/GUSTO/private/typed.vcf.gz
data/GUSTO/private/personalization.vcf.gz
```

Typed mother/father/child genotypes supply parental anchors. Center numeric
covariates for the configured model.

## H3. Use the existing human references

The configuration already points to these source files:

```text
/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon.fasta
/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/reference/HPLRv2.0.annotated.PANTRY.gtf
```

Build is `GRCh38`. Reference preparation copies the genome into the cohort
directory and writes derived references there.

## H4. Build the catalog and child indices

```bash
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

For n libraries, this submits 1 + 2n jobs: one catalog and preparation/index
jobs for each child. When all finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/annotation_audit.tsv"
read -r -p "sample_id to inspect: " SAMPLE_ID
ls -lh "$COHORT_DIR/private/$SAMPLE_ID/personal.index"
```

Expected status: `reference complete`, with completed child/index jobs.

## H5. Run the pilot

```bash
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

When the n + 1 pilot/collector jobs finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/pilot_report.tsv"
```

Expected status: `pilot complete`.

## H6. Quantify the libraries

```bash
bin/submit_pipeline.sh "$CONFIG" scale --submit
```

When the 2n + 1 quantification/assignment/collector jobs finish:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
read -r -p "sample_id to inspect: " SAMPLE_ID
read -r -p "Configured strand (fr, rf, unstranded): " STRAND_MODE
less -S "$COHORT_DIR/work/$SAMPLE_ID/complete/$STRAND_MODE/assignment_qc.tsv"
```

Expected status: `quantification complete`.

## H7. Fit the promoter models

```bash
bin/submit_pipeline.sh "$CONFIG" models --submit
```

When the model job finishes:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/model_promoter.tsv"
less -S "$COHORT_DIR/model_pairs.tsv"
```

Expected status: `models complete`. Human gene output records
`deferred_unimplemented`.

## H8. Optional aggregate export

Set `exports.min_cell` to the positive integer threshold for aggregate output.
Prepare the candidate:

```bash
nano "$CONFIG"
bin/Rscript_sif 07_exports/prepare_candidate.R "$CONFIG"
less -S "$COHORT_DIR/candidate.tsv"
```

Set `exports.enabled: true`, then write the aggregate file:

```bash
nano "$CONFIG"
bin/Rscript_sif 07_exports/release_aggregates.R "$CONFIG"
less -S "$COHORT_DIR/approved_exports/promoter_aggregates.tsv"
```

Export validates eligible model status, annotation identities, allowed columns,
and the cell/family threshold.

# Part 7. Results

| File | Contents |
|---|---|
| `annotation_audit.tsv` | Reference and transcript checks |
| `events.tsv`, `pairs.tsv` | Promoter groups and comparisons |
| `pilot_report.tsv` | Per-library diagnostics, measured usage, runtime projection |
| `assignment_qc.tsv` | Coordinate and parental-assignment status by library |
| `model_gene.tsv`, `model_promoter.tsv` | Effect estimates, uncertainty, q-values, eligibility status |
| `model_pairs.tsv` | Paired observations used by the model |

Primary promoter groups use fixed-anchor 50-bp windows; 25/100-bp catalogs are
also written. Public models use coverage-selected sentinels, dosage-adjusted
log2 ratios, and reciprocal-cross orientation. Missing parental evidence is NA.

Human promoter models use natural-log contrasts, cis orientation, and configured
covariates. Eligibility requires depth ≥20, ≥10 families, ≥3 maternal-REF and
maternal-ALT orientations at both sentinels, and ≥5 residual degrees of freedom.
Inference uses HC3 t statistics and 2,000 family bootstraps, seed 1147, with
≥90% valid fits. The q≤0.05 and absolute contrast≥log(2) flag denotes conditional
discordance.
