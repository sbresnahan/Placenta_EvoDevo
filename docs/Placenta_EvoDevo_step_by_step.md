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

Run once per environment (not per cohort):

```bash
tests/test_setup.sh
```

Expected last line: `Setup test passed: all required tools are available.`
It checks the host commands, native tools, `sratoolkit`, the R image, library
and packages. See [tests](../tests/README.md).

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
| Reference, pilot, models | short/medium | 2 / 8–24 | Match measured RSS and wall time |
| Quantification, assignment | medium | 2 / 8 | `12:00` |

For existing private configurations, choose the site-allowed queue and use
observed peak RSS to set memory. The Babr/CAST index required >8 GB on its
first attempt; the successful retry used a larger allocation.

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

## 6. Refit already quantified maize and mouse cohorts (November 2026 update)

The public-cross modeling change removes the **cluster-independence
eligibility gate** and relaxes **coverage callability**. It still requires
both reciprocal directions and positive residual degrees of freedom.
Sentinels are chosen by total depth and the number of callable units
(without inspecting parental imbalance); three callable units are sufficient.
Promoter event pairs need three *overlapping* callable units.

The model additionally rejects effectively zero residual variance. In
particular, tied allele ratios at three units can yield residual SEs around
`1e-16` and spuriously tiny p-values. Such fits now have status
`degenerate_residual_variance`: their estimates are descriptive, with
`NA` SE, confidence interval, p-value and BH q-value. Models with genuine
residual variation continue to be evaluated. These are still conditional
unit-level fits; the numerical safeguard does not establish biological
independence, validate genotype proxies, or resolve sex-specific effects.

The new public example YAML settings are:

```yaml
models:
  min_depth: 10.0
  min_callable_units: 3
```

Existing `*.private.yaml` files do **not** inherit these changes. Inspect
and edit the existing `models:` block in each private YAML, preserving
the rest of the file. Human GUSTO thresholds are unchanged.

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export DATA_ROOT="$PROJECT/data"
export SCRIPTS_DIR="$PROJECT/scripts"
cd "$SCRIPTS_DIR"

for COHORT in maize_B73_Ki11 maize_B73_W22 mouse_Babr_CAST; do
  CONFIG="$SCRIPTS_DIR/config/$COHORT.private.yaml"
  echo "=== $COHORT ==="
  grep -A7 '^models:' "$CONFIG"
  bin/Rscript_sif poe.R "$CONFIG" status
done
```

After editing these private files, archive the previous results:

```bash
for COHORT in maize_B73_Ki11 maize_B73_W22 mouse_Babr_CAST; do
  COHORT_DIR="$DATA_ROOT/$COHORT"
  ARCHIVE="$COHORT_DIR/archive/models_before_relaxed_callability_20261109"
  mkdir -p "$ARCHIVE"
  for FILE in model_gene.tsv model_promoter.tsv model_pairs.tsv; do
    if [ -f "$COHORT_DIR/$FILE" ]; then cp -n -- "$COHORT_DIR/$FILE" "$ARCHIVE/$FILE"; fi
  done
  if [ -f "$COHORT_DIR/receipts/model_receipt.rds" ]; then
    cp -n -- "$COHORT_DIR/receipts/model_receipt.rds" "$ARCHIVE/model_receipt.rds"
  fi
done
```

**Rerun models only**, assuming `quantification complete` for each cohort
and unchanged upstream inputs. There is no need to retrieve reads, rebuild
the reference or markers, rerun the pilot, or repeat quantification.

```bash
for COHORT in maize_B73_Ki11 maize_B73_W22 mouse_Babr_CAST; do
  CONFIG="$SCRIPTS_DIR/config/$COHORT.private.yaml"
  bin/submit_pipeline.sh "$CONFIG" models --submit
done
```

When the new jobs finish:

```bash
for COHORT in maize_B73_Ki11 maize_B73_W22 mouse_Babr_CAST; do
  CONFIG="$SCRIPTS_DIR/config/$COHORT.private.yaml"
  COHORT_DIR="$DATA_ROOT/$COHORT"
  echo "=== $COHORT ==="
  bin/Rscript_sif poe.R "$CONFIG" status
  for FILE in model_gene.tsv model_promoter.tsv; do
    awk -F '\t' 'NR>1 && $9 !~ /NA/ && $9 ~ /[0-9]/ {n++} END {print FILENAME, "p-values:", n+0}' "$COHORT_DIR/$FILE"
  done
done
```

If `mouse_FVB_CAST_E12` has also been quantified, perform the same
private-YAML edit, archive, and models-only rerun. The model treats each
`unit_id` as one observation; missing/repeated litter or ear IDs no longer
prevent fitting but may make uncertainty and p-values anti-conservative.
This is exploratory evidence, not independently validated imprinting.

### Refit after the residual-variance safeguard

For the three completed cohorts, install the patched code on Seadragon, inspect
whether previous model jobs have finished, and rerun **models only**. Before
submitting, save the old model tables and receipt in a new timestamped folder;
the worker will replace the tables and receipt even if `status` was already
`complete`. Obtain the update with `git pull --ff-only` after it has been
published to GitHub, or apply the reviewed commit patch locally using
`git am /path/to/patch` after updating to its base commit.

```bash
cd "$SCRIPTS_DIR"
git status --short  # Review tracked changes before updating code.
git pull --ff-only
# If this fix has NOT been published to GitHub, instead first apply its
# provided patch from a local file: git am /path/to/patch
STAMP="$(date +%Y%m%dT%H%M%S)"
for COHORT in maize_B73_Ki11 maize_B73_W22 mouse_Babr_CAST; do
  CONFIG="$SCRIPTS_DIR/config/$COHORT.private.yaml"
  COHORT_DIR="$DATA_ROOT/$COHORT"
  ARCHIVE="$COHORT_DIR/archive/models_before_variance_fix_$STAMP"
  mkdir -p "$ARCHIVE"
  for FILE in model_gene.tsv model_promoter.tsv model_pairs.tsv; do
    if [ -f "$COHORT_DIR/$FILE" ]; then cp -p -- "$COHORT_DIR/$FILE" "$ARCHIVE/$FILE"; fi
  done
  if [ -f "$COHORT_DIR/receipts/model_receipt.rds" ]; then
    cp -p -- "$COHORT_DIR/receipts/model_receipt.rds" "$ARCHIVE/model_receipt.rds"
  fi
  grep -A3 '^models:' "$CONFIG"
  bin/Rscript_sif poe.R "$CONFIG" status
done

# Submit after verifying no existing model jobs are still running or pending.
for COHORT in maize_B73_Ki11 maize_B73_W22 mouse_Babr_CAST; do
  bin/submit_pipeline.sh "$SCRIPTS_DIR/config/$COHORT.private.yaml" models --submit
done
```

The current `models.min_callable_units: 3` policy is retained for comparability
with exploratory results. For a more reliable **primary sensitivity analysis**,
archive the new three-unit tables, change `models.min_callable_units` to `4`
in each public cohort's private YAML, and submit `models --submit` again.
This changes sentinel/event-pair selection and may reduce the number of
testable features. Preserve both result sets and label the four-unit
results separately; changing a private YAML does not require upstream reruns.

After fitting, inspect the numerical guard and the number of testable rows:

```bash
for COHORT in maize_B73_Ki11 maize_B73_W22 mouse_Babr_CAST; do
  CONFIG="$SCRIPTS_DIR/config/$COHORT.private.yaml"
  COHORT_DIR="$DATA_ROOT/$COHORT"
  echo "=== $COHORT ==="
  bin/Rscript_sif poe.R "$CONFIG" status
  for TYPE in gene promoter; do
    awk -F '\t' 'NR>1 {
      status=$10; gsub(/"/,"",status); counts[status]++
      if($9!="NA" && $9!="") n++
    } END {
      printf "%s: testable p-values=%d; degenerate=%d\n", FILENAME,n+0,counts["degenerate_residual_variance"]+0
    }' "$COHORT_DIR/model_$TYPE.tsv"
  done
done
```

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

The submission creates eight SRA-run jobs, thirteen reference-resource jobs
(genome, `F1maize.FINAL.gff` from Zenodo record 2611319, ten AGPv3
chromosome VCFs, and an AGPv3-to-AGPv4 liftover chain), plus one final receipt job.
At most two retrieval jobs run
concurrently by default (`retrieval.max_parallel_jobs`). When all finish:

```bash
less -S "$COHORT_DIR/resource_receipt.tsv"
ls -lh "$COHORT_DIR/raw/"
ls -lh "$COHORT_DIR/reference/"
```

HTTPS resources are checked against expected byte sizes during download.
CyVerse VCFs instead undergo gzip-integrity checks and can be staged in
`reference/cyverse_downloads/` under their original filenames. SRA conversion
writes size-only receipts and skips expensive post-conversion FASTQ scans.

## K3. Build the reference

```bash
tests/test_inputs.sh "$CONFIG"
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
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

This cohort has six paired-end pools from PRJNA623806:
B73 maternal SRR11509920/SRR11509919/SRR11509910 and W22 maternal
SRR11509916/SRR11509915/SRR11509914. Compressed RNA totals 18,265,330,625 bytes.
The supplied strand is `rf` and dosage is 2:1.

Independent transfer jobs retrieve the configured RNA and B73v4 genome;
a dependent final job writes the receipt. When they finish:

```bash
less -S "$COHORT_DIR/resource_receipt.tsv"
ls -lh "$COHORT_DIR/reference/genome.fa"
```

## W2. Supply the annotation and reuse the HapMap3 genotypes

Provide a B73 RefGen_v4 annotation with transcript/exon relationships:

```bash
read -r -p "B73 RefGen_v4 annotation filename: " ANNOTATION_INPUT
cp -n -- "$ANNOTATION_INPUT" "$COHORT_DIR/reference/source_annotation.gtf"
```

The B73/Ki11 acquisition already downloaded ten AGPv3 HapMap3 imputed VCFs
containing B73 and W22 and the official AGPv3-to-B73v4 liftover chain. **Do not
download these genotype files again.** The B73/W22 marker-source manifest
references the original VCFs read-only, and writes its own processed B73/W22
marker tables under `maize_B73_W22/reference/`.

Verify those shared inputs before submitting the B73/W22 reference stage:

```bash
KI11_REF="$DATA_ROOT/maize_B73_Ki11/reference"
for CHR in {1..10}; do
  test -s "$KI11_REF/hmp321_imputed_agpv3_chr${CHR}.vcf.gz" || \
    echo "MISSING: HapMap3 chromosome $CHR"
done
test -s "$KI11_REF/AGPv3_to_B73_RefGen_v4.chain.gz" || \
  echo "MISSING: AGPv3-to-B73v4 chain"
zgrep -m 1 '^#CHROM' "$KI11_REF/hmp321_imputed_agpv3_chr1.vcf.gz" \
  | tr '\t' '\n' | grep -E '^(B73|W22)$'
```

The configured `panzea` adapter filters opposite homozygous B73/W22 SNPs,
enforces GQ where available, lifts them to B73 RefGen_v4, and checks the
resulting REF alleles against the B73 genome. The new B73/W22 marker tables
use one-based positions and `chrom`, `genomic_pos`, `ref`, `alt`,
`focal_index`, `alternate_index` columns (0=REF, 1=ALT). B73/Ki11's existing
marker tables are **not** reused because they compare different parents.

If you already created `maize_B73_W22.private.yaml` or a working copy of the
marker manifest before this update, refresh their **marker settings only**:
the `use_cohort` helper intentionally does not overwrite private files.

Sample settings are in:

```bash
nano "$COHORT_DIR/manifests/maize_B73_W22_libraries.tsv"
nano "$COHORT_DIR/manifests/maize_B73_W22_units.tsv"
```

`cluster_id` records the biological sampling cluster.

## W3. Build the reference

```bash
tests/test_inputs.sh "$CONFIG"
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
Alternatively, enter the publisher's actual byte size and `required=TRUE`
in the resources table to download this source during retrieval.

The actual focal parent is `C57BL_6J_Babr`; markers use the
`REFERENCE_C57BL_6J` proxy and `CAST_EiJ`. This proxy is recorded in the output.
MGP preparation uses called parental `FI=1` evidence.
Fill `cluster_id` from the biological/litter sampling records when available.

## B3. Build the reference

```bash
tests/test_inputs.sh "$CONFIG"
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

Expected status: `models complete`. The public model uses all callable
units, even if litter identities are unavailable. Interpret its standard
errors conditional on possible biological clustering.

# Part 5. Mouse FVB/CAST — E12.5 placenta

## F1. Select the cohort and retrieve RNA

```bash
use_cohort mouse_FVB_CAST_E12
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
For public acquisition of this VCF, fill its actual publisher byte size
and set `required=TRUE` in the resources table.

The configured MGP samples are `FVB_NJ` and `CAST_EiJ`, using called `FI=1`
parental evidence. Enter actual litter/biological clusters in `cluster_id`
when available.

## F3. Build the reference

```bash
tests/test_inputs.sh "$CONFIG"
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

Expected status: `models complete`. Biological clustering remains a
provenance limitation, not a fitting gate.

# Part 6. Placental fish — *Poeciliopsis prolifica* Line A × Line B

**Current readiness: staged, not yet executable.** The public RNA study
is [GSE318355 / SRP674561](https://www.ncbi.nlm.nih.gov/Traces/study?acc=SRP674561).
The study [Zhang et al., *Science Advances*, 2026](https://doi.org/10.1126/sciadv.aeh5487)
examined reciprocal crosses between semi-inbred Line A and Line B.
The target fetal placental tissue is the **pericardial sac**, distinct from
the embryonic tissue and maternal follicle libraries in the same study.
The authors reported 84 strongly supported imprinted protein-coding genes.
These are published results, **not** findings generated by this pipeline.

[Dryad 10.5061/dryad.ns1rn8q8r](https://doi.org/10.5061/dryad.ns1rn8q8r)
provides the authors' STAR/TopHat alignment, GATK SNP calling, and allele
count scripts plus `DataS1.csv` (placenta), `DataS2.csv` (embryo),
`DataS3.csv` (maternal follicle) and other analysis tables. It does **not**
list a downloadable ready-to-use joint parental genotype VCF or aligned BAM.
The reference is the [Line A genome/RefSeq annotation GCF_027474105.1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10655195/).

The separate [Line A PacBio HiFi project PRJNA1509661](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1509661)
contains a genomic Revio run **SRR40078298** (SRX34719858,
SAMN62272763; ~12.3 Gb / ~3.2 GB SRA download). It can independently
validate Line A variation and support methylation analysis after alignment
and variant calling. It is **Line A only**; it does not provide Line B DNA
genotypes. The identifier is recorded in
`manifests/fish_Pprolifica_AB_pacbio.tsv` and is not included among
placental RNA libraries.

## P1. Create working configuration and retrieve the run catalog

```bash
use_cohort fish_Pprolifica_AB
mkdir -p "$COHORT_DIR/private"
```

On a Seadragon **transfer** node, fetch ENA run-level metadata for this RNA study.
Submit this small metadata transfer from a login node:

```bash
bsub -q transfer -n 1 -M 4 -R 'rusage[mem=4]' -W 01:00 \
  -J fish_ena_catalog \
  -o "$COHORT_DIR/logs/fish_ena.%J.out" \
  -e "$COHORT_DIR/logs/fish_ena.%J.err" \
  /usr/bin/curl -fSL --retry 5 --connect-timeout 30 --max-time 600 \
  -o "$COHORT_DIR/private/fish_Pprolifica_AB_ENA.tsv" \
  'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP674561&result=read_run&fields=run_accession,study_accession,sample_accession,experiment_accession,library_name,library_layout,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes,sample_alias&format=tsv'
```

Check the LSF job is `DONE` and the TSV has nonzero records. The catalog
is **not** automatically classified into fetal placenta/embryo/maternal tissue
or reciprocal direction; the original GEO sample-level metadata must be
cross-checked against the paper's sample table.

## P2. Explicitly curate biological units and reciprocal directions

```bash
python3 bin/prepare_fish_manifests.py \
  --runinfo "$COHORT_DIR/private/fish_Pprolifica_AB_ENA.tsv" \
  --scaffold "$COHORT_DIR/manifests/fish_Pprolifica_AB_curation.tsv"
nano "$COHORT_DIR/manifests/fish_Pprolifica_AB_curation.tsv"
```

Set `include=YES` **only** for placental **pericardial sac** RNA runs,
with `tissue=pericardial_sac`. Assign the biological `unit_id`, actual
`maternal_parent` and `paternal_parent` (`Line_A` or `Line_B`),
the single analyzed developmental `stage`, validated `strand_mode`,
and the `cluster_id` if available. Technical runs from one biological
sample must share a unit ID. Everything else remains `include=NO`.
Unknown maternal parent is a blocker to parental assignment.

After verifying the sample crosswalk:

```bash
python3 bin/prepare_fish_manifests.py \
  --runinfo "$COHORT_DIR/private/fish_Pprolifica_AB_ENA.tsv" \
  --build "$COHORT_DIR/manifests/fish_Pprolifica_AB_curation.tsv" \
  --output-dir "$COHORT_DIR/manifests"
bin/submit_pipeline.sh "$CONFIG" reference  # validation and dry-run plan
```

The helper refuses to replace previously populated manifests unless
`--force` is supplied. It refuses nonplacental tissues, missing directions,
or fewer than three curated units. The included `fish_Pprolifica_AB_sources.tsv`,
`libraries.tsv` and `units.tsv` in the repository start **header-only on
purpose**; no unverified SRR cross assignment is invented.

## P3. Stage the reference and obtain a verified parental marker table

Acquire the matching RefSeq genomic FASTA and GTF from assembly
`GCF_027474105.1`. When NCBI Datasets CLI is available on a transfer node:

```bash
datasets download genome accession GCF_027474105.1 \
  --include genome,gtf --filename "$COHORT_DIR/reference/ncbi_dataset.zip"
unzip -l "$COHORT_DIR/reference/ncbi_dataset.zip" | grep -E 'genomic.fna|genomic.gtf'
```

Unpack and stage the genomic FASTA and GTF as
`$COHORT_DIR/reference/genome.fa` and
`$COHORT_DIR/reference/source_annotation.gtf`; check that FASTA and GTF
contig identifiers agree. If the CLI is unavailable, download the same
NCBI Datasets package elsewhere and upload it to the cohort reference folder.

The `prepared` marker adapter expects
`$COHORT_DIR/reference/lineA_lineB.fixed_exonic.markers.tsv`, with columns:

```text
chrom  genomic_pos  ref  alt  focal_index  alternate_index
```

Each site must be a **verified fixed difference** between Line A and Line B,
with `0=REF`, `1=ALT`, opposite homozygous indices, one-based positions,
and the same RefSeq assembly as the FASTA. The article called these from
Line A/Line B RNA-seq using GATK 3.6; the Dryad SNP scripts and DataS1–S3
can guide reconstruction, but RNA-based genotype calls must be clearly
distinguished from independently genotyped parental DNA.

For stronger Line A validation, retrieve **SRR40078298** on a transfer node,
align HiFi reads with an appropriate PacBio-compatible mapper, and call
SNVs using a validated HiFi caller on a compute queue. No processed VCF/BAM
was identified in the cited BioProject or Dryad file lists. Reconcile Line A
HiFi calls with the published fixed-SNP table and obtain Line B calls from
the original RNA genotyping workflow or separate Line B DNA. Do not treat
Line A HiFi alone as a two-parent genotype panel. Genotype-calling tools
and suitable LSF resources must be verified before creating those jobs.

## P4. Run the standard stages after the readiness gates

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
tests/test_inputs.sh "$CONFIG"   # after the retrieval jobs finish
bin/submit_pipeline.sh "$CONFIG" reference --submit
bin/submit_pipeline.sh "$CONFIG" pilot --submit
bin/submit_pipeline.sh "$CONFIG" scale --submit
bin/submit_pipeline.sh "$CONFIG" models --submit
```

Wait for each stage's jobs and validate its outputs before starting the next.
The fish's two-parent prepared markers, reference, RNA sample identities,
and tissue labels are prerequisites; the code does not infer them.
For the fetal pericardial sac, the starting allelic dosage ratio is 1:1
unless the paper or ploidy data demonstrate a different expectation.

# Part 7. Human GUSTO — precomputed transcript VCF and parental-origin matrix

## H1. Select the cohort

```bash
use_cohort GUSTO
```

Use the new precomputed GUSTO adapter in `config/GUSTO.example.yaml`.
It uses the original ESPRESSO transcript reference, already-phased child
transcriptomic VCF (including biallelic SNPs and indels for personalization),
chromosome-specific Left/Right maternal/paternal assignments,
and the private `J` FASTQ-to-`B` infant correspondence. The old trio genotype
extraction and chromosome-origin inference are bypassed. See
[`REVISED_RUNBOOK.md`](REVISED_RUNBOOK.md#6-human-gusto--precomputed-phased-transcript-variants-and-chromosome-parental-origin-private)
for full configuration paths and QC.

## H2. Populate the input tables

```bash
python3 bin/prepare_gusto_manifests.py --help
# Run the helper with the full VCF, matrix, ID map, FASTQ directory and
# verified strand mode as shown in docs/REVISED_RUNBOOK.md.
# Existing tables are protected unless --overwrite is explicitly supplied.
```

| Table | Contents |
|---|---|
| Sources | One row per FASTQ ID; `run_id=J1010`, `sample_id=J1010`, local mate paths and file sizes; total read count may be blank |
| Libraries | `sample_id=J1010`, `child_id=B010-04002`, layout and verified strand mode |
| Units | One biological unit per child; stage, tissue and source labels |
| Trios | `child_id`, `family_id`; parent genotypes are not reprocessed |
| Covariates | `sample_id` plus the columns named in `models.covariates` |

Verify one independent family per `family_id` before models. A generated
family ID initially equals the infant ID and must be corrected if siblings
or other non-independent biological observations are present. PE fragment
counts are validated during quantification. The pilot runtime projection is
`NA` if full-file `read_count` is unavailable, which does not block the GUSTO
pilot. Keep the sample linkage and counts private. Center numeric configured
covariates as needed.

## H3. Use the existing human references

The updated GUSTO example points to these inputs:

```text
/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon.fasta
/rsrch5/home/epi/bhattacharya_lab/data/Placenta_scRNAseq/salmon/ESPRESSO_corrected_SC_filtered.fasta
/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered/ESPRESSO_corrected_SC_filtered.gtf
```

The precomputed VCF ID (column 3) encodes hg38 coordinates. VCF `CHROM` and
`POS` are transcript coordinates. The VCF header uses a delimiter, so
`010-04002_B010-04002` denotes infant `B010-04002`. The reference stage copies the curated
FASTA into the cohort directory and verifies transcript IDs and spliced
lengths against the ESPRESSO GTF. The per-child map also checks genomic
positions against the exon coordinates during assignment.

## H4. Build the catalog and child indices

```bash
tests/test_inputs.sh "$CONFIG"
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

For n libraries, this submits 1 + 2n jobs: one catalog and per-child VCF
selection/personalization/index jobs. It uses the **full** private VCF path,
not the locally uploaded excerpt. When all finish:

If the shared reference completed but any child/index job failed, follow the
[GUSTO recovery instructions](REVISED_RUNBOOK.md#6-human-gusto--precomputed-phased-transcript-variants-and-chromosome-parental-origin-private):
audit malformed boundary deletion placeholders and use
`reference --resume --sample J1001 --submit` for one test child, then
`reference --resume --submit` for the remaining incomplete jobs. This skips
the completed shared reference and retains sequence-resolved indels.

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

# Part 8. Results

| File | Contents |
|---|---|
| `annotation_audit.tsv` | Reference and transcript checks |
| `events.tsv`, `pairs.tsv` | Promoter groups and comparisons |
| `pilot_report.tsv` | Per-library diagnostics, measured usage, runtime projection |
| `assignment_qc.tsv` | Coordinate and parental-assignment status by library |
| `model_gene.tsv`, `model_promoter.tsv` | Effect estimates, uncertainty, q-values, eligibility status |
| `model_pairs.tsv` | Paired observations used by the model |

Primary promoter groups use fixed-anchor 50-bp windows; 25/100-bp catalogs are
also written. Public models use coverage-selected sentinels (≥10 assigned
counts in ≥3 units), dosage-adjusted log2 ratios, and reciprocal-cross
orientation. Both orientations remain required. Unknown litter/ear clustering
no longer blocks exploratory fits. Missing parental evidence is NA.

Human promoter models use natural-log contrasts, cis orientation, and configured
covariates. Eligibility requires depth ≥20, ≥10 families, ≥3 maternal-REF and
maternal-ALT orientations at both sentinels, and ≥5 residual degrees of freedom.
Inference uses HC3 t statistics and 2,000 family bootstraps, seed 1147, with
≥90% valid fits. The q≤0.05 and absolute contrast≥log(2) flag denotes conditional
discordance.
