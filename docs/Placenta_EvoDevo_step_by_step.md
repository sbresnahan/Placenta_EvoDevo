# Placenta_EvoDevo: one step-by-step operating guide

## What to do with this guide

Complete **Part 1 once**, then use the complete recipe for the cohort you want to run. The five recipes do not depend on one another. You do not need to open another Markdown document.

This guide describes the existing package; it does not change its code, scientific methods, configuration examples, manifests, or ZIP. The commands are instructions for a future run on Seadragon. They have not been run against the cohorts as part of writing this guide.

**Do not paste the entire document into a terminal.** Run one numbered step at a time. A submission command only puts jobs in the queue: wait for completion and inspect the results before proceeding. Wherever the guide says **STOP**, the next stage is not justified until the stated evidence or file exists.

The examples contain real project paths. Text inside a non-executable YAML example that says `INSERT_...` is a description of information you must supply, not a usable value.

### A few terms used below

- **Terminal:** the command line on Seadragon.
- **Configuration:** a YAML file containing settings. Indentation matters; use spaces, not tabs.
- **Manifest:** a tab-separated text file listing inputs or samples. Keep actual tabs between columns, not spaces.
- **Library:** one primary RNA sample to quantify. Several sequencing runs or lanes can belong to the same library.
- **Biological unit:** the independently sampled organism, family, ear pool, or other justified replicate. Sequencing a unit twice does not make two biological replicates.
- **Reference:** the genome and transcript annotation used to interpret reads.
- **Marker:** a genomic SNP whose parental alleles are supported by the specified genotype evidence.
- **Receipt:** a small file recording successful completion and hashes of inputs or outputs. A file merely existing is not the same as a valid receipt.
- **SHA256 / MD5:** checksums identifying exact file bytes. A checksum proves byte identity, not biological suitability or the source history of software.
- **LSF:** Seadragon's job scheduler. `bsub` submits; `bjobs` inspects jobs.
- **Pilot:** the small technical run required before full quantification.

---

# Part 1. Shared setup — do this once

## 1. Open a Seadragon terminal and set the project location

Use your normal approved login method. Run:

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export DATA_ROOT="$PROJECT/data"
umask 077
ls -ld "$PROJECT" "$PROJECT/scripts"
```

The scripts must be in:

`/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/scripts`

If that directory already contains `poe.R`, `run_pipeline.R`, `bin/`, and `config/`, do not unpack over it.

**Only if the package has not been deployed:** place the previously supplied `Placenta_promoter_POE.zip` somewhere you can read on Seadragon. Its top-level folder is `scripts/`. The following commands ask for the actual ZIP location, show its checksum, and check the archive; they do not unpack it yet.

```bash
read -r -p "Full path to Placenta_promoter_POE.zip: " PACKAGE_ZIP
sha256sum "$PACKAGE_ZIP"
unzip -t "$PACKAGE_ZIP"
```

The checksum of the supplied ZIP is:

`e6aa306b4cdb47b7e0f6d846dec9abaa7f86100beb467a36e8c64ba686b45de0`

**STOP** if it differs or the archive check fails. Do not replace an existing installation. If the archive is correct and `scripts/` does not exist, unpack once:

```bash
if [[ ! -e "$PROJECT/scripts" ]]; then
  unzip "$PACKAGE_ZIP" -d "$PROJECT"
else
  printf 'Existing scripts directory left untouched.\n'
fi
```

Now enter the script directory:

```bash
cd "$PROJECT/scripts"
export SCRIPTS_DIR="$PWD"
export POE_SITE_FILE="$PWD/config/seadragon.sh"
ls -l poe.R run_pipeline.R bin/submit_pipeline.sh bin/Rscript_sif
```

Repeat these environment-setting commands whenever you open a new terminal. Later recipe commands assume this working directory. Setting `DATA_ROOT` does **not** rewrite paths inside a YAML file.

## 2. Make the working directories and editable copies

Run once. This makes separate data directories and copies the supplied tables and YAML examples only when an editable copy is absent.

```bash
mkdir -p "$DATA_ROOT/shared"
for COHORT in maize_B73_Ki11 maize_B73_W22 mouse_Babr_CAST mouse_FVB_CAST_E12 GUSTO; do
  mkdir -p "$DATA_ROOT/$COHORT"/{manifests,raw,reference,private,work,scratch,logs,receipts}
  for SOURCE in "$SCRIPTS_DIR/manifests/${COHORT}_"*.tsv; do
    DEST="$DATA_ROOT/$COHORT/manifests/$(basename "$SOURCE")"
    if [[ ! -e "$DEST" ]]; then
      cp -- "$SOURCE" "$DEST"
    fi
  done
  PRIVATE="$SCRIPTS_DIR/config/$COHORT.private.yaml"
  if [[ ! -e "$PRIVATE" ]]; then
    cp -- "$SCRIPTS_DIR/config/$COHORT.example.yaml" "$PRIVATE"
  fi
done
```

Inspect the result:

```bash
ls "$SCRIPTS_DIR/config/"*.private.yaml
ls "$DATA_ROOT/maize_B73_Ki11/manifests"
ls "$DATA_ROOT/GUSTO/manifests"
```

Edit **only the operational copies** in `data/<cohort>/manifests/` and the corresponding `.private.yaml` file. Do not edit the supplied manifest templates or `.example.yaml` files. If an editable copy already existed, inspect it rather than assuming it matches a fresh example.

Keep the private YAML files in `scripts/config/`. Their relative wrapper paths, such as `../bin/native/ornaments`, depend on that location.

For GUSTO, verify that institutional access controls restrict the data directory and private settings to authorized users. `umask 077` protects newly created files; it does not fix permissions on files that already existed. Do not make private data world-readable.

To edit a file, use `nano` followed by its path. In nano, **Ctrl+O**, then **Enter** saves; **Ctrl+X** exits. For example:

```bash
nano "$SCRIPTS_DIR/config/maize_B73_Ki11.private.yaml"
```

## 3. Check the existing software — do not reinstall it

Display the site settings:

```bash
less "$SCRIPTS_DIR/config/seadragon.sh"
```

Press `q` to leave `less`. The existing settings should point to:

- Module initialization: `/etc/profile.d/modules.sh`.
- Conda executable: `/risapps/rhel8/miniforge3/24.5.0-0/bin/conda`.
- samtools environment: `samtools-1.16.1`.
- bcftools environment: `bcftools-1.16`.
- Ornaments activation script: `/rsrch5/home/epi/bhattacharya_lab/software/ornaments/bin/activate`.
- Installed native commands: `ornaments` and `create_personalized_transcriptome`.
- R container: `/risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif`.
- R library: `/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1`.

Do not assume the Ornaments installation is a Conda environment: its configured activation is separate. Do not install or rebuild anything to bypass a failed check.

R runs inside Singularity through `bin/Rscript_sif`; native samtools, bcftools, Python/personalizer, Ornaments, and `bsub` run on the host. Both sides need access to `/rsrch5`, `/rsrch9`, and `/tmp`.

Check the R side:

```bash
bin/Rscript_sif -e 'p <- c("data.table","digest","yaml","curl","GenomicRanges","IRanges","S4Vectors","Rsamtools","Biostrings"); ok <- vapply(p, requireNamespace, logical(1), quietly=TRUE); print(data.frame(package=p, available=ok)); stopifnot(getRversion() >= "4.3.1", all(ok))'
```

**Success:** every package prints `TRUE`, R is at least 4.3.1, and the command exits without an error. **STOP** for a missing image, library, package, or native environment; ask the site maintainer to resolve the existing installation.

Collect the host software identities:

```bash
bin/host_inventory.sh "$DATA_ROOT/shared/software_inventory.tsv"
less -S "$DATA_ROOT/shared/software_inventory.tsv"
```

This inventories the actual executables, site file, R image, and Python dependency identity. It does not download cohort data or submit jobs. Hashing the container image can take time.

The output columns are `component`, `path`, `version`, and `sha256`. Locate the `ornaments` and `personalizer` rows. Their `sha256` values identify the installed executable and console entry point, not just this package's wrappers.

## 4. Obtain actual evidence for the installed Ornaments source

The package is pinned to source commit:

`2c4b1a2c18616f9de1cb307749688daaf3313d17`

The software directory is not automatically a source checkout. If `source_revision` in the inventory is unresolved, a binary hash alone is insufficient. Ask the person who installed Ornaments for build/install records linking the **actual binary and personalizer** to this exact commit.

There are two supported evidence routes:

1. A genuine clean source checkout at the pinned revision, identified through `POE_ORNAMENTS_SOURCE` in the site environment. This requires an actual checkout, not the environment directory. Site-environment changes require review; do not guess a source path.
2. A reviewed installation-provenance file tying the observed executable hashes to the pinned commit.

For the second route, open a new local record:

```bash
nano "$DATA_ROOT/shared/installation_provenance.yaml"
```

Use this structure, substituting the real inventory values and real evidence. This is a **template to edit**, not an approval to copy blindly:

```yaml
reviewed: true
commit: 2c4b1a2c18616f9de1cb307749688daaf3313d17
binary_sha256: INSERT_REVIEWED_ORNAMENTS_SHA256
personalizer_sha256: INSERT_REVIEWED_PERSONALIZER_SHA256
evidence: INSERT_ACTUAL_LOCAL_BUILD_OR_INSTALL_EVIDENCE
```

The evidence should identify retrievable local records explaining how those installed artifacts were linked to that commit. Merely writing “reviewed” is not evidence. **STOP** if the source-to-installation link cannot be established.

After the record is genuinely approved, copy it to each cohort's configured private location without overwriting an existing record:

```bash
for COHORT in maize_B73_Ki11 maize_B73_W22 mouse_Babr_CAST mouse_FVB_CAST_E12 GUSTO; do
  DEST="$DATA_ROOT/$COHORT/private/installation_provenance.yaml"
  if [[ ! -e "$DEST" ]]; then
    cp -- "$DATA_ROOT/shared/installation_provenance.yaml" "$DEST"
  fi
  sha256sum "$DEST"
done
```

Existing records must be reviewed against the current inventory too. A record from another installed build is not interchangeable.

In each private YAML, fill:

- `ornaments.binary_sha256`: the actual approved `ornaments` executable checksum.
- `ornaments.installation_provenance`: keep the cohort's `private/installation_provenance.yaml` path.
- `ornaments.installation_provenance_sha256`: the checksum printed for that exact record.

Keep `ornaments.commit` pinned. With a verified clean-checkout route, the runtime can establish source identity directly; do not fabricate an installation record when that route is used. In either case, the binary checksum must be reviewed and filled.

## 5. Establish the storage and scheduler facts before changing approval flags

Check available space and your scheduler access:

```bash
df -h "$DATA_ROOT"
bqueues
bjobs
```

The cohort-specific RNA byte totals below are **compressed RNA only**. They exclude genome/VCF downloads, personalized references, indexes, temporary files, and outputs. B73/Ki11 also has ten large chromosome VCFs. Allow space for all of these, not just FASTQs.

Ask the Seadragon administrator or use the approved local scheduler documentation to establish:

1. What units `bsub -M` uses at this site.
2. What units `rusage[mem=...]` uses.
3. Whether the requested memory is interpreted per job or per CPU slot.
4. Which compute nodes can access the project, existing tools, and Singularity.
5. Which queue may access the internet for public downloads.
6. Which queues, CPU allocations, memory amounts, and wall times are allowed for this project.
7. Whether work, scratch, logs, and receipts support normal filesystem operations, including atomic rename and exclusive directory locks. Do not point them at an object-storage mount.

**STOP** if these facts are unknown. The supplied values are starting requests, not measured requirements: most analysis stages request 8-class memory values and two CPUs, and retrieval requests 4-class values and one CPU. The raw numbers are unsafe until the site's units and accounting are confirmed.

## 6. Enter the common settings in each private YAML you intend to use

Open each intended file individually:

```bash
nano "$SCRIPTS_DIR/config/maize_B73_Ki11.private.yaml"
nano "$SCRIPTS_DIR/config/maize_B73_W22.private.yaml"
nano "$SCRIPTS_DIR/config/mouse_Babr_CAST.private.yaml"
nano "$SCRIPTS_DIR/config/mouse_FVB_CAST_E12.private.yaml"
nano "$SCRIPTS_DIR/config/GUSTO.private.yaml"
```

You can prepare only the cohort you intend to run now; do not approve unreviewed cohorts. Make the following concrete edits using the facts just established:

1. Verify that every cohort-owned path points to `data/<that cohort>/`. Keep the documented human source-reference exceptions in the GUSTO recipe.
2. Keep `execution.site_file: seadragon.sh` consistent with `POE_SITE_FILE`. Keep the configured image and R library consistent with the site settings.
3. Fill the approved Ornaments identity values from step 4.
4. Set `execution.posix_work_approved: true` only after the storage check.
5. Replace `lsf.memory_units`' unresolved description with the actual local units description. Set `lsf.memory_units_verified: true` and `lsf.memory_per_job_verified: true` only after confirming both limit/reservation units and per-job/per-slot behavior.
6. Set `lsf.execution_environment_verified: true` only after confirming compute-node access to tools, mounts, container, and environment propagation.
7. Under `lsf.stages`, inspect **each** of `retrieve`, `reference`, `pilot`, `scale`, and `models`. Enter locally permitted `memory_value`, `reservation_value`, `cpus`, `queue`, and `wall_time`. Memory and reservation need not be numerically identical if the site uses different units. Wall time is hours and minutes, for example `'03:00'`.
8. Keep `ornaments.threads` no greater than the CPU request for pilot and scale. The supplied value is two.
9. Review `pilot.max_projected_seconds` and `pilot.max_work_bytes` against local limits. These are technical gates, not measured resource needs.
10. Set `lsf.resource_limits_approved: true` only when the proposed requests, storage allowance, and project budget are locally authorized. This approves the limits, **not** the eventual pilot result. You must still review the pilot and explicitly accept its exact hash before scale.
11. For public cohorts only, set `retrieval.internet_host_approved: true` after confirming the configured retrieval queue runs on an approved internet-capable host.
12. Leave `exports.enabled: false` throughout analysis.

Do not change thresholds, parent identities, reference builds, or scientific choices merely to make a check pass. Cohort-specific annotation, marker, strand, independence, and human QC approvals are addressed in the recipes.

**Finish all analysis-setting and manifest edits before reference construction.** The package fingerprints the entire configuration, manifests, software identity, and relevant code. Even a resource or approval edit after the pilot can invalidate earlier receipts. If evidence forces a change, make it, then regenerate the affected prerequisites and pilot; never reuse an old pilot acceptance.

## 7. Understand submission and completion before running any cohort

The correct submission interface is:

`bin/submit_pipeline.sh "$CONFIG" STAGE --submit`

Without `--submit`, it is a **dry run**: it can create plans, registry files, and directories, but does not call `bsub`. R-only entry points can prepare a host plan without actually submitting; do not use them as a replacement for the shell submission driver.

Do not separately submit the individual `.lsf` modules as well as the driver's graph. That would duplicate work.

After each actual submission:

1. Read the terminal's `Submitted ... as ...` lines.
2. List the submission records and logs:

```bash
ls -lt "$DATA_ROOT/$COHORT/private/"jobs_*.submitted.tsv
ls -lt "$DATA_ROOT/$COHORT/logs/lsf/"
bjobs
bjobs -a
```

3. Open the **actual** record printed by the driver. The following prompts avoid pretending an example job ID or filename is real:

```bash
read -r -p "Full path to the actual .submitted.tsv record just printed: " SUBMITTED
less -S "$SUBMITTED"
read -r -p "Actual numeric LSF job ID to inspect: " JOB_ID
bjobs -l "$JOB_ID"
read -r -p "Full path to that job's actual .out or .err log: " JOB_LOG
less -S "$JOB_LOG"
```

The real job record is under `data/<cohort>/private/jobs_<unique>.tsv.submitted.tsv`. The planning registry is `jobs_<unique>.tsv`. Labels beginning `DRYRUN_` or `HOST_PLANNED_` are **not** scheduler IDs. Logs are under `data/<cohort>/logs/lsf/`; the filename contains the cohort, action, and actual job ID.

Inspect **every node** of the submitted graph, not only the first job. `PEND` means waiting; `RUN` means active; `DONE` indicates scheduler success; `EXIT` indicates failure. Disappearance from the default `bjobs` view is not proof of success. Check the retained scheduler history or site accounting if needed.

Then run:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$DATA_ROOT/$COHORT/status.tsv"
```

This checks on-disk completion receipts and hashes. It does **not** query LSF. Before the first stages run, later stages are expected to be `blocked`.

- After reference: expect `reference` to be `complete`.
- After pilot: expect `pilot` to be `complete`.
- After accepted scale: `scale` means `authorized`, not that quantification finished.
- After full quantification and assignment: expect `quantification` to be `complete`.
- After models: expect `models` to be `complete`.

For human reference jobs, shared `reference: complete` alone does not certify every child's personalized index; inspect all child/index jobs and their receipts.

### If anything fails

**Stop before the next stage.** Read the failed job's `.err` and `.out` logs, its `.submitted.tsv` record, and `status.tsv`.

- If submission stopped partway through, some earlier jobs may already exist. Reconcile the real submitted IDs with LSF **before** rerunning. Do not replay a `.host.sh` plan or delete its `.started` guard.
- Retrieval can resume valid partial downloads under its checked HTTP-range rules. It does not use an SRA fallback, trim reads, discover missing sources, or call new variants. A checksum/source problem requires resolving the real source, not bypassing validation.
- A stale receipt means inputs, settings, software, or outputs changed. Correct the cause and rerun the affected stages in order.
- Locks under `receipts/locks/` are not automatically discarded. Confirm the owning job and live processes are gone before any locally supervised removal.
- Never rebuild or replace a reference/index while downstream jobs are reading it.
- If a pilot requires a different strand or resource setting, change the configuration or manifest, regenerate stale prerequisites, and repeat the pilot. Do not accept the old report.

---

# Part 2. Maize B73/Ki11 — 20 DAP endosperm

## K1. Select the cohort

```bash
cd "$PROJECT/scripts"
export COHORT=maize_B73_Ki11
export CONFIG="$SCRIPTS_DIR/config/maize_B73_Ki11.private.yaml"
export COHORT_DIR="$DATA_ROOT/$COHORT"
```

This recipe uses four biological endosperm pools and eight technical sequencing runs:

- B73 maternal: `Sample_15`, runs ERR2889292/ERR2889293; `Sample_16`, ERR2889297/ERR2889298.
- Ki11 maternal: `Sample_13`, ERR2889302/ERR2889303; `Sample_14`, ERR2889307/ERR2889308.

The compressed RNA total in the supplied source table is **52,116,556,804 bytes**. The ten parental VCF downloads are additional large inputs. Do not count the eight runs as eight replicates.

## K2. Review the four pools and the protocol

```bash
less -S "$COHORT_DIR/manifests/maize_B73_Ki11_sources.tsv"
nano "$COHORT_DIR/manifests/maize_B73_Ki11_libraries.tsv"
nano "$COHORT_DIR/manifests/maize_B73_Ki11_units.tsv"
nano "$CONFIG"
```

Keep the supplied run-to-library grouping. The libraries are paired-end (`PE`). Their supplied `strand_mode` is `fr`, and the template marks it verified. Confirm that designation against the actual library protocol; a prefilled `TRUE` is not a replacement for your review. If unsupported, set `strand_verified` to `FALSE` and stop before pilot until resolved. The pilot also checks the selected orientation.

In the units table, confirm the four pool identities and their independent biological sampling. `cluster_id` identifies the biological sampling cluster, not a sequencing run. Retain `independence_verified: TRUE` only where supported by the study records.

Keep B73 as focal and Ki11 as alternate. Keep the maternal:paternal dosage at **2:1**, with `dosage.reviewed: true`, only for the normal triploid endosperm samples in this recipe.

## K3. Check setup, then retrieve the public inputs

The supplied resources table specifies B73 RefGen_v4/AGPv4, `F1maize.FINAL.gff` from Zenodo record 2611319, and ten raw AGPv4 chromosome VCFs. Inspect it:

```bash
less -S "$COHORT_DIR/manifests/maize_B73_Ki11_resources.tsv"
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
bin/submit_pipeline.sh "$CONFIG" retrieve
```

The first preflight validates configuration/manifests. The `--execute` preflight additionally checks host identities and the R environment; it does not download or submit. It should accept the reviewed installed software. **STOP** if either preflight fails. The dry run should describe one retrieval job with the approved queue and resource requests.

Then submit:

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

Wait for the one retrieval job to finish, using the scheduler and log checks in shared step 7. It retrieves the RNA and all resources marked `required=TRUE`, checks downloaded bytes and MD5, tests gzip inputs and FASTQ integrity, and creates uncompressed outputs where `decompress_to` is specified.

Inspect the recorded resources and actual files:

```bash
less -S "$COHORT_DIR/resource_receipt.tsv"
ls -lh "$COHORT_DIR/raw/"
ls -lh "$COHORT_DIR/reference/"
```

**STOP** for a failed download, checksum, gzip, or paired-read check. Do not proceed with partial inputs.

## K4. Enter the reference and marker checksums

The genome must now be `reference/genome.fa`. The source annotation must be `reference/source_annotation.gtf`. The latter filename does not change the source's GFF content; keep the exact configured annotation rather than substituting a convenient GTF.

Run:

```bash
sha256sum "$COHORT_DIR/reference/genome.fa" "$COHORT_DIR/reference/source_annotation.gtf"
sha256sum "$COHORT_DIR/reference/"hmp321_raw_agpv4_chr*.vcf.gz
nano "$COHORT_DIR/manifests/maize_B73_Ki11_marker_sources.tsv"
nano "$CONFIG"
```

In the marker-source table, replace each placeholder checksum with the matching chromosome VCF's actual SHA256. There are ten rows, for chromosomes 1 through 10. Do not put a checksum from one chromosome in another row.

In the YAML:

1. Put the **uncompressed** `genome.fa` checksum in `reference.genome_sha256`.
2. Put the exact source annotation checksum in `reference.annotation_sha256`.
3. Set `reference.annotation_approved: true` only after confirming this annotation is compatible with B73 RefGen_v4 and is the intended annotation.
4. Keep `markers.mode: panzea`, B73/Ki11 sample names, and the configured duplicate-line checks.
5. Set `markers.provenance_approved: true` only after reviewing the ten raw, nonimputed AGPv4 sources and accepting that they are external lines, not genotypes measured in the original experimental plants.

The existing marker preparation requires supported opposite homozygous parental SNPs, depth at least five, called-allele fraction at least 0.9, allowed site filters, and at least 99% REF agreement with the actual genome. It excludes NI5 and does not infer missing genotypes.

Finish the common scheduler/software settings now. Keep exports off. Rerun both preflights after these edits:

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
```

## K5. Build the reference and index

```bash
bin/submit_pipeline.sh "$CONFIG" reference
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

Expect **three jobs**: catalog/reference, parental markers, and personalized index, with dependencies. Wait for all three to succeed. Then:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/status.tsv"
less -S "$COHORT_DIR/annotation_audit.tsv"
less -S "$COHORT_DIR/events.tsv"
less -S "$COHORT_DIR/markers.tsv"
ls -lh "$COHORT_DIR/reference/personal.index"
```

**Proceed only if** `reference` is `complete`, all three jobs succeeded, the annotation audit is acceptable, and the index exists. The reference stage checks exon structure, transcript orientation, and coordinate mapping; it builds fixed-anchor 25/50/100 bp catalogs with 50 bp as primary.

## K6. Run and inspect the pilot

```bash
bin/submit_pipeline.sh "$CONFIG" pilot
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

Expect **five jobs**: four library pilots and one collector. Wait for all five. Then:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/pilot_report.tsv"
```

Require `pilot: complete`. The default pilot uses 100,000 fragments per strand diagnostic and a one-million-fragment throughput prefix per library. These are nonrandom prefixes, not representative random samples of the whole experiment.

Check `pseudoalignment_fraction`, strand comparisons, `elapsed_seconds`, `peak_rss_kb`, `index_bytes`, `projected_quant_seconds`, and `work_bytes_endpoint_lower_bound`. The configured pseudoalignment floor is 0.2 and strand tolerance is 0.05. The collector rejects configured compute/storage-limit violations. The work-byte measurement is an endpoint lower bound, **not** peak temporary disk usage; projected quantification time is not total pipeline wall time.

**STOP** if the protocol or resources are unresolved. If anything must change, make the change and rerun affected prerequisites and pilot before accepting.

## K7. Accept the exact pilot and quantify all four libraries

Only after reviewing the report and confirming the authorized resource limits:

```bash
PILOT_SHA=$(sha256sum "$COHORT_DIR/pilot_report.tsv" | cut -d ' ' -f 1)
printf 'Pilot selected for acceptance: %s\n' "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --accept "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --submit --accept "$PILOT_SHA"
```

The first scale command authorizes that exact report and prints a plan; it does not submit. Expect **nine submitted jobs** from the second command: four quantifications, four dependent parent assignments, and a collector.

Wait for all nine, then:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/status.tsv"
ls "$COHORT_DIR/work/Sample_15/complete/fr/"
less -S "$COHORT_DIR/work/Sample_15/complete/fr/assignment_qc.tsv"
```

The example inspection path assumes the reviewed `fr` orientation was retained. If it was changed, use that actual strand directory. Inspect assignment QC for all four libraries, not only Sample_15. Require `quantification: complete`; `scale: authorized` alone is insufficient.

## K8. Fit the models and inspect the results

```bash
bin/submit_pipeline.sh "$CONFIG" models
bin/submit_pipeline.sh "$CONFIG" models --submit
```

Wait for the **one model job**, then:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/model_gene.tsv"
less -S "$COHORT_DIR/model_promoter.tsv"
less -S "$COHORT_DIR/model_pairs.tsv"
```

Require `models: complete`. Retain ineligible/deferred statuses rather than treating them as zero effect. Interpret the public contrasts relative to the 2:1 maternal:paternal dosage, not a diploid 1:1 null. The estimated allocations and promoter-coordinate groups are not direct observations of individual jointly parent-and-event-assigned molecules, and they are not by themselves proof of imprinting or conflict.

---

# Part 3. Maize B73/W22 — 14 DAP endosperm

## W1. Select the cohort

```bash
cd "$PROJECT/scripts"
export COHORT=maize_B73_W22
export CONFIG="$SCRIPTS_DIR/config/maize_B73_W22.private.yaml"
export COHORT_DIR="$DATA_ROOT/$COHORT"
```

This recipe uses six paired-end ear pools from PRJNA623806:

- B73 maternal: SRR11509920, SRR11509919, SRR11509910.
- W22 maternal: SRR11509916, SRR11509915, SRR11509914.

The compressed RNA total is **18,265,330,625 bytes**. Keep the 14 DAP samples separate from the 20 DAP B73/Ki11 cohort.

## W2. Resolve the annotation and parental-marker blockers

**This profile is not reference-ready as supplied.** The B73v4 genome download is specified, but the compatible annotation and prepared parental-marker input are not.

Obtain an approved annotation compatible with `B73_RefGen_v4`, with transcript/exon relationships needed by the existing reference builder. Place the exact approved source at:

`data/maize_B73_W22/reference/source_annotation.gtf`

Do not borrow a differently built annotation or silently switch this cohort to another reference. Its biological suitability requires review, not merely a matching filename.

If the approved source is already available locally, copy it without overwriting an existing project file:

```bash
read -r -p "Full path to the approved B73 RefGen_v4 source annotation: " ANNOTATION_INPUT
cp -n -- "$ANNOTATION_INPUT" "$COHORT_DIR/reference/source_annotation.gtf"
```

Verify an existing destination against that approved source; `cp -n` leaves it untouched.

Obtain a real prepared marker table for B73 versus W22 on the same build. Its required columns are:

```text
chrom	genomic_pos	ref	alt	focal_index	alternate_index
```

These are one-based genomic positions, single-base A/C/G/T REF and ALT alleles, and opposite homozygous allele indices: `0` means REF and `1` means ALT. Focal is B73; alternate is W22. Do not create example SNP rows, assume absent variant records mean REF, or substitute imputed genotypes without an approved change of method.

Put that actual table under this cohort's `reference/` directory. Record its source and derivation locally. Then:

```bash
read -r -p "Full path to the real approved B73/W22 prepared marker table: " MARKER_FILE
sha256sum "$MARKER_FILE"
nano "$COHORT_DIR/manifests/maize_B73_W22_marker_sources.tsv"
nano "$CONFIG"
```

The marker-source file currently has only the header `source_id`, `path`, `build`, `sha256`. Add a row naming your real source, its actual absolute path, build **`B73_RefGen_v4`**, and the printed SHA256. Add more rows only if there are real additional sources. Keep `markers.mode: prepared`.

**STOP before reference** if either the approved annotation or the supported marker table is missing. Setting an approval to true does not supply those inputs. You may retrieve the specified public RNA/genome before resolving this blocker, but that does not make reference construction ready.

## W3. Review biological units, protocol, and setup

```bash
less -S "$COHORT_DIR/manifests/maize_B73_W22_sources.tsv"
nano "$COHORT_DIR/manifests/maize_B73_W22_libraries.tsv"
nano "$COHORT_DIR/manifests/maize_B73_W22_units.tsv"
less -S "$COHORT_DIR/manifests/maize_B73_W22_resources.tsv"
nano "$CONFIG"
```

Keep six biological pools and the supplied run-to-library relationships. Confirm independent ear sampling in `cluster_id` and `independence_verified`. The supplied strand is `rf` with `strand_verified=FALSE`. Verify the library protocol; set `strand_verified=TRUE` only with actual protocol evidence. Resolve any orientation disagreement before pilot.

Keep B73/W22 parental identities and the normal-endosperm maternal:paternal dosage **2:1**. Confirm the common software identity, storage, scheduler, and public internet approvals. Keep exports off.

## W4. Retrieve the specified RNA and genome

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
bin/submit_pipeline.sh "$CONFIG" retrieve
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

Run the submission line only after both preflights and the dry run succeed. Expect **one retrieval job**. Wait for completion and inspect its actual scheduler state and logs.

```bash
less -S "$COHORT_DIR/resource_receipt.tsv"
ls -lh "$COHORT_DIR/raw/"
ls -lh "$COHORT_DIR/reference/"
```

The supplied resources table downloads the B73v4 genome. It does not obtain an annotation or prepared marker table for you.

## W5. Finish the exact input approvals

When all three real inputs exist:

If this is a new terminal, first set `MARKER_FILE` again to the exact path recorded in the marker-source table.

```bash
sha256sum "$COHORT_DIR/reference/genome.fa" "$COHORT_DIR/reference/source_annotation.gtf"
sha256sum "$MARKER_FILE"
nano "$COHORT_DIR/manifests/maize_B73_W22_marker_sources.tsv"
nano "$CONFIG"
```

Enter the uncompressed genome and exact source annotation hashes in `reference.genome_sha256` and `reference.annotation_sha256`. Ensure the table's marker checksum matches the real file. Set `reference.annotation_approved: true` only after annotation/build review and `markers.provenance_approved: true` only after parental-source review. Retain the external-parent limitation where the genotype source is not the original experimental plants.

The marker preparation still validates the actual genome REF and requires at least 99% agreement. It does not repair a mismatched build.

Complete all settings now, then:

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
```

## W6. Build the reference

```bash
bin/submit_pipeline.sh "$CONFIG" reference
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

Wait for **three jobs**: reference catalog, markers, and index.

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/status.tsv"
less -S "$COHORT_DIR/annotation_audit.tsv"
less -S "$COHORT_DIR/markers.tsv"
ls -lh "$COHORT_DIR/reference/personal.index"
```

Require `reference: complete`, successful jobs, and an acceptable annotation audit. Stop on missing assets, genome REF disagreement, or an invalid source table.

## W7. Run the pilot and review all six libraries

```bash
bin/submit_pipeline.sh "$CONFIG" pilot
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

Wait for **seven jobs**: six pilots and the collector.

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/pilot_report.tsv"
```

Require `pilot: complete`. Review every library's strand diagnostics, pseudoalignment fraction, measured elapsed time and peak RSS, index size, projected quantification time, and work-byte lower bound. Defaults are 100,000 fragments for strand diagnostics and a one-million-fragment nonrandom throughput prefix. The selected strand must pass the configured 0.2 pseudoalignment floor and orientation checks, including the 0.05 tolerance where applicable.

Confirm available storage and the authorized budget; endpoint work size is not peak disk demand. If a protocol or resource change is necessary, edit and regenerate stale prerequisites/pilot. **Do not accept the previous pilot.**

## W8. Accept and run full quantification

```bash
PILOT_SHA=$(sha256sum "$COHORT_DIR/pilot_report.tsv" | cut -d ' ' -f 1)
printf 'Pilot selected for acceptance: %s\n' "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --accept "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --submit --accept "$PILOT_SHA"
```

Expect **13 jobs**: six quantifications, six parent assignments, and a collector. Wait for all to succeed.

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/status.tsv"
less -S "$COHORT_DIR/work/BW_endo1/complete/rf/assignment_qc.tsv"
```

Use the actual approved strand directory if not `rf`. Inspect assignment QC for all six libraries. Require `quantification: complete`, not only `scale: authorized`.

## W9. Fit and inspect the models

```bash
bin/submit_pipeline.sh "$CONFIG" models
bin/submit_pipeline.sh "$CONFIG" models --submit
```

Wait for **one model job**, then:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/model_gene.tsv"
less -S "$COHORT_DIR/model_promoter.tsv"
less -S "$COHORT_DIR/model_pairs.tsv"
```

Require `models: complete`. Interpret effects relative to the 2:1 dosage and preserve deferred/ineligible results. This is a separate cohort analysis; do not merge its raw contrast or multiple-testing denominator with B73/Ki11.

---

# Part 4. Mouse Babr/CAST — E7.5 extraembryonic ectoderm

## B1. Select the cohort

```bash
cd "$PROJECT/scripts"
export COHORT=mouse_Babr_CAST
export CONFIG="$SCRIPTS_DIR/config/mouse_Babr_CAST.private.yaml"
export COHORT_DIR="$DATA_ROOT/$COHORT"
```

This recipe includes five single-end wild-type extraembryonic ectoderm units:

- Babr maternal: SRR8360867 and SRR8360868.
- CAST maternal: SRR8360857, SRR8360873, and SRR8360872.

The compressed RNA total is **2,162,452,856 bytes**.

## B2. Resolve the genome, annotation, and MGP source

**The reference inputs are not complete as supplied.** Obtain a mutually compatible, approved GRCm38 genome and annotation. Place the actual uncompressed files at:

- `data/mouse_Babr_CAST/reference/genome.fa`
- `data/mouse_Babr_CAST/reference/source_annotation.gtf`

Keep `reference.build: GRCm38`. Do not substitute GRCm39 or an annotation on another assembly.

The marker resource is the joint MGP v5 GRCm38 VCF:

`https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz`

Its expected project location is:

`data/mouse_Babr_CAST/reference/mgp.v5.merged.snps_all.dbSNP142.vcf.gz`

Open the operational resource file:

```bash
nano "$COHORT_DIR/manifests/mouse_Babr_CAST_resources.tsv"
```

The supplied MGP row has `bytes=NA`, `md5=NA`, and `required=FALSE`. **It will not be downloaded in that state.** Choose one of these actual acquisition routes:

1. Obtain its genuine byte count and MD5 from a trusted source record, enter those values, keep the correct URL/destination, and set `required=TRUE` before retrieval; or
2. Supply an already approved, source-verified copy at the expected local path and leave download disabled for that row.

If using local inputs, use your actual approved paths rather than a guessed location. For example, after obtaining the files, the following prompts copy them into this cohort without overwriting existing files:

```bash
read -r -p "Full path to the approved uncompressed GRCm38 FASTA: " GENOME_INPUT
read -r -p "Full path to the approved GRCm38 source annotation: " ANNOTATION_INPUT
cp -n -- "$GENOME_INPUT" "$COHORT_DIR/reference/genome.fa"
cp -n -- "$ANNOTATION_INPUT" "$COHORT_DIR/reference/source_annotation.gtf"
```

For the local-VCF route only:

```bash
read -r -p "Full path to the approved joint MGP v5 GRCm38 VCF.gz: " MGP_INPUT
cp -n -- "$MGP_INPUT" "$COHORT_DIR/reference/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
```

If a destination already exists, verify it against the approved source; `cp -n` leaves it untouched. A self-computed hash must be tied to the reviewed source, not treated as proof of provenance. **STOP before reference** until all three inputs are real and reviewed.

## B3. Review the Babr proxy, independence, and single-end assumptions

```bash
less -S "$COHORT_DIR/manifests/mouse_Babr_CAST_sources.tsv"
nano "$COHORT_DIR/manifests/mouse_Babr_CAST_libraries.tsv"
nano "$COHORT_DIR/manifests/mouse_Babr_CAST_units.tsv"
nano "$CONFIG"
```

The actual focal parent is **`C57BL_6J_Babr`**. The configured marker proxy is **`REFERENCE_C57BL_6J`**, paired with `CAST_EiJ`. It is not a measured Babr genotype and must not be replaced by B6N.

In the YAML, retain `markers.proxy: true` only if you accept this explicit limitation. The example's `proxy_approved` is already true; retain it only after real local approval, and retain an accurate `proxy_reason`. Do not claim the proxy captures Babr-specific variants.

The units table starts with missing `cluster_id` and `independence_verified=FALSE`. Obtain litter/biological sampling information. Fill the actual clusters and set independence true only where justified. Do not assign unique cluster IDs merely to force five independent observations. If independence cannot be established, stop before inferential models; technical processing does not solve the design.

The supplied `unstranded` setting also starts unverified. Review the protocol and set `strand_verified=TRUE` only when supported. Keep layout `SE`.

This recipe uses **one** single-end fragment assumption: `fragment_mean_nt: 250`, `fragment_sd_nt: 50`, and `fragment_parameter_origin: assumed_user_defined`. These are not measurements. Do not describe them as measured or claim robustness across a sensitivity grid; no such grid is run. Keep dosage **1:1** for this diploid tissue.

Confirm all common software, scheduler, storage, resource-limit, and public internet approvals. Keep exports off.

## B4. Retrieve the public RNA and any enabled resource

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
bin/submit_pipeline.sh "$CONFIG" retrieve
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

Submit only after both preflights and the dry run succeed. Wait for **one retrieval job** and inspect its actual scheduler state and logs.

```bash
less -S "$COHORT_DIR/resource_receipt.tsv"
ls -lh "$COHORT_DIR/raw/"
ls -lh "$COHORT_DIR/reference/"
```

If the MGP row stayed `required=FALSE`, this job did not acquire it. A successful RNA download does not resolve missing genotype or reference assets.

## B5. Record the exact reference and marker identities

```bash
sha256sum "$COHORT_DIR/reference/genome.fa" "$COHORT_DIR/reference/source_annotation.gtf"
sha256sum "$COHORT_DIR/reference/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
nano "$COHORT_DIR/manifests/mouse_Babr_CAST_marker_sources.tsv"
nano "$CONFIG"
```

Replace the MGP source checksum placeholder with the actual approved VCF checksum. Its build must be `GRCm38`. Fill the YAML's uncompressed genome and exact annotation checksums.

Set `reference.annotation_approved: true` after the compatibility review, and `markers.provenance_approved: true` after reviewing the joint MGP source and B6J-reference proxy. Retain `markers.mode: mgp`.

The existing MGP preparation uses called parental `FI=1` evidence, not a blanket site-PASS rule. Missing records in a variant-only strain file are not inferred as REF; this is why the specified joint source matters.

Finish all settings now and rerun:

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
```

## B6. Build the reference

```bash
bin/submit_pipeline.sh "$CONFIG" reference
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

Wait for **three jobs**: catalog, markers, and index.

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/status.tsv"
less -S "$COHORT_DIR/annotation_audit.tsv"
less -S "$COHORT_DIR/markers.tsv"
ls -lh "$COHORT_DIR/reference/personal.index"
```

Require successful jobs and `reference: complete`. Stop for incompatible annotation, unsupported genotypes, or failed REF agreement.

## B7. Run the five-library pilot

```bash
bin/submit_pipeline.sh "$CONFIG" pilot
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

Wait for **six jobs**, including the collector.

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/pilot_report.tsv"
```

Require `pilot: complete`. Review protocol/strand checks, pseudoalignment fraction, elapsed time, peak RSS, index size, projected quantification time, and work-byte lower bound for every library. Confirm the report records the single-end 250/50 assumption and marker-proxy limitation.

The one-million-fragment throughput prefixes are nonrandom. Do not treat their performance or callability as a whole-cohort biological result. Endpoint work bytes do not measure peak disk usage. Resolve any protocol/resource failure and rerun stale prerequisites/pilot before accepting.

## B8. Accept the pilot and quantify all five units

```bash
PILOT_SHA=$(sha256sum "$COHORT_DIR/pilot_report.tsv" | cut -d ' ' -f 1)
printf 'Pilot selected for acceptance: %s\n' "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --accept "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --submit --accept "$PILOT_SHA"
```

Expect **11 jobs**: five quantifications, five assignments, and a collector. Wait for all.

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/status.tsv"
less -S "$COHORT_DIR/work/SRR8360867/complete/unstranded/assignment_qc.tsv"
```

Use the actual strand directory if the reviewed protocol changed it. Inspect all five assignment-QC files. Require `quantification: complete`.

## B9. Fit models only with justified biological independence

Recheck that the units table has genuine cluster/independence evidence, not convenience labels. If it is unresolved, stop here. If resolving it changes a manifest after pilot, regenerate the stale stages first.

```bash
bin/submit_pipeline.sh "$CONFIG" models
bin/submit_pipeline.sh "$CONFIG" models --submit
```

Wait for **one model job**, then:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/model_gene.tsv"
less -S "$COHORT_DIR/model_promoter.tsv"
less -S "$COHORT_DIR/model_pairs.tsv"
```

Require `models: complete`, and distinguish eligible results from design/depth failures. Retain the Babr proxy, assumed fragment distribution, small biological sample, and tissue-mixture limitations in any interpretation.

---

# Part 5. Mouse FVB/CAST — E12.5 placenta

## F1. Select the cohort

```bash
cd "$PROJECT/scripts"
export COHORT=mouse_FVB_CAST_E12
export CONFIG="$SCRIPTS_DIR/config/mouse_FVB_CAST_E12.private.yaml"
export COHORT_DIR="$DATA_ROOT/$COHORT"
```

This recipe uses four female paired-end placenta units from GSE75957 / PRJNA305855 / SRP067322:

- CAST maternal: SRR3085950 and SRR3085951, libraries `Pl_E12_5_CF_1` and `Pl_E12_5_CF_2`.
- FVB maternal: SRR3085952 and SRR3085953, libraries `Pl_E12_5_FC_1` and `Pl_E12_5_FC_2`.

The compressed RNA total is **22,013,620,572 bytes**. E16.5 runs SRR3085954–SRR3085957 are **not** extra E12.5 replicates; do not add them here.

## F2. Resolve the compatible reference and MGP resource

**The profile needs approved GRCm38 genome/annotation inputs and a completed MGP source.** Obtain mutually compatible GRCm38 files and put them at:

- `data/mouse_FVB_CAST_E12/reference/genome.fa`
- `data/mouse_FVB_CAST_E12/reference/source_annotation.gtf`

For actual local approved sources:

```bash
read -r -p "Full path to the approved uncompressed GRCm38 FASTA: " GENOME_INPUT
read -r -p "Full path to the approved GRCm38 source annotation: " ANNOTATION_INPUT
cp -n -- "$GENOME_INPUT" "$COHORT_DIR/reference/genome.fa"
cp -n -- "$ANNOTATION_INPUT" "$COHORT_DIR/reference/source_annotation.gtf"
nano "$COHORT_DIR/manifests/mouse_FVB_CAST_E12_resources.tsv"
```

Verify any existing destination against its approved source because `cp -n` does not overwrite it.

The resource table names:

`https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz`

The row currently has unknown bytes/MD5 and `required=FALSE`. It will be **skipped**, not automatically resolved. Either:

1. Obtain genuine source size and MD5, fill `bytes` and `md5`, and set `required=TRUE`; or
2. Supply the approved joint VCF locally at the configured destination and leave that download disabled.

For the second route only:

```bash
read -r -p "Full path to the approved joint MGP v5 GRCm38 VCF.gz: " MGP_INPUT
cp -n -- "$MGP_INPUT" "$COHORT_DIR/reference/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
```

Do not use separate variant-only files to infer absent parental calls as REF. Keep the joint v5 GRCm38 source. **STOP before reference** if any required real input or source evidence is missing.

## F3. Verify litter independence and the library protocol

```bash
less -S "$COHORT_DIR/manifests/mouse_FVB_CAST_E12_sources.tsv"
nano "$COHORT_DIR/manifests/mouse_FVB_CAST_E12_libraries.tsv"
nano "$COHORT_DIR/manifests/mouse_FVB_CAST_E12_units.tsv"
nano "$CONFIG"
```

Keep four E12.5 units, two per reciprocal direction. The supplied units have `cluster_id=NA` and `independence_verified=FALSE`. Obtain actual litter/biological sampling records. Fill cluster identities and independence approvals only when supported. Multiple placentas from one litter cannot be made independent by assigning new IDs.

The libraries are `PE`, supplied `rf`, initially unverified. Confirm orientation from protocol evidence, then set `strand_verified=TRUE`. Paired-end quantification does **not** use the YAML's single-end fragment mean/SD.

Keep the configured FVB_NJ and CAST_EiJ MGP sample identities and actual maternal/paternal labels. Keep diploid dosage **1:1**. This is not the Babr proxy analysis.

Confirm common software identity, local filesystem, scheduler, resource-limit, and internet approvals. Leave exports off.

## F4. Retrieve the public RNA and any enabled resources

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
bin/submit_pipeline.sh "$CONFIG" retrieve
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

Submit only after both preflights and dry run succeed. Wait for **one retrieval job**, inspect its actual logs and scheduler state, then:

```bash
less -S "$COHORT_DIR/resource_receipt.tsv"
ls -lh "$COHORT_DIR/raw/"
ls -lh "$COHORT_DIR/reference/"
```

If the MGP row remained `required=FALSE`, it was not downloaded. Confirm the locally supplied VCF exists before continuing.

## F5. Fill checksums and approve the actual sources

```bash
sha256sum "$COHORT_DIR/reference/genome.fa" "$COHORT_DIR/reference/source_annotation.gtf"
sha256sum "$COHORT_DIR/reference/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
nano "$COHORT_DIR/manifests/mouse_FVB_CAST_E12_marker_sources.tsv"
nano "$CONFIG"
```

Replace the marker-source SHA256 placeholder with the approved joint VCF checksum. Keep build `GRCm38`. Put the uncompressed genome and exact source annotation checksums into the YAML.

Set `reference.annotation_approved: true` only after compatibility review. Set `markers.provenance_approved: true` only after reviewing the joint MGP v5 source and FVB_NJ/CAST_EiJ genotype identities. Keep `markers.mode: mgp`. The existing preparation uses called parental `FI=1` evidence and checks genome REF agreement; it does not require every joint site to be blanket PASS.

Finish all settings and independence/protocol decisions now. Then:

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
```

## F6. Build the reference and index

```bash
bin/submit_pipeline.sh "$CONFIG" reference
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

Wait for **three jobs**: reference catalog, markers, and index.

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/status.tsv"
less -S "$COHORT_DIR/annotation_audit.tsv"
less -S "$COHORT_DIR/markers.tsv"
ls -lh "$COHORT_DIR/reference/personal.index"
```

Require successful jobs, acceptable annotation audit, and `reference: complete`.

## F7. Run and inspect the four-library pilot

```bash
bin/submit_pipeline.sh "$CONFIG" pilot
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

Wait for **five jobs**, including the collector.

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/pilot_report.tsv"
```

Require `pilot: complete`. Review all four libraries' strand checks, pseudoalignment, measured elapsed time/peak RSS, index size, projected quantification time, and work-byte lower bound. Defaults use 100,000-fragment strand diagnostics and a one-million-fragment nonrandom throughput prefix per library.

Confirm the approved resources can accommodate the actual measurements and projected work. Peak temporary storage is not measured by the endpoint work-byte lower bound. A strand or resource change invalidates the old acceptance; regenerate stale stages/pilot before proceeding.

## F8. Accept the exact report and run full quantification

```bash
PILOT_SHA=$(sha256sum "$COHORT_DIR/pilot_report.tsv" | cut -d ' ' -f 1)
printf 'Pilot selected for acceptance: %s\n' "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --accept "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --submit --accept "$PILOT_SHA"
```

Expect **nine jobs**: four quantifications, four assignments, and one collector. Wait for all:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/status.tsv"
less -S "$COHORT_DIR/work/Pl_E12_5_CF_1/complete/rf/assignment_qc.tsv"
```

Use the actual strand directory if different. Inspect assignment QC for every library. Require `quantification: complete`, not only scale authorization.

## F9. Fit and inspect models only after independence is justified

If litter independence remains unresolved, stop before modeling; retain that limitation explicitly. Any late manifest change requires regeneration of stale stages.

```bash
bin/submit_pipeline.sh "$CONFIG" models
bin/submit_pipeline.sh "$CONFIG" models --submit
```

Wait for **one model job**, then:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/model_gene.tsv"
less -S "$COHORT_DIR/model_promoter.tsv"
less -S "$COHORT_DIR/model_pairs.tsv"
```

Require `models: complete`. Preserve ineligible or deferred rows. These four female E12.5 placenta units do not support claims about E16.5, males, or placenta cell types separately.

---

# Part 6. Human GUSTO — private trio-linked placenta analysis

## H1. Select the cohort and keep it local

```bash
cd "$PROJECT/scripts"
export COHORT=GUSTO
export CONFIG="$SCRIPTS_DIR/config/GUSTO.private.yaml"
export COHORT_DIR="$DATA_ROOT/$COHORT"
```

**All supplied GUSTO tables contain headers only. There is no actual sample count or analysis denominator yet.** Do not infer participants from filenames or use another cohort's links. Do not run public retrieval for human data.

Work only in approved private storage. Do not transfer participant-level tables, identifiers, logs, counts, genotype files, personalized references/indexes, receipts, or review candidates outside the approved environment.

## H2. Assemble the real private inputs and verified identity crosswalk

Using your authorized local records, identify:

1. The exact RNA FASTQ files and lane/run grouping for each library.
2. The child, mother, father, and family identities.
3. Independent RNA/DNA identity QC and pedigree/genotype QC.
4. One independently selected primary RNA library per child and one selected child per family.
5. The typed trio VCF and the separate personalization VCF.
6. The prespecified covariates and their verified links to the RNA samples.

Place the actual RNA files under `data/GUSTO/raw/`, the typed VCF at `data/GUSTO/private/typed.vcf.gz`, and personalization VCF at `data/GUSTO/private/personalization.vcf.gz`, or explicitly update the corresponding private configuration paths to approved exact locations. Keep all cohort-owned inputs under the stated project data layout.

If the approved VCFs still need copying into this layout, use their real paths:

```bash
read -r -p "Full path to the approved typed trio VCF.gz: " TYPED_INPUT
read -r -p "Full path to the approved personalization VCF.gz: " PERSONAL_INPUT
cp -n -- "$TYPED_INPUT" "$COHORT_DIR/private/typed.vcf.gz"
cp -n -- "$PERSONAL_INPUT" "$COHORT_DIR/private/personalization.vcf.gz"
sha256sum "$COHORT_DIR/private/typed.vcf.gz" "$COHORT_DIR/private/personalization.vcf.gz"
```

Compare those exact project files with the approved source records. Existing destinations are not replaced by `cp -n`. For RNA, copy only files in the verified crosswalk, keeping mate identities intact; do not use a guessed wildcard over a private directory. If a file needs copying:

```bash
read -r -p "Full path to one verified RNA FASTQ.gz: " RNA_INPUT
read -r -p "Approved destination filename within GUSTO/raw, including .fastq.gz: " RNA_NAME
cp -n -- "$RNA_INPUT" "$COHORT_DIR/raw/$RNA_NAME"
ls -l "$COHORT_DIR/raw/$RNA_NAME"
md5sum "$COHORT_DIR/raw/$RNA_NAME"
```

Repeat for each actual mate/run file and record these exact destinations in `GUSTO_sources.tsv`. Confirm that each copy's size and checksum match its approved source before relying on it.

Do not overwrite files merely because their names resemble an example. The typed VCF supplies mother/father/child anchor evidence. The personalization VCF is a separate input; imputed personalization information must not replace typed parental anchor evidence.

**STOP** if the identity crosswalk, genotype QC, authorized access, or file provenance is missing. This package does not establish those links by guessing filenames.

## H3. Populate the three linked RNA/sample tables

Open the actual operational tables:

```bash
nano "$COHORT_DIR/manifests/GUSTO_sources.tsv"
nano "$COHORT_DIR/manifests/GUSTO_libraries.tsv"
nano "$COHORT_DIR/manifests/GUSTO_units.tsv"
```

### In `GUSTO_sources.tsv`

Keep these columns:

```text
run_id	sample_id	layout	url_1	url_2	md5_1	md5_2	bytes_1	bytes_2	read_count	local_1	local_2
```

Add one row per actual sequencing run/lane group:

- `run_id`: unique, stable run identifier.
- `sample_id`: exact library identifier used in the next table.
- `layout`: `PE` or `SE`, matching the real data.
- `local_1` and `local_2`: exact absolute paths to the actual mate files; mate 2 is absent for SE.
- `bytes_1`/`bytes_2` and `md5_1`/`md5_2`: actual local file sizes and checksums.
- `read_count`: actual complete-file read count in the source-table convention; PE totals count both mates, and the pilot projection divides that total by two to obtain fragments. SE totals are reads/fragments. Do not enter a prefix count or a guessed ENA-like total.
- Keep URL fields empty where there is no public source; human data are not fetched through public retrieval.

Obtain counts and hashes from verified local sequencing/QC records or an approved full-file check. Several run rows can point to one `sample_id`; they must not become extra independent libraries.

### In `GUSTO_libraries.tsv`

Keep:

```text
sample_id	unit_id	child_id	identity_verified	layout	strand_mode	strand_verified
```

Add one row per primary library:

- `sample_id` matches its source rows.
- `unit_id` links to the biological unit table.
- `child_id` matches the exact genotype/pedigree identifier.
- `identity_verified=TRUE` only after independent RNA/DNA identity QC.
- `layout` matches all runs in that library.
- `strand_mode` is the protocol-supported `fr`, `rf`, or `unstranded`.
- `strand_verified=TRUE` only after protocol review.

If any library is SE, the YAML must contain one finite positive fragment mean and SD and an honest `fragment_parameter_origin` (`measured` or `assumed_user_defined`). Do not label assumed 250/50 values as measured. PE libraries receive no SE fragment arguments.

### In `GUSTO_units.tsv`

Keep:

```text
unit_id	cluster_id	independence_verified	maternal_parent	paternal_parent	stage	tissue
```

Add one row per selected biological unit. Link `unit_id` exactly to the library table. Use the actual family/independence cluster, not the sequencing run. Record actual stage, tissue, and maternal/paternal identities consistently with the verified trio records. Set `independence_verified=TRUE` only after the one-child-per-family selection and biological independence review.

Identifiers used for run/sample/unit keys must use letters, numbers, underscores, periods, or hyphens. Keep one primary library per unit; do not use extra lanes as new units.

## H4. Populate pedigree and covariate tables

```bash
nano "$COHORT_DIR/manifests/GUSTO_trios.tsv"
nano "$COHORT_DIR/manifests/GUSTO_covariates.tsv"
nano "$CONFIG"
```

`GUSTO_trios.tsv` requires:

```text
child_id	mother_id	father_id	family_id
```

Add the verified trio for each selected child. The child, mother, and father IDs must be distinct and must identify the actual corresponding genotype samples. Keep one selected child per family and no duplicate child rows.

`GUSTO_covariates.tsv` starts with `sample_id`. Add the **approved covariate columns** and one correctly linked row per selected RNA sample. Center numeric covariates as required by the approved model; do not add post hoc covariates because they improve a result.

In the YAML:

- Put the exact approved column names in `models.covariates`. An empty list means a deliberately reviewed no-covariate model, not “unknown.”
- Set `human.covariates_approved: true` only after approving the model list, centered values, complete linked rows, and design.
- Set `human.genotype_qc_approved: true` only after actual genotype/pedigree QC.
- Keep dosage **1:1**, with the diploid-tissue review.
- Leave exports disabled.

Do not fill the public `GUSTO_marker_sources.tsv` with invented shared parental markers. Human anchors are prepared per child from typed trio evidence, not from a public fixed pair of strains. The resource table can remain header-only when no downloads are required.

## H5. Verify the exact existing human references

Keep these two read-only source references:

**Genome FASTA**

`/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon.fasta`

Expected SHA256:

`fb4243ebb014caf27111f24dd62b7ce42160f28581da6f8fcd6cba5977778d02`

**HPLRv2 annotation**

`/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/reference/HPLRv2.0.annotated.PANTRY.gtf`

Expected SHA256:

`88ac601ad75968427e4985b628f3eae4945e145e60ce1d37e1adb3a7b33884cc`

The configured GTF is on `/rsrch9`, not `/rsrch5`. Do not invent a different GTF path to make both references appear to reside on one filesystem.

Check the actual files:

```bash
sha256sum \
  /rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set_cleaned_ready_for_salmon.fasta \
  /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/reference/HPLRv2.0.annotated.PANTRY.gtf
nano "$CONFIG"
```

**STOP** if a file is unreadable or a checksum differs. Resolve the reference identity rather than replacing the expected hash without review.

Keep `reference.build: GRCh38`, `reference.genome_source` pointing to the shared FASTA, `paths.genome` pointing to this cohort's `reference/genome.fa`, and `paths.annotation` pointing to the exact shared HPLRv2 GTF. The reference job copies the shared genome into the project before indexing; it must not write an index beside the shared source.

If a project `reference/genome.fa` already exists, verify its checksum too; the reference job does not blindly replace an existing project genome. Keep the supplied expected genome/annotation hashes when the source files match. Set `reference.annotation_approved: true` only after confirming the intended HPLRv2/GRCh38 combination.

## H6. Complete preflight before scheduling

Finish the common software identity, scheduler, storage, resources, and all private table/model decisions before reference. Confirm that no human sample identifiers or paths have been placed in public documents.

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
```

Both must succeed. The initially empty GUSTO tables would fail; a successful preflight now means the required metadata and environment checks passed. It is not a replacement for independent RNA/DNA/genotype QC and does not certify every biological input.

**There is no human `retrieve` step.**

## H7. Build the shared catalog and every selected child's index

```bash
bin/submit_pipeline.sh "$CONFIG" reference
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

For **n selected libraries**, expect **1 + 2n jobs**: one shared reference catalog, one typed/personalization preparation per child/library, and one personalized index per child/library.

Wait for every job and inspect the actual `.submitted.tsv` record and LSF logs, not just the shared reference job.

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/status.tsv"
less -S "$COHORT_DIR/annotation_audit.tsv"
read -r -p "Actual sample_id from GUSTO_libraries.tsv to inspect: " SAMPLE_ID
ls -lh "$COHORT_DIR/private/$SAMPLE_ID/personal.fa" "$COHORT_DIR/private/$SAMPLE_ID/personal.index"
ls -lh "$COHORT_DIR/receipts/$SAMPLE_ID/child.rds" "$COHORT_DIR/receipts/$SAMPLE_ID/human.rds"
```

Repeat the child/index checks for **every selected library**. Shared `reference: complete` does not alone certify the child indexes; their jobs and receipts are required before pilot.

## H8. Run and inspect all human library pilots

```bash
bin/submit_pipeline.sh "$CONFIG" pilot
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

Expect **n + 1 jobs**, including the collector. Wait for all.

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/pilot_report.tsv"
```

Require `pilot: complete`. Review every selected library's strand evidence, pseudoalignment fraction, elapsed time/peak RSS, child-index size, projected quantification time, and work-byte lower bound. Keep all reports local.

The default diagnostics use 100,000 fragments per strand and a one-million-fragment nonrandom throughput prefix per library. A completed technical pilot is not a human association result or a sample-size guarantee. Endpoint work bytes underestimate peak storage; projected quantification time excludes other pipeline work.

If protocol, covariates, sample links, software, or resources need changing, make the approved correction and regenerate stale prerequisites and pilot. Do not use an old pilot hash after edits.

## H9. Accept the report and run full private quantification

Only after local resource review:

```bash
PILOT_SHA=$(sha256sum "$COHORT_DIR/pilot_report.tsv" | cut -d ' ' -f 1)
printf 'Pilot selected for acceptance: %s\n' "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --accept "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --submit --accept "$PILOT_SHA"
```

Expect **2n + 1 jobs**: quantification and dependent parent assignment for each library, plus a collector. Wait for all.

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/status.tsv"
read -r -p "Actual sample_id whose assignment QC you will inspect: " SAMPLE_ID
read -r -p "That library's approved strand (fr, rf, or unstranded): " STRAND_MODE
less -S "$COHORT_DIR/work/$SAMPLE_ID/complete/$STRAND_MODE/assignment_qc.tsv"
```

Inspect every library's assignment QC and require `quantification: complete`. Scale authorization alone is not completion. Keep private outputs and receipts in the approved environment.

## H10. Fit the human promoter models

```bash
bin/submit_pipeline.sh "$CONFIG" models
bin/submit_pipeline.sh "$CONFIG" models --submit
```

Wait for **one model job**, then:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
less -S "$COHORT_DIR/model_promoter.tsv"
less -S "$COHORT_DIR/model_pairs.tsv"
less -S "$COHORT_DIR/model_gene.tsv"
```

Require `models: complete`. A successfully completed model job can still report no eligible promoter comparisons. Do not turn ineligible or absent results into zero effects.

The existing human model requires, per comparison:

- At least depth 20 for each event.
- At least ten paired families.
- At least three maternal-REF and three maternal-ALT families at each sentinel.
- A full-rank model with at least five residual degrees of freedom after the approved covariates.

It fits a natural-log paired-promoter contrast with cis orientation and approved centered covariates. Only an exactly equal/opposite redundant cis term is removed. Inference uses HC3 t statistics and 2,000 family bootstraps, seed 1147, with at least 90% valid bootstraps required. The configured q ≤ 0.05 and absolute contrast ≥ log(2) flag is **conditional discordance**, not definitive imprinting.

`model_gene.tsv` is deliberately a `deferred_unimplemented` status for human gene modeling, not a table of human gene findings. Whole-gene/masking and TE comparisons are deferred; this run does not silently provide them.

## H11. Optional: prepare a separately reviewed aggregate release

**Do not do this as part of ordinary analysis.** Keep `GUSTO.private.yaml` fixed. First verify the completed analysis receipts, private result integrity, and eligibility. The release scripts check aggregate schema, policy, and approved hash; they do not independently validate the entire analysis receipt chain.

Obtain the institution's approved disclosure policy, actual integer minimum cell size, complementary-suppression review, and authorization to prepare aggregates. There is no universal default threshold supplied here.

Create a separate release configuration without replacing an existing one:

```bash
export RELEASE_CONFIG="$SCRIPTS_DIR/config/GUSTO.release.private.yaml"
if [[ ! -e "$RELEASE_CONFIG" ]]; then
  cp -- "$CONFIG" "$RELEASE_CONFIG"
fi
nano "$RELEASE_CONFIG"
```

Keep the same analysis/result paths. Edit only the release configuration's `exports` section:

1. Set `policy_approved: true` only after the local policy approval.
2. Set `min_cell` to the actual locally approved positive integer.
3. Keep `enabled: false` while preparing the review candidate.

Create the local candidate:

```bash
bin/Rscript_sif 07_exports/prepare_candidate.R "$RELEASE_CONFIG"
less -S "$COHORT_DIR/candidate.tsv"
sha256sum "$COHORT_DIR/candidate.tsv"
```

The candidate remains private. Review its exact contents and possible complementary disclosure against other releasable tables. Do not assume threshold filtering alone proves disclosure safety.

After the actual reviewer approves this exact candidate, edit the release configuration:

```bash
nano "$RELEASE_CONFIG"
```

Enter:

- `approved_candidate_sha256`: the exact approved candidate's printed hash.
- `review_token`: the real local approval/review identifier.
- `complementary_suppression_reviewed: true`: only after that review.
- `enabled: true`: only when local release-file creation is authorized.

Then:

```bash
bin/Rscript_sif 07_exports/release_aggregates.R "$RELEASE_CONFIG"
less -S "$COHORT_DIR/approved_exports/promoter_aggregates.tsv"
```

This writes the locally approved aggregate file. It does **not** authorize transfer, publication, or release of the private candidate or underlying participant-level material. Follow the institution's separate transfer/publication process.

---

# Part 7. What to conclude — and what not to conclude

After any recipe, the practical endpoint is a successful model job **plus** `models: complete`, with the input provenance, pilot review, assignment QC, eligibility statuses, and biological design limitations retained.

The package's existing analysis rules remain unchanged:

1. Missing parental evidence stays missing (`NA`), not zero. An uncallable event is not evidence of balanced expression.
2. Assignment uses same-event, same-genomic-locus expected allocations. SNP rows are not independent RNA molecules.
3. Sentinel and promoter-pair selection uses coverage rather than observed imbalance across all units. The configured depth threshold is 20.
4. The primary promoter grouping is fixed-anchor 50 bp, with 25/100 bp structural alternatives and first-donor compatibility; it is not a transitive chain of nearby starts.
5. Public models use `log2((M + 0.5)/(P + 0.5))` minus the log2 dosage ratio, with an unweighted intercept and reciprocal orientation term. Both directions, model rank, positive residual degrees of freedom, and justified independence matter.
6. A single cross direction cannot separate strain/cis effects from parent of origin. A technical replicate does not fix that problem.
7. Multiple-testing adjustment stays separate by cohort, stage, species, and endpoint. Do not combine raw counts or q-value denominators across these recipes.
8. Maternal seed-coat, decidual, blood, or other tissue admixture; cis effects; parental proxies; and cellular composition can affect interpretation. The package does not eliminate them by completing successfully.
9. Expected allocations and promoter-coordinate contrasts are candidate signals. They are not, alone, definitive imprinting or evidence that intragenomic conflict caused the signal.

If all requested stages are complete, keep the analysis configuration and operational manifests fixed with their receipts. For a later approved change, document the change and rerun the stages whose fingerprints are invalidated rather than modifying completed outputs by hand.
