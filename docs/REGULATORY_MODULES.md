# Opt-in TE–TSS POE and human regulatory-context modules

These modules address immediately testable hypotheses **#1** and **#5**.
No primary or `additional_v1` tables, gates or receipts are rewritten. All new
receipts and analysis tables default to `additional_v2` beside `model_gene`.
Network transfers and LSF submission run on the **host**, not in the R image.

## Configure once, then run

Synchronize missing regulatory keys from the species-aware cohort examples
into existing schema-2 private configs with `python3 bin/sync_regulatory_configs.py --create-missing` (preview first with `--dry-run`).
See `docs/REGULATORY_CONFIG_SYNC.md`. The standalone
`config/regulatory.example.yaml` remains an illustrative fragment, not a full config.
Set paths and the complete human input bundle before starting. Regulatory
receipts track successful completion and required output presence; they do not
rehash inputs, scripts or configuration between stages. The workflow assumes
inputs are frozen for an analysis. Use a separate versioned `output_dir` only
for an intentionally distinct analysis, not for routine code updates.
Paths in manifests are relative to the cohort config directory, not the
manifest directory. Explicit absolute paths are recommended.

```bash
# From the actual checkout, with CONFIG already set to the existing cohort YAML:
: "${CONFIG:?Set CONFIG to your actual configured cohort YAML}"
export SCRIPTS_DIR="$(pwd)"
# Inspect (prepares scripts only); append --submit to submit.
bin/submit_pipeline.sh "$CONFIG" regulatory-resources
bin/submit_pipeline.sh "$CONFIG" regulatory-resources --submit
# After successful resource completion and frozen primary models:
bin/submit_pipeline.sh "$CONFIG" te-tss --submit
# Human only; after te-tss completes:
bin/submit_pipeline.sh "$CONFIG" human-context --submit
bin/Rscript_sif poe_regulatory.R "$CONFIG" regulatory-status
```

`poe.R` and `run_pipeline.R` both route these opt-in actions. The standalone
`poe_regulatory.R` also accepts them. Only `bin/submit_pipeline.sh ... --submit`
executes the returned host submission plan; R never calls `bsub`.
`--after NUMERIC_JOB_ID` adds a strict `done()` dependency. Jobs have separate
stdout/stderr files and scheduler IDs under `additional_v2/jobs`.
A host plan can be executed only once to prevent accidental duplicate jobs.
Submission failure does not cancel jobs already submitted: inspect their IDs.
An incomplete module or missing/empty required output fails without marking the
primary pipeline failed. Existing partial outputs are never overwritten implicitly.
Completed receipts remain usable after a code change.

### Minimum bundles and acquisition

**#1, all six cohorts:** completed frozen model, quantification and assignment
receipts; exact transcript/event maps and reference FAI; maternal/paternal
expected event allocations at shared loci; unit/cluster labels; human family,
cis-orientation and covariate linkage; full-reference repeat insertions.
The modules do not rerun RNA retrieval, assignment or quantification.

Species defaults download UCSC hg38/mm10 `rmsk.txt.gz` and its schema, or the
native B73v4 TE GFF3. Human resources additionally download the three **human
hTSC** GSE200761 CUT&Tag peak sets (H3K27ac/H3K4me1/H3K4me3), not mTSC peaks.
Fish defaults download only the RefSeq GCF_027474105.1 FASTA/assembly report.
The existing cluster image is reused; no image or repeat-classification
database is downloaded or installed automatically.
No whole-atlas, raw RNA, long-read or chromosome-wide LD download is implied.

The `regulatory_*example.tsv` files contain headers only, with no invented
local paths or identifiers. Populate them with actual reviewed input rows.
The resources manifest enables local-first replacement by `resource_id`;
omitted defaults remain enabled. Add URL and
optional MD5 URL to download a replacement. Supported transfer formats are
rmsk, gff3, bed, sql, text, fasta_gz and h5_gz. Plain or gzip BED/GFF are
accepted. Missing declared local files are errors, not fallback downloads.
Transfers use bounded retries, partial files, gzip/HTML checks and upstream
MD5 when available. Final manifests retain URLs, bytes, SHA256 and retrieval
time. Existing bad files are not overwritten. Local overrides require
`scope=full_reference` or `subset`; unhit subset sites are **not_assessed**.
Contig aliases are conservative (`chr`/case normalization or an explicit
one-to-one alias table); ambiguous aliases and off-reference coordinates fail.

**#5, human:** #1 outputs plus authorized local atlas lead/PIP/credible-set
or supplied colocalization tables, an explicitly reviewed phenotype crosswalk,
and optional CAGE/ATAC/chromatin BED inputs. Use the `regulatory_*example.tsv`
manifests. Export atlas outputs as TSV/gzip or supply Parquet with R `arrow`.
No stable public URL or permission to private atlas data is assumed.
Requested missing files block execution with the missing path. Public peaks
do not substitute for QTL. To explicitly request chromatin-only output, set
`require_qtl: false` and omit QTL/crosswalk keys.

QTL manifest kinds: `lead`, `pip`, `credible_set`, `coloc`. Variant rows require
`coordinate_system=one_based`, explicit `chrom`, `pos` and variant ID columns.
The adapter recognizes atlas `pip_all`, `pip_cs`, `cs_index`, `sig_in` and
`lead_snp`/`lead_pos`; retain all original fields as `source_*`. Optional
`modality_column`, `phenotype_column`, `ancestry_column`, `variant_column`,
`chrom_column`, `pos_column`, `pip_column`, `pip_cs_column`, `cs_index_column`,
`beta_column`, `se_column`, `ref_column`, `alt_column`, `cs_start_column`,
`cs_end_column` map other schemas. All tables require modality/phenotype IDs
and ancestry (or a manifest label `joint`/`unspecified`). `a0/a1` are preserved
as allele0/allele1, **not** assumed genomic REF/ALT. No variant-ID parsing,
liftover or allelic effect flipping is guessed.

Crosswalk columns are modality, phenotype_id, gene_id, event_id, mapping_level
(`event`/`gene`) and mapping_status (`verified`/`ambiguous`/`unmapped`). Exact
gene/event agreement is enforced. Nonunique phenotype mappings are withheld,
not resolved arbitrarily. Gene-level links are explicitly gene context, not
validation of a TSS phenotype. Peak manifests require a TRUE/FALSE
`selection_used` flag to expose annotation reuse, including CAGE selection.

### GUSTO: PANTRY placenta cis-QTL inputs

`bin/prepare_gusto_qtl_context.py` turns the PANTRY tensorQTL results into the
QTL tables, manifest and phenotype crosswalk above. Per ancestry (EAS, EUR) and
modality (`alt_TSS`, `isoforms`, `expression`) it keeps every significant
permutation lead (`qval < 0.05` in `*_top.tsv`) and every conditionally
independent signal (`*_independent_top.tsv`), from both grouped and ungrouped
runs where PANTRY has both, de-duplicated by phenotype and variant. `chrom` and
`pos` are parsed from PANTRY's `CHROM:POS:REF:ALT` variant IDs; effects are not
flipped. The crosswalk uses GUSTO's own `event_map`:

| Outcome | Rule |
|---|---|
| event, verified | the PANTRY transcript (alt_TSS representative or isoform) is a GUSTO transcript of the same gene; linked to its 50 bp TSS event |
| gene, verified | otherwise the PANTRY gene matches a GUSTO gene exactly, or by unversioned Ensembl ID when unique (`mapping_basis` records which) |
| ambiguous / unmapped | the unversioned ID matches several GUSTO genes / none |

PANTRY's novel `HPLRT_*` transcripts are not GUSTO ESPRESSO transcripts, so
their phenotypes fall back to gene context.

```bash
use_cohort GUSTO
QTL_DIR="$COHORT_DIR/private/qtl_context_v1"     # must not exist yet
python3 bin/prepare_gusto_qtl_context.py \
  --pantry-dir /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/qtl_results \
  --event-map "$(grep -m1 '^  event_map:' "$CONFIG" | awk '{print $2}')" \
  --out-dir "$QTL_DIR"
```

It prints rows, phenotypes and mapping outcome per table (also in
`summary.tsv`). Then set, under `regulatory.human_context` in the private YAML:

```yaml
  human_context:
    require_qtl: true
    qtl_tables: <QTL_DIR>/qtl_inputs.tsv
    phenotype_map: <QTL_DIR>/phenotype_map.tsv
```

and submit `human-context` after `te-tss` completes. `--qval`, `--ancestries`
and `--modalities` change the selection; write each variant to a new
`--out-dir`.

### GUSTO: ENCODE placenta cCREs as chromatin inputs

ENCODE annotation ENCSR057ORE (cCREs in placenta, female embryo 113 days,
GRCh38; DNase-seq, H3K27ac, H3K4me3) is distributed as a bigBed that lists
every registry cCRE with its activity group in this tissue.
`bin/prepare_encode_ccre_beds.py` reads the bigBed directly, drops `Low-DNase`
(inactive in this tissue) and writes one BED per active group (`PLS`, `pELS`,
`dELS`, `DNase-H3K4me3`, `DNase-only`) plus a `beds` manifest, so each group
is reported as its own evidence type alongside the hTSC CUT&Tag peaks.

```bash
use_cohort GUSTO
ENC="$COHORT_DIR/reference/regulatory_resources/encode"
# ENCFF603ENH.bigBed from the ENCODE portal (md5 b719530058318231d99fa45c8ed05264)
md5sum "$ENC/ENCFF603ENH.bigBed"
python3 bin/prepare_encode_ccre_beds.py --bigbed "$ENC/ENCFF603ENH.bigBed" \
  --accession ENCSR057ORE --out-dir "$ENC/ENCSR057ORE_v1"
```

Then add `beds: <ENC>/ENCSR057ORE_v1/beds_manifest.tsv` under
`regulatory.human_context`.

## Fish repeat discovery is separately runnable

`fish-repeats` runs BuildDatabase/RepeatModeler (discovery) and then
RepeatMasker (masking) through `apptainer run` on the configured TE Tools image
(`/risapps/apptainer/repo/tetools/20240523/tetools-202405023.sif`;
RepeatModeler 2.0.5, RepeatMasker 4.1.6). Set `regulatory.fish.cpus` to the
allocated slots (at least four; the RMBlast masking worker uses four cores per
parallel job) and configure `lsf.stages.fish_repeats` with an approved queue,
CPU count, memory and wall time. Discovery requires a completed
`regulatory-resources` receipt; masking requires completed discovery.

### Submit

```bash
# First finish regulatory-resources; no curated insertion track is presumed.
bin/submit_pipeline.sh "$CONFIG" fish-repeats
bin/submit_pipeline.sh "$CONFIG" fish-repeats --submit
# Explicit recovery after inspecting a failed discovery log:
: "${RECOVERY_DIR:?Set RECOVERY_DIR to the actual existing RM_* directory}"
bin/submit_pipeline.sh "$CONFIG" fish-repeats \
  --resume "$RECOVERY_DIR" --submit
```

Two separate jobs use strict discovery `done()` → masking dependencies.
The configuration generator fills missing `lsf.stages.fish_repeats` settings
with an initial `long` queue, 4 CPUs, 64 GB memory/reservation, and a 120-hour
wall limit. It retains customized private values. These settings are provisional,
not measurements of full-genome time, memory or scratch demand, and must be
reviewed against site-approved LSF limits before submission. The default
`short` queue is not reused. No production fish job was run here.

Discovery stages a separate unmasked reference, verifies the expected assembly
length/contig count and sequence identity to an existing primary FASTA, and
runs BuildDatabase → RepeatModeler (`-threads`, `-srand`, `-LTRStruct`).
Masking uses its family library with RMBlast; RepeatMasker `-pa=floor(cpus/4)`
accounts for four cores per worker. Unknown families remain unknown. Tool
versions, actual image/database identities and seeds are logged. The
image and primary FASTA are never modified. No native route is offered until
such an installation is explicitly inspected and supported.
Only an exact RM_* subdirectory of this run's work directory, with an existing
running discovery receipt, can be recovered. There is no general masking recovery flag.
The resulting track represents reference insertions, **not** Line B insertion
polymorphism, TE expression or a completed evolutionary homology analysis.

## Endpoints and outputs

### #1: paired within-gene effects, conditional on first donor

Exact transcript TSS bases determine TE membership; first-exon and ±500/100 bp
strand-aware windows are separate descriptive annotations. Mixed or unknown
events do not become pure TE/nonTE events. Pairs have the same gene,
chromosome, strand and known first donor; distinct 50-bp events must be >50 bp
apart with nonzero distinguishing sequence. Single-exon/unknown donors are
not silently matched.

Use a single shared callable locus, selected only by biological-unit coverage
and depth. One representative per family is kept. The two event responses are
log2((maternal+0.5)/(paternal+0.5)) minus the parental dosage offset. Human
models include cis orientation and configured centered covariates; reciprocal
models include cross direction. Frozen human family/orientation minima and the
configured public callable-unit minimum still apply. One-direction public
crosses explicitly report cis/POE confounding. No new 30-gene power gate is
introduced.

Report **signed difference** beta_TE − beta_nonTE and separately **magnitude
difference** |beta_TE| − |beta_nonTE|. Joint biological-cluster bootstrap keeps
the paired covariance and preserves reciprocal strata. Signed tests use
HC3 (independent clusters) or CR1 cluster sandwich errors; magnitude tests use
a centered paired bootstrap. Magnitude inference is exploratory and
nonregular near zero, not a definitive test that TEs cause stronger imprinting.
Rank/cluster/leverage failures and <90% valid bootstrap replicates are explicit.
BH families are separate for selected-pair signed/magnitude tests, and for one
callability-selected representative per gene. Other pairs/statuses are retained.
Cohort summaries use gene-level resampling without pooling species or claiming
independent replication from several analyses of the same allocations.

Outputs: `transcript_TE_annotation.tsv`, `event_TE_annotation.tsv`,
`TE_annotation_coverage.tsv`, `te_tss_event_pairs.tsv`,
`te_tss_gene_summary.tsv`, `te_tss_cohort_summary.tsv`,
`te_tss_provenance.tsv`. Empty/untestable results remain auditable.

### #5: human evidence integration, not causal inference

`human_regulatory_context.tsv` retains event TE labels, the representative
paired POE fit and its status, event-level/gene-level QTL counts, exact-TSS
peak overlaps, maximum directly linked PIP, and variants overlapping the same
TE insertion as the event TSS. Long tables retain all mapped QTL/peak evidence;
coverage tables distinguish unavailable data from assessed tables with no
observed link. An assessed QTL resource is **not** proof every event was
measured; inspect the phenotype mapping coverage.
Variant overlaps use the reported 1-based position, not the full span of an
indel; structural/insertion-variant mapping needs a separate explicit adapter.

Additional files: `human_qtl_links.tsv.gz`, `human_peak_overlaps.tsv.gz`,
`human_mapping_coverage.tsv`, `human_resource_coverage.tsv`,
`human_context_provenance.tsv`. Missing QTL is not converted to zero evidence
unless chromatin-only scope was explicitly requested. There is no PIP or POE
significance candidate cutoff, new fine-mapping, new colocalization fitting,
parent-of-origin QTL claim or developmental trajectory model. Shared GUSTO
samples and reused annotations are not independent validation.
