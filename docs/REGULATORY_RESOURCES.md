# Repeat and regulatory resources

These host-submitted stages prepare the repeat tracks and resource manifests
that the baseline TE stages read (`te-inputs`, `between-genes`, `context`).
They write to `regulatory.output_dir`, which defaults to `regulatory/` beside
`model_gene`; scratch work goes to `regulatory.work_dir` (default
`<paths.scratch>/regulatory_work`). Network transfers and LSF submission run on
the **host**, not in the R image.

## Configure

Copy missing regulatory keys from the species-aware cohort examples into
schema-2 private configs with `python3 bin/sync_regulatory_configs.py --create-missing`
(preview first with `--dry-run`). `config/regulatory.example.yaml` is an
illustrative fragment, not a full config. Receipts track successful completion
and required output presence, not input, script or configuration content
identity; inputs are assumed frozen for an analysis. Paths in manifests are
relative to the cohort config directory, not the manifest directory. Explicit
absolute paths are recommended.

```bash
# From the scripts directory, with CONFIG set to the cohort's private YAML:
: "${CONFIG:?Set CONFIG to your actual configured cohort YAML}"
export SCRIPTS_DIR="$(pwd)"
# Inspect (prepares scripts only); append --submit to submit.
bin/submit_pipeline.sh "$CONFIG" regulatory-resources
bin/submit_pipeline.sh "$CONFIG" regulatory-resources --submit
bin/Rscript_sif poe_regulatory.R "$CONFIG" regulatory-status
```

`poe.R`, `run_pipeline.R` and `poe_regulatory.R` route `regulatory-resources`,
`fish-repeats` and `regulatory-status`. Only `bin/submit_pipeline.sh ... --submit`
executes the returned host submission plan; R never calls `bsub`.
`--after NUMERIC_JOB_ID` adds a strict `done()` dependency. Jobs have separate
stdout/stderr files and scheduler IDs under `<regulatory.output_dir>/jobs`.
A host plan can be executed only once to prevent accidental duplicate jobs.
Submission failure does not cancel jobs already submitted: inspect their IDs.
Existing partial outputs are never overwritten implicitly.

`regulatory-status` reports `regulatory-resources`, `fish-discovery` and
`fish-masking` as `not_run`, `complete_outputs_present` or
`incomplete_or_outputs_missing`.

## Resource acquisition

Species defaults download UCSC hg38/mm10 `rmsk.txt.gz` and its schema, or the
native B73v4 TE GFF3. Fish defaults download only the RefSeq GCF_027474105.1
FASTA/assembly report; fish repeats come from `fish-repeats` below. Human
chromatin peaks are not downloaded; supply them through
`regulatory.human_context.beds`. The cluster TE Tools image is used as
installed; no image or repeat-classification database is downloaded or
installed automatically. No whole-atlas, raw RNA, long-read or
chromosome-wide LD download is implied.

The `regulatory_*example.tsv` files contain headers only, with no invented
local paths or identifiers. Populate them with actual reviewed input rows.
The resources manifest enables local-first replacement by `resource_id`;
omitted defaults remain enabled. Add a URL to download a replacement.
Supported transfer formats are rmsk, gff3, bed, sql, text, fasta_gz and h5_gz.
Plain or gzip BED/GFF are accepted. Missing declared local files are errors,
not fallback downloads. Transfers use bounded retries, partial files and
gzip/HTML format checks. Final manifests record URLs, bytes and retrieval time
only; these do not prove content identity. Existing bad files are not
overwritten. Local overrides require `scope=full_reference` or `subset`; unhit
subset sites are **not_assessed**. Contig aliases are conservative
(`chr`/case normalization or an explicit one-to-one alias table); ambiguous
aliases and off-reference coordinates fail.

Outputs: `resource_manifest.tsv` and, when repeats are declared,
`repeat_intervals.tsv.gz` (normalized insertions with `is_te`, `is_unknown`,
class, family and identifier). Only a `full_reference` repeat scope supports
nonTE labels; a partial reference leaves unhit transcripts `not_assessed`.

## Fish repeat discovery and masking

`fish-repeats` runs BuildDatabase/RepeatModeler (discovery) and then
RepeatMasker (masking) through `apptainer run` on the configured TE Tools image
(`/risapps/apptainer/repo/tetools/20240523/tetools-202405023.sif`;
RepeatModeler 2.0.5, RepeatMasker 4.1.6). Set `regulatory.fish.cpus` to the
allocated slots (at least four; the RMBlast masking worker uses four cores per
parallel job) and configure `lsf.stages.fish_repeats` with an approved queue,
CPU count, memory and wall time. Discovery requires a completed
`regulatory-resources` receipt; masking requires completed discovery.

```bash
bin/submit_pipeline.sh "$CONFIG" fish-repeats
bin/submit_pipeline.sh "$CONFIG" fish-repeats --submit
# Explicit recovery after inspecting a failed discovery log:
: "${RECOVERY_DIR:?Set RECOVERY_DIR to the actual existing RM_* directory}"
bin/submit_pipeline.sh "$CONFIG" fish-repeats \
  --resume "$RECOVERY_DIR" --submit
```

Two separate jobs use strict discovery `done()` → masking dependencies.
The configuration generator fills missing `lsf.stages.fish_repeats` settings
with an initial `long` queue, 4 CPUs, 64 GB memory/reservation and a 120-hour
wall limit, and keeps customized private values. These starting values are not
measurements of full-genome time, memory or scratch demand; review them
against site-approved LSF limits before submission.

Discovery stages a separate unmasked reference, verifies the expected assembly
length/contig count and sequence identity to the primary FASTA, and runs
BuildDatabase → RepeatModeler (`-threads`, `-srand`, `-LTRStruct`). Masking
uses its family library with RMBlast; RepeatMasker `-pa=floor(cpus/4)`
accounts for four cores per worker. Unknown families remain unknown. Tool
versions, image/database identities and seeds are logged. The image and
primary FASTA are never modified. Only an exact RM_* subdirectory of this
run's work directory, with an existing running discovery receipt, can be
recovered; there is no general masking recovery flag.

Outputs: `fish_repeat_intervals.tsv.gz` and `fish_repeat_provenance.tsv`. A
completed masking receipt marks the fish track as a full reference. The track
represents reference insertions, **not** insertion polymorphism, TE
expression or a completed evolutionary homology analysis.

## Human context inputs (GUSTO)

The baseline `context` stage (H5) reads cis-QTL and chromatin inputs declared
under `regulatory.human_context`: authorized local atlas lead/PIP/credible-set
or supplied colocalization tables, an explicitly reviewed phenotype crosswalk,
and optional CAGE/ATAC/chromatin BED inputs. Use the `regulatory_*example.tsv`
manifests. Export atlas outputs as TSV/gzip or supply Parquet with R `arrow`.
No stable public URL or permission to private atlas data is assumed.
Requested missing files block execution with the missing path. Public peaks
do not substitute for QTL. To request chromatin-only output explicitly, set
`require_qtl: false` and omit QTL/crosswalk keys.

QTL manifest kinds: `lead`, `pip`, `credible_set`, `coloc`. Variant rows require
`coordinate_system=one_based`, explicit `chrom`, `pos` and variant ID columns.
The adapter recognizes atlas `pip_all`, `pip_cs`, `cs_index`, `sig_in` and
`lead_snp`/`lead_pos`; all source fields are kept as `source_*`. Optional
`modality_column`, `phenotype_column`, `ancestry_column`, `variant_column`,
`chrom_column`, `pos_column`, `pip_column`, `pip_cs_column`, `cs_index_column`,
`beta_column`, `se_column`, `ref_column`, `alt_column`, `cs_start_column`,
`cs_end_column` map other schemas. All tables require modality/phenotype IDs
and ancestry (or a manifest label `joint`/`unspecified`). `a0/a1` are kept
as allele0/allele1, **not** assumed genomic REF/ALT. No variant-ID parsing,
liftover or allelic effect flipping is guessed.

Crosswalk columns are modality, phenotype_id, gene_id, event_id, mapping_level
(`event`/`gene`) and mapping_status (`verified`/`ambiguous`/`unmapped`). Exact
gene/event agreement is enforced. Nonunique phenotype mappings are withheld,
with an ambiguous mapping status. Gene-level links provide gene context.
Transcript validation requires a verified event/transcript mapping.
Peak manifests require a TRUE/FALSE `selection_used` flag to expose annotation
reuse, including CAGE selection.

### PANTRY placenta cis-QTL tables

`bin/prepare_gusto_qtl_context.py` turns the PANTRY tensorQTL results into the
QTL tables, manifest and phenotype crosswalk above. Per ancestry (EAS, EUR) and
modality it keeps every significant permutation lead (`qval < 0.05` in
`*_top.tsv`) and every conditionally independent signal
(`*_independent_top.tsv`), from both grouped and ungrouped runs where PANTRY
has both, de-duplicated by phenotype and variant. `chrom` and `pos` are parsed
from PANTRY's `CHROM:POS:REF:ALT` variant IDs; effects are not flipped. All
nine PANTRY modalities are supported; `combined` is excluded because it pools
the others. The crosswalk uses GUSTO's own `event_map`:

| Modalities | Link rule |
|---|---|
| `alt_TSS`, `isoforms`, `isoform_expression` | TSS event when the PANTRY transcript is a GUSTO transcript of the same gene; otherwise gene context |
| `expression`, `stability`, `alt_polyA`, `splicing`, `intron_retention`, `RNA_editing` | gene context only; these phenotypes describe a gene, 3′ end, intron or editing site, not a TSS |

Genes match exactly, or by unversioned Ensembl ID when that match is unique;
several matches are `ambiguous`, none `unmapped`. `mapping_basis` records the
rule used for every phenotype. PANTRY's novel `HPLRT_*` transcripts and
`HPLRG_*` genes are not GUSTO ESPRESSO IDs, so they fall back to gene context
or remain unmapped.

```bash
QTL_DIR="$PROJECT/data/GUSTO/private/qtl_context_v1"     # must not exist yet
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

`--qval`, `--ancestries` and `--modalities` change the selection; write each
selection to a new `--out-dir`.

### ENCODE placenta cCREs as chromatin inputs

ENCODE annotation ENCSR057ORE (cCREs in placenta, female embryo 113 days,
GRCh38; DNase-seq, H3K27ac, H3K4me3) is distributed as a bigBed that lists
every registry cCRE with its activity group in this tissue.
`bin/prepare_encode_ccre_beds.py` reads the bigBed directly, drops `Low-DNase`
(inactive in this tissue) and writes one BED per active group (`PLS`, `pELS`,
`dELS`, `DNase-H3K4me3`, `DNase-only`) plus a `beds` manifest, so each group
is reported as its own evidence type. These are the only peak inputs.

```bash
ENC="$PROJECT/data/GUSTO/reference/regulatory_resources/encode"
# Obtain ENCFF603ENH.bigBed from the ENCODE portal.
python3 bin/prepare_encode_ccre_beds.py --bigbed "$ENC/ENCFF603ENH.bigBed" \
  --accession ENCSR057ORE --out-dir "$ENC/ENCSR057ORE"
```

Then add `beds: <ENC>/ENCSR057ORE/beds_manifest.tsv` under
`regulatory.human_context`.

## Between-gene comparison method

The baseline `between-genes` stage complements the within-gene H1 contrast
with all genes that have an estimable baseline gene POE: do genes with at least
one multi-exon transcript TSS window intersecting a TE differ from genes whose
windows are clear of TEs and unknown repeats? Genes whose only repeat overlaps
are unknown repeats are excluded. Gene POE is the baseline all-SNP
`gene_POE.tsv` estimate (ln paternal/maternal, dosage-adjusted). The active
promoter-window radius defines exposure and the shifted-track placebo. The
comparison requires a full-reference repeat track; otherwise it is reported as
unavailable.

Between genes, expression, informativeness, gene structure and genomic
neighbourhood can confound the contrast, so they are adjusted explicitly:

| Confounder | Covariate |
|---|---|
| Depth / precision | log total allelic count, log SE, number of samples |
| SNP informativeness | log number of informative loci with observed samples |
| Gene structure | log number of 50-bp TSS events, log number of transcripts, log gene span |
| Neighbourhood | TE bp fraction within ±100 kb of the gene; promoter GC (±1 kb of the 5'-most TSS) |

Noise inflates |b|, so magnitude uses b^2 - SE^2 (unbiased for the true squared
POE) and the signed contrast uses b; significance calls are never the outcome.
A logistic propensity model on these covariates defines common support
(overlapping propensity range) and ten strata on the TE genes' propensity
quantiles; non-TE genes are weighted by the TE:non-TE ratio in their stratum
(average effect in TE genes). `te_genes_balance.tsv` reports standardized mean
differences before and after weighting; any |SMD| > 0.1 is flagged in the status.

Uncertainty: genes in the same region are not independent, so the bootstrap
resamples 2 Mb genomic blocks and refits the propensity model each time
(`regulatory.bootstrap` replicates; normal-theory p and percentile CI). A
within-stratum permutation of TE labels gives a second p-value. A placebo
circularly shifts the TE track along each chromosome (keeping TE clustering and
density, breaking the TSS–TE link), recomputes exposure and the full analysis,
and reports how often a shifted track gives an effect at least as large.
Reciprocal cohorts also write a sensitivity analysis without genes flagged
`not_strictly_POE` under `between_genes/no_detected_cross_interaction/`.

Outputs in `<baseline output_dir>/between_genes/`: `te_genes_result.tsv` (one
row: contrasts, CIs, three p-values, counts, status), `te_genes_balance.tsv`,
`te_genes_gene_table.tsv` (per-gene exposure, estimate, covariates, propensity,
stratum, weight), `te_genes_placebo.tsv`, `te_genes_provenance.tsv`. Optional
settings under `regulatory.te_genes`: `min_genes_per_group` (5),
`propensity_strata` (10), `block_bp` (2e6), `te_density_flank_bp` (1e5),
`gc_flank_bp` (1000), `permutations` (2000), `placebo_shifts` (200),
`min_shift_bp` (1e6). Interpretation: a within-cohort, between-gene
association conditional on the measured covariates; unmeasured differences
(for example gene age) remain.
