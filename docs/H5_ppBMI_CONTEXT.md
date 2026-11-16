# Baseline and ppBMI-dependent QTL context

## Input contract

Producer attribution: `sbresnahan/Placenta_xQTL_atlas`,
revision `6f7435a334ed2d9d7f36ecd737c46efede46dd87`, directory
`07_gxe_mapping/`. Relevant producers are `50_build_gxe_inputs.py`,
`52_gxe_scan.py`, `54_tier1_meta.py`, `54_tier2_stratified.py` and
`57_aggregate_gxe.py`. Source:
https://github.com/sbresnahan/Placenta_xQTL_atlas/tree/6f7435a334ed2d9d7f36ecd737c46efede46dd87/07_gxe_mapping

Read the three authoritative aggregate filenames in
`manifests/gxe_tables.template.tsv`. Avoid duplicate per-modality files,
significant-only tables, parquet shards and exhaustive Tier2 variant scans.
Every manifest row declares GRCh38, one-based positions and
`CHROM:POS:REF:ALT`. Supplied coordinate and allele columns are validated.
Effects and alleles retain the producer's orientation.

For the established Seadragon deployment, run `bin/setup_external_ppbmi.R`
through `bash bin/Rscript_sif`, supplying the project root. The populated
`config/external_ppbmi.yaml` selects the external result root and actual
producer config. The helper derives type from `inputs/exposure_manifest.tsv`
and transformation/range from the actual producer config. Supplied records
declare continuous ppBMI, no transformation, source limits 14–70, 417 nonmissing
observations and external `enabled=1`. The reviewed input builder scales
continuous exposures using pooled retained ancestries and population SD
(`ddof=0`). [1] The versioned exposure manifest records the mean and SD.

Local `enabled=TRUE` means include these existing external results. It is not
evidence of a prior ppBMI analysis in Placenta_EvoDevo. A generic checked-in
producer default cannot establish a production run's enabling state. The
small-file receipt records supplied metadata and source sizes, not an immutable
production config snapshot. The adapter does not open participant-level
exposures, covariates or genotype inputs.

The setup helper writes versioned private manifests and an operator copy,
preserving existing H1 settings, optional exposure settings and QTL/BED overrides.
The H1 count-provenance and pipeline annotations are optional: missing labels
do not block configuration generation or H1 initialization. No new pipeline
version audit is required for unchanged counts and no H1 provenance is invented.
Genuinely required unresolved settings are still reported. All output conflicts are checked
before writes. Identical reruns reuse outputs; differing files are not replaced.
After completed setup, run the generator without rerunning setup. Previously
generated operators containing the two shipped H1 annotation placeholders can be
used directly; they and their historical setup receipts need not be rewritten.
The generic `ppBMI_exposure.template.tsv` remains a schema guide, not a completed
production declaration.

The supplied summary contains 18 merged groups labelled exactly `ppBMI` as well
as chromosome-labelled entries. The aggregator can include both. [2] The adapter's
exact exposure filter must remain: never strip chromosome suffixes or combine
shard rows with merged rows. Setup reads aggregate headers only. The later
context stage still validates actual rows and mapping coverage.

## Statistical meaning

- Tier1: ancestry-specific full-feature cis scans. `b_gi` and `b_gi_se`
  describe genotype × ppBMI interaction; `b_g` and `b_e` are main effects.
  Variant-level nominal p, feature-level `pval_beta` and feature `qval`
  remain distinct. Primary source support is `qval <= 0.05` within
  ancestry × modality × exposure. All tested features remain available.
- Tier1 synthesis: feature-level ACAT with ancestry-only features retained.
  `meta_variant_id`, `beta_meta` and `se_meta` require concordant lead variants.
  Discordant leads support feature-level synthesis without a shared variant
  effect. Synthesis flags remain separate from ancestry-specific primary support.
- Tier2: upstream coloc/TWAS-prioritized features; meta-p is minimized across
  cis variants before BH. `q_meta` lacks within-cis search calibration in
  this producer and is retained as an exploratory source flag.

The adapter performs no allelic harmonization, fine mapping, credible-set
construction, PIP estimation or independent-replication test.

## Mapping and missingness

The verified phenotype map records event/gene mappings and ambiguous/unmapped
features. Exact-transcript evidence additionally requires the embedded PANTRY
transcript ID to match the native transcript. Other event members retain
event-level scope. Gene mappings retain gene context.

Baseline and `gxe_ppBMI` outputs remain separate. GxE strata retain ancestry,
tier, modality and exposure. Assessed nonsupport is FALSE; untested,
unmapped and missing-statistic support stays NA. Mapping coverage preserves
ambiguous and absent crosswalk states.

Exact TSS–cCRE overlaps and baseline exact counts keep their definitions.
`variant_overlaps_TSS_TE` requires the variant inside an insertion overlapping
that transcript's exact TSS. `variant_overlaps_proximal_only_TE` requires
the variant inside an insertion linked only through the active promoter
window. A variant elsewhere in the window satisfies neither requirement.

## H6/H7 follow-up

Baseline support enrichment and ppBMI support enrichment are separate outputs.
ppBMI follow-up has separate ancestry × modality × tier × support contexts.
Within each context, multiplicity spans support types and modeled exposures
for each module/cohort/feature level/predictor. Uncertainty clusters on gene.
Sparse or separated fits remain unavailable.

ppBMI context neither replaces maternal PFAS nor enters the H7 mediation
covariates. Shared atlas cohorts and annotation selection require provenance
review before any independent-validation claim.

## References

1. https://raw.githubusercontent.com/sbresnahan/Placenta_xQTL_atlas/main/07_gxe_mapping/50_build_gxe_inputs.py. [cited 2026 Nov 16]. Available from: https://raw.githubusercontent.com/sbresnahan/Placenta_xQTL_atlas/main/07_gxe_mapping/50_build_gxe_inputs.py
2. https://raw.githubusercontent.com/sbresnahan/Placenta_xQTL_atlas/main/07_gxe_mapping/57_aggregate_gxe.py. [cited 2026 Nov 16]. Available from: https://raw.githubusercontent.com/sbresnahan/Placenta_xQTL_atlas/main/07_gxe_mapping/57_aggregate_gxe.py
