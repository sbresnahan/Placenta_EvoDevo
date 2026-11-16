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

`ppBMI_exposure.template.tsv` requires an enabled ppBMI row, the actual
continuous/binary definition, transform, scale and range, and paths to the
producing configuration and runtime provenance. Template declarations are
incomplete until reviewed. The producer's checked-in config disables ppBMI.
Production enabling and scaling require an explicit declaration. The adapter
does not open participant-level exposures, covariates or genotype inputs.

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
