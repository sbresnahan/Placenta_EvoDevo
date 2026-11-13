# All-SNP baseline results

`baseline_nb_v1` is a separate result namespace. Legacy primary, additional_v1
and additional_v2 POE files are not valid inputs to its model or TE contrasts.

## Inputs

Use completed Ornaments `allele_counts.txt` and `tpms.txt` for every configured
library, the original annotation/exon mapping and validated parental assignment.
REF/ALT must remain REF/ALT. Undefined parental origin remains missing.

Transcript records are keyed by sample, transcript and full genomic SNP
`CHROM:POS:REF:ALT`. Gene records are keyed by sample, gene and SNP. At each
gene-SNP, each parental count is the sum of transcript count times that
sample's within-gene TPM usage. Weights are not renormalized among the subset
of isoforms carrying a SNP. Keep fractional values until model-input rounding.

Missing feature-SNP rows have zero contributed counts, as in the reference
matrix construction. Explicit unknown parental assignments remain missing
and are excluded as paired maternal/paternal observations. A zero is not
evidence that a sample is heterozygous or informative.

## Model

Human baseline

```r
count ~ allele + sex + ancestry + (1 | SNPpos) + (1 | sample_id)
```

Reciprocal baseline

```r
count ~ allele * cross + available_sex + available_ancestry +
  (1 | SNPpos) + (1 | sample_id)
```

Both use glmmTMB `nbinom2`, a shared maternal/paternal sample size-factor
offset and reference optimizer settings. Sex and ancestry are mandatory
for the literal GUSTO model. Public constant or unavailable factors are
omitted explicitly, not replaced by invented labels. Missing values in a
requested nonconstant factor cause an explicit ineligible result.

Reference-style gates are 3 positive-count informative SNPs and 50 total
counts after fit-input rounding. GUSTO requires 15 families, with one imbalance-blind representative
per family. The public exploratory gate is 4 biological units/clusters and
2 per cross direction. Unknown cluster IDs are flagged. They are never
silently certified independent by assigning sample IDs.

## Main outputs

- `gene_POE.tsv`, `transcript_POE.tsv`. One row per feature, including failures.
- `*_inventory.tsv.gz`. Every contributed feature-SNP, positive sample count,
  total expected count and count of explicit unknown parental records.
- `*_availability.tsv`. Pre-fit feature coverage and count filters.
- `transcript_TE_annotation.tsv`. Exact TSS TE labels for all reference transcripts.
- `transcript_TE_POE.tsv`. Joined annotations and baseline fits.
- `within_gene_jobs.tsv`. Mixed genes eligible for the own-marker-set comparison.
- `within_gene_contrasts.tsv.gz`. Each TE/non-TE transcript contrast and one
  equal-transcript-weight `class_mean` row per gene.
- `within_gene_cohort_summary.tsv`. Equal-gene-weight summaries, not isoform-pair replication.
- `between_genes/`. Gene-level weighted comparisons and balance/placebo diagnostics.
- `transcript_between_gene_*.tsv`. Distinct descriptive transcript-class estimand.
- `context/`. New baseline transcript estimates with separately qualified xQTL/chromatin context.
- `nominal_transcript_candidates.tsv`. New nominal candidates, not the previous eleven.
- `TE_availability.tsv`, `input_audit_summary.tsv`, `provenance.tsv`,
  `package_versions.tsv`, `baseline_settings.tsv`, `report_summary.tsv`.

### Effect fields

`beta` is the raw ln(paternal/maternal) effect. For reciprocal models it is an
equal-cross average, with variance calculated from the full fixed-effect
covariance matrix. `beta_reference_cross` is retained separately.

`beta_adjusted = beta - ln(dosage_paternal/dosage_maternal)`. This centers
maize at its expected 1:2 paternal:maternal ratio without changing the fitted
formula. `p_raw` tests raw beta against zero. `p`, `ci_low`, `ci_high` and
`log2_ratio_p_over_m` refer to the dosage-adjusted effect. For diploid tissues
raw and adjusted estimates coincide.

`beta_cross0` and `beta_cross1` are dosage-adjusted effects for focal-mother
and alternate-mother crosses. `interaction_p` tests their difference.
`cross_main_p` concerns abundance, not parental allelic bias.
`not_strictly_POE` is true when the interaction raw p is below 0.05.
Non-significance does not prove strict POE. All flags remain in the primary
analyses, with separate no-detected-interaction sensitivity outputs.

### Status fields

`eligible_NB` and `exploratory_NB` mean a fit passed numerical checks, not that
POE is established. Other rows specify insufficient SNPs, samples, biological
clusters, reciprocal coverage, counts, rank, convergence or Hessian quality.
`pdHess`, optimizer status, warnings and actual formula are retained.

`q` and `interaction_q` are BH values retained for audit. They are not feature
inclusion gates. Raw p < 0.05 is a nominal reference, not false-discovery control.

Within-gene bootstrap CIs require at least 90% valid refits for the same
transcript set. Absolute-effect tests near zero and all small-cohort results
remain exploratory. A class contrast is TE minus nonTE on the ln scale.

The adjusted gene-level exposure comparison retains the existing definition:
any multi-exon TSS event in a TE versus no such TSS in a TE or unknown repeat.
Its magnitude outcome is unbiased squared POE, `beta_adjusted^2 - se^2`,
with squared-ln units. This is distinct from the absolute-beta magnitude
contrasts in the within-gene and descriptive transcript-class outputs.

## Privacy and completeness

Individual-level records, metadata, size factors, model batches and bootstrap
draws live under `private/`. They are not included by `--baseline-nb` export.
Only allowlisted aggregate tables are exported. Identifier columns cause
export failure. File hashes, versioned receipts and source/input fingerprints
prevent silent reuse of old or changed products.

Missing repeat resources produce `not_assessed` labels and explicit unavailable
TE comparisons, not negative TE results. An explicitly unavailable optional
resource can be archived. A missing or failed required baseline stage cannot.
