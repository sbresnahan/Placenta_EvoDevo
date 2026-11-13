# Baseline POE model comparison

## Pinned sources

The patched pipeline is `sbresnahan/Placenta_EvoDevo` at
`ccfa3c825d941b08c6941b92fd7990996bb4cd74`. The requested reference is
`sbresnahan/placenta-poe-pfas` at
`c285f131071becec7347c830d263e78a88ac45de`. Its negative-binomial model
and gene aggregation code are the implementation targets, not a new
independent scientific validation. [1, 2]

The previously supplied results archive is a legacy output snapshot. It
contains aggregate/sentinel outputs, not the private transcript-SNP products
required for this rerun. Consequently no real cohort-wide revised model was
fitted in this patch-development session. [3]

## What changes

| Aspect | Existing selected-SNP/legacy route | New `baseline_nb_v1` route |
|---|---|---|
| POE inputs | Legacy gene/event fits and selected/shared-marker comparisons | Reconstructed parental transcript-SNP counts from raw completed quantification |
| Gene counts | Existing pipeline aggregation | Within-sample transcript-TPM-weighted gene-SNP union |
| Transcript unit | Event/promoter and selected isoform comparisons | One model for each individual transcript, with its own informative SNP set |
| Human model | Legacy allelic log-ratio route | Literal reference NB mixed model, sex and ancestry required |
| Reciprocal model | Legacy route | NB model with `allele * cross`, cross-specific and equal-cross effects |
| Shared SNP gate | Used in the superseded selected-SNP comparison | Removed |
| Same first donor gate | Used in the superseded comparison | Removed |
| TE outcome | Legacy model summaries | New all-SNP gene and individual-transcript effects only |
| Failure behavior | Legacy outputs may exist | Explicit statuses; no fallback to legacy output |
| Export | Existing report archive route | Separate strict `--baseline-nb` allowlist |

The old pipeline entry points and files are retained for reproducibility.
The patch adds an opt-in analysis rather than silently relabeling old fits.

## Literal human reference alignment

The reference fits `glmmTMB::nbinom2` with
`count ~ allele + sex + ancestry + (1|SNPpos) + (1|sample_id)`,
a log median-of-ratios size-factor offset, and `nlminb` with relative
tolerance `1e-5`, at most 300 iterations and 300 evaluations.
Expected counts are rounded at model input. Maternal is the reference
allele, so the allele coefficient is ln(paternal/maternal). The minimum
reference gates are three informative SNPs, fifteen samples and fifty
counts. The reference formula does not explicitly adjust cis orientation;
the new human route deliberately preserves that choice. [1]

The numerical regression check compares the patch directly with that formula
on a software fixture using the same R packages and input. Agreement tests
implementation equivalence under those conditions, not cross-package/version
identity, parent-assignment correctness on private samples, or biological POE.

The pipeline uses its existing REF/ALT parental-assignment path. It does not
copy the reference repository's L/R count-column assumptions. Undefined
parental origin remains missing; either missing parent removes the paired
observation instead of creating artificial allelic imbalance.

## Gene aggregation

The reference calculates each transcript's fraction of the gene's total
transcript TPM for a sample, uses equal fractions if all TPMs are zero, and
sums weighted transcript-level parental counts at each gene-SNP.
The patch follows this weighting, including the union of transcript SNPs.
Weights are not renormalized among only the SNP-carrying transcripts. [2]

For transcript t, gene g, sample s, and SNP k:

`w(t,s) = TPM(t,s) / sum(TPM(g,s))`

`count(g,k,s) = sum_{t carrying k} w(t,s) * count(t,k,s)`

Fractional values are retained through aggregation and rounded only when
forming the NB model input. This intentionally applies usage weighting even
though Ornaments expected counts have already been allocated. An unreadable
TPM file is an error, not justification for substituting equal weights.
Only a genuinely all-zero gene TPM receives equal weights.

## Approved extensions beyond the human reference

### Reciprocal crosses

Public cohorts use
`count ~ allele * cross + available_covariates + (1|SNPpos) + (1|sample_id)`.
The two levels are focal mother and alternate mother. The overall log effect
is the equal-cross contrast `b_allele + 0.5*b_interaction`, with its full
fixed-effect covariance. This is deliberately not the reference-cross allele
coefficient mislabeled as a global POE estimate.

The outputs retain the reference-cross coefficient, both cross-specific
effects and SEs, the allele-by-cross interaction, and cross main effect.
Nominal interaction p < 0.05 sets `not_strictly_POE=TRUE`.
Cross main effects concern abundance, not cross-dependent allelic bias.
A nonsignificant interaction is not proof of strict POE, especially with
two or three units per direction.

The raw fitted coefficient is retained. Dosage adjustment is the derived
contrast `beta_adjusted = beta - log(paternal_dosage/maternal_dosage)`.
For 2 maternal:1 paternal endosperm the log null is `log(0.5)`.
This changes the interpreted null and derived p value, not the fitted model.
Positive adjusted beta means paternal bias.

### Biological replication and small cohorts

Human preparation chooses one child per family globally by callability,
then depth, then sample ID, never by imbalance. Public exploratory models
require at least four units/clusters, with at least two per reciprocal
direction. SNPs and libraries are not biological replicates.
Missing cluster identities are explicitly independence-unverified; generated
computational labels do not establish independent organisms.

The old archive records 192 GUSTO libraries/families, but lacks the required
model covariates. It records only five fish units with both reciprocal
directions, mixed stages and forced unadjusted pooling, and contains no fish
TE-TSS analysis. These are reasons to obtain missing inputs and preserve
cohort-specific limitations, not to fabricate a revised estimate. [3]

### TE comparisons

Exact transcript TSS overlap determines TE versus non-TE classification.
Unknown repeats and unassessed references remain distinct. Within a gene,
each eligible transcript keeps its own SNPs. The main comparison includes
all estimable mixed TE/non-TE genes, not only nominally significant POE
transcripts.

Each gene exports all TE/non-TE transcript-pair contrasts and the difference
between equal-transcript class means, for both signed adjusted beta and
absolute adjusted beta. Joint biological-cluster bootstrap refits retain
within-gene dependence and reciprocal strata. The default is 2000 draws,
with at least 90% valid complete-transcript draws required for uncertainty.
Bootstrap summaries in very small cohorts remain exploratory.

Across genes, the existing adjusted gene-level exposure comparison is
reused with new gene NB outcomes and new depth inventories. A separate
descriptive transcript-class comparison uses equal gene weights, not one
replicate per transcript. Cross-dependent flags and no-detected-interaction
sensitivity results are retained; they do not silently replace the main set.
Missing full repeat references make the gene-level exposure comparison
unavailable, not negative.

The adjusted gene-level exposure remains any TE-overlapping multi-exon
TSS event versus no such TSS in a TE or unknown repeat. Its magnitude
outcome is `beta_adjusted^2 - SE^2`, in squared-ln units. This differs
from the absolute-beta contrasts in within-gene and descriptive
all-transcript class comparisons. They must not be pooled as one estimand.

### Context and handoff

Human xQTL/chromatin overlap is annotation context, not evidence of
colocalization, TE causality, or strict imprinting. Event/gene mappings are
distinguished from matching individual transcript identifiers. Missing
resources remain unassessed.

The final archive includes versioned gene/transcript results, availability,
cross flags, comparisons, context, nominal candidates, provenance and
checksums. It excludes individual-level inputs, private metadata, model
objects, bootstrap records and old sentinel endpoints. The user must run
this workflow on their private deployment and return that archive before
any revised biological claims can be assessed.

## References

1. R/poe_nbmm.R. [cited 2026 Nov 13]. Available from: https://github.com/sbresnahan/placenta-poe-pfas/blob/c285f131071becec7347c830d263e78a88ac45de/R/poe_nbmm.R
2. 06_feature_matrices/ornaments_summarizeToGene.R. [cited 2026 Nov 13]. Available from: https://github.com/sbresnahan/placenta-poe-pfas/blob/c285f131071becec7347c830d263e78a88ac45de/06_feature_matrices/ornaments_summarizeToGene.R
3. results_20261113_1251.tar.gz. File provided by the user: results_20261113_1251.tar.gz.
