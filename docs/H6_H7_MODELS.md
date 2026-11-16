# H6/H7 model specification

## Status and scope

These modules specify within-cohort H6/H7 association and mediation models.
Production execution requires cohort-native counts, reviewed biological clusters
and measured resource pilots.

H6 tests whether **sample-varying RNA abundance of the matching TE family** is
associated with continuous POE beta. H7 tests **maternal PFAS → matching-family
TE RNA → signed POE beta**, with mediation primary and modification optional.
H7 is restricted to GUSTO. No baseline POE, H1, H5, H6, total-effect or individual
path significance is an entry criterion.

POE beta is the natural-log paternal/maternal expression ratio minus
`log(paternal_dosage/maternal_dosage)`. Positive means paternal bias, negative
means maternal bias. A positive exposure slope shifts expression toward paternal
expression; it need not increase the absolute magnitude of bias.

## Inputs and TE exposure

Use the existing private prepared gene/transcript SNP-count batches. Keep
fractional expected counts until rounding at the NB fit boundary. Reuse full
genomic SNP keys, offsets, parental orientation and biological-family definitions.
Require an explicit RNA/unit/family/SubjectID crosswalk.

Native TE loci are 1-based inclusive intervals with separate locus, family,
repeat-name/subfamily and class fields. Insertions intersecting the inclusive
native TSS ±3000-bp window determine feature-family pairs. Windows are clipped
to native contigs. Radius-0 and radius-1000 runs are separate sensitivities.
Each transcript-family target retains insertion IDs, exact/proximal-only flags,
radius and nearest distance. Gene tests deduplicate families across transcripts.
Cross-taxon family equivalence requires an explicit mapping.

The default predictor excludes all loci overlapping the target gene span ±1 kb,
sums the remaining matching-family locus counts, divides by mapped fragments,
then uses `log2(CPM + 0.5)` and cohort z-scores. Scaling uses quantified,
crosswalk-matched samples **before feature-specific allele-count exclusions**.
Store centers/scales. Require ≥1 CPM in ≥20% of matched samples and ≥3 samples,
with positive variance. Sensitivities: inclusive, intergenic-only, and
unique-mapper target-excluded. Family-only imports cannot supply the primary
local-excluded predictor.

Strict unique-mapper counts must derive from genomic NH=1 alignments and a
separate locus-quantification report, declared in QC as `genomic_NH1`.
Telescope's ordinary annotation-unique count is retained separately and never
silently substituted for genomic mapping uniqueness.

The measurement is **TE-sequence-containing RNA abundance**, not proof of
autonomous transcription. Local exclusion reduces direct self-inclusion; it
does not remove all read-assignment uncertainty or establish causal ordering.
Locus models use disjoint segments: conflicting families/strands are excluded,
and the shortest same-family locus owns a nested segment. The original catalogue
and GTF/QC retain the provenance of this conservative annotation choice.

## H6: matching TE–POE association

GUSTO:

```r
count ~ allele * E + sex + ancestry + (1 | SNPpos) + (1 | sample_id)
```

Reciprocal cohorts:

```r
count ~ allele * cross + allele * E_within_cross +
  available_sex_and_ancestry + (1 | SNPpos) + (1 | sample_id)
```

Both use NB2, existing count offsets, `nlminb`, one thread and existing optimizer
tolerances. Public E is centered within reciprocal direction on the informative
sample set. Between-direction exposure differences alone cannot identify H6.
The common within-direction slope remains in cohort TE-z units; there is no
second within-direction SD standardization.

Gates: ≥3 informative SNPs; ≥50 total rounded counts; GUSTO ≥15 informative
families; public ≥4 biological units with ≥2 per direction; varying exposure
(at least two distinct values), varying exposure within both reciprocal
directions, full fixed-effect rank, biological units greater than sample-design
rank, condition number ≤1e8, convergence code zero, positive-definite Hessian and
finite coefficients/covariance. Human sex/ancestry must remain estimable.
All identifiable small public cohorts are attempted; no failed model is
automatically simplified. Four-to-six-unit estimates remain exploratory:
asymptotic Wald intervals may be poorly calibrated at that size.

At exposure x, the human signed beta is `beta_A + beta_AE*x - dosage_null`.
Public predictions equally average the two reciprocal-direction allele
contrasts. Derive p10/p50/p90 predictions and signed p90–p10 change from the
full covariance. Absolute-magnitude change is `abs(beta90)-abs(beta10)`,
not `abs(beta90-beta10)`. A simultaneous two-endpoint confidence region gives a
conservative magnitude-change interval. Fitted and individually interval-supported
direction reversals are distinct output flags; the latter is not a separate
multiple-testing-adjusted discovery claim.

## H7: primary additive mediation

Use the supplied maternal z-scores unchanged:
`PFBA_mat, PFOA_mat, PFNA_mat, PFDA_mat, PFUnDA_mat, PFBS_mat, PFHxS_mat, PFOS_mat`.
An explicit private preprocessing dictionary is mandatory. Each compound uses
its available complete cases; optional common-panel mode excludes any participant
missing any of the eight. Both paths use exactly the outcome-informative sample
set for that feature/family/compound.

```r
# E: supplied PFAS z-score; M: matching-family TE predictor
M ~ E + sex + ancestry
count ~ allele * (E + M) + sex + ancestry + (1 | SNPpos) + (1 | sample_id)
```

The mediator is linear with family-cluster HC3 uncertainty and a t reference
using `n_families - mediator_design_rank` degrees of freedom. In the count model,
`b = allele:M`, `c_prime = allele:E`; `a` is the mediator E slope. Signed beta
has no E×M interaction, so per one supplied PFAS-z unit the indirect effect is
`a*b`, the direct effect is `c_prime`, and the modeled total is `c_prime+a*b`.
This decomposition concerns the modeled **log allelic ratio**, not an effect on
total RNA counts. The separately fitted PFAS-only total model is not forced
to equal the model-implied total, especially if its sample set differs.

Primary indirect-path p-value: `max(p_a,p_b)`, an intersection-union test of
the null that at least one path is zero. It needs no path-test independence.
For the 95% product interval, use two 97.5% coefficient intervals and extrema
of all four endpoint products. A three-coefficient simultaneous rectangle
gives a conservative interval for `c_prime+a*b`. Coverage remains contingent on
the underlying coefficient-interval approximations.

Empirical-residual standardization contrasts PFAS p25 and p75:
predict `M(e0)`/`M(e1)` from the mediator mean plus each observed residual,
then form `beta(e0,M(e0))`, `beta(e1,M(e0))`, `beta(e1,M(e1))`.
Direct and indirect signed changes add to total; the corresponding differences
of **absolute beta** also telescope to total magnitude change. Absolute indirect
effect is not `abs(a*b)`. Report scenario values, reversals and the fraction of
mediator counterfactuals outside observed mediator support.

Secondary paired bootstrap: resample biological families, retain every SNP row,
relabel repeated blocks, refit both models, and repeat the standardization.
Default: 2,000 draws for the 20 lowest adjusted indirect-path p-values, with
deterministic tie-breaking. Intervals require ≥90% valid paired refits.
Exposure quantiles and the empirical standardization population are re-estimated
within each draw. The bootstrap is selected-result uncertainty checking, not
independent replication, post-selection correction or a replacement for primary
all-pair inference. `bootstrap_all: true` requires a measured cost pilot.

## Assumptions and secondary modification

Retain sex and ancestry; **assume no remaining confounding**, as requested.
No confounder selection or confounding-sensitivity module is included. A causal
interpretation additionally assumes PFAS precedes TE, TE precedes POE, consistency,
adequate support, correct functional forms, and no exposure-induced
mediator–outcome confounding. These assumptions are not demonstrated by the fit.
TE and POE from the same RNA sample do not establish temporal order. Natural
direct/indirect effects require these identification conditions. [1, 2]

Optional modification expands the outcome to `allele * E * M`, retaining all
lower-order terms, with ≥30 families, ≥5 families in each E/M median support
quadrant, full rank and conditioning ≤1e8. The exposures stay continuous.
Standardization, not a simple product, defines its effect decomposition.
The modification interaction p-value is secondary and is not indirect-path
evidence.

## Multiplicity, diagnostics and follow-up

BH discovery threshold: 0.05. H6 correction spans feature-family tests within
cohort/level/predictor. `ex_pool_h6()` performs an additional cross-cohort
correction without pooling coefficients; the pooled q-value remains unavailable
in a single-cohort run. H7 indirect correction spans all eight PFAS and tested
feature-family pairs within level/predictor. Direct, standalone total and
modification p-values have separate families. Per-compound q-values are secondary.
Unavailable planned rows remain visible; all planned rows count in the new
module's conservative BH denominator.

Diagnostics include informative counts/units/SNPs, residual unit df, exposure
support, conditioning, convergence/Hessian, PFAS–TE correlation, Pearson residual
RMS/extremes, maximum family count/residual concentration, and mediator residual
RMS, leverage and Cook distance. A full leave-one-family-out influence analysis
requires additional fits. Count-model Wald tests use the specified NB model covariance;
family bootstrap provides a secondary dependence/uncertainty check for selected
pairs. Technical replicates cannot be treated as independent units.

Freeze functional memberships and versions before new fits. One gene-level
score Bonferroni-combines all planned within-gene transcript/family tests.
Competitive enrichment compares unsigned evidence ranks against matched-gene
permutations, stratified on counts, biological units, SNPs, length and number of
tests; sparse strata are deterministically coarsened without set membership.
Require ≥10 tested genes, set/background ≥5 and set ≤2,000. Report coarsening.
Keep cohort-specific native gene universes; orthology mappings must be documented.

H5 follow-up is GUSTO-only and retains exact-transcript QTL, 50-bp TSS-event QTL,
lead-in-same-TSS-TE and exact-TSS cCRE definitions. Unassessed annotations stay
missing. Logistic discovery-label enrichment adjusts for counts/SNPs/units with
gene-cluster HC1 uncertainty; sparse/separated models are flagged. Continuous
evidence summaries accompany it. These feature-detectability adjustments are
not extra causal confounder adjustment in H7. Neither functional nor H5 enrichment
validates mediation.

## Verified implementation resources

The strand-aware Telescope source is pinned at
`26a6c38df6d143387799edc2f3e9f976505f093d`, declaring version `1.0.3.1`.
The compatible report
uses named locus/count columns. STAR 2.7.11b is a verified release. [3–6]

The frozen **category selection** includes verified GO terms for amino-acid/
carbohydrate transport, lipid/cholesterol metabolism, hormone signaling,
placenta development and immune regulation; it does not claim to provide a
downloaded annotation release or completed memberships. Reactome transport and
lipid-metabolism pathway records provide complementary human categories.
Geneimprint supplies an independent, status-specific imprinting catalogue.
Actual species-native memberships, source versions and mapping coverage must
be frozen privately before fitting. [7–10]

## References

1. Valeri L, Vanderweele T. Mediation analysis allowing for exposure-mediator interactions and causal interpretation: theoretical assumptions and implementation with SAS and SPSS macros. Psychological methods. 2013. doi:10.1037/a0031034. PMCID: PMC3659198.
2. Vanderweele T, Vansteelandt S, Robins J. Effect decomposition in the presence of an exposure-induced mediator-outcome confounder. Epidemiology (Cambridge, Mass.). 2014. doi:10.1097/ede.0000000000000034. PMCID: PMC4214081.
3. https://raw.githubusercontent.com/mlbendall/telescope/v1.0.3/telescope/utils/model.py. [cited 2026 Nov 15]. Available from: https://raw.githubusercontent.com/mlbendall/telescope/v1.0.3/telescope/utils/model.py
4. 2.7.11b. [cited 2026 Nov 15]. Available from: https://api.github.com/repos/alexdobin/STAR/releases/tags/2.7.11b
5. 26a6c38 fixed logic in determining fragment strands. [cited 2026 Nov 15]. Available from: https://github.com/mlbendall/telescope/commit/26a6c38df6d143387799edc2f3e9f976505f093d
6. https://raw.githubusercontent.com/mlbendall/telescope/26a6c38df6d143387799edc2f3e9f976505f093d/telescope/_version.py. [cited 2026 Nov 15]. Available from: https://raw.githubusercontent.com/mlbendall/telescope/26a6c38df6d143387799edc2f3e9f976505f093d/telescope/_version.py
7. Geneimprint : Genes. [cited 2026 Nov 15]. Available from: https://www.geneimprint.com/site/genes-by-species
8. R-HSA-556833. [cited 2026 Nov 15]. Available from: https://reactome.org/ContentService/data/query/R-HSA-556833
9. GO%3A0006865%2CGO%3A0008643%2CGO%3A0006629%2CGO%3A0008203%2CGO%3A0009755%2CGO%3A0002682%2CGO%3A0001890. [cited 2026 Nov 15]. Available from: https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A0006865%2CGO%3A0008643%2CGO%3A0006629%2CGO%3A0008203%2CGO%3A0009755%2CGO%3A0002682%2CGO%3A0001890
10. R-HSA-382551. [cited 2026 Nov 15]. Available from: https://reactome.org/ContentService/data/query/R-HSA-382551
