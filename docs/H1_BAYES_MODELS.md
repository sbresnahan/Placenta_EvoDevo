# H1 conditional parental-count model

## Scope and estimands

Fit each reciprocal-cohort gene separately. A transcript's beta is its
dosage-adjusted natural-log paternal/maternal effect, balanced equally across
the two established reciprocal directions.

Primary: `Delta_mag = mean_TE(abs(beta)) - mean_nonTE(abs(beta))`.
Secondary: `Delta_signed = mean_TE(beta) - mean_nonTE(beta)`.
Report both TE and nonTE class means. A positive signed difference can reflect
weaker maternal bias. Cross-specific effects are `beta - kappa/2` and
`beta + kappa/2`. Cancellation is the mean cross-specific magnitude minus the
magnitude of balanced beta.

The adjusted between-gene analysis uses `beta^2 - SE^2` and retains its
propensity, block-bootstrap and permutation design.

## Likelihood

For each original transcript–SNP–sample row:

```text
P ~ beta_binomial(N, p, phi_transcript)
logit(p) = log(dosage_p/dosage_m) + beta_transcript +
  cross_c * kappa_transcript + u_cluster + v_transcript_cluster + w_marker
cross_c = -0.5 (focal_mother), +0.5 (alternate_mother)
```

Hierarchical beta and kappa priors share information across transcripts within
the gene. TE labels enter posterior contrasts. The log-precision has a gene
intercept and transcript random effects. Cluster, transcript-cluster and full
genomic-marker effects account for the specified dependence structure.
The population intercept uses brms `0 + Intercept` coding so its prior applies
directly to balanced beta even with unequal reciprocal sample sizes.

Conditioning on total parental counts cancels allele-independent abundance,
library and covariate factors. The upstream NB POE and H6/H7 likelihoods
remain separate.

## Support and count handling

Read immutable original transcript-SNP batches through their preparation
manifests. Use ≥3 informative SNPs, ≥50 rounded total parental counts,
≥4 observed samples and ≥2 verified independent clusters in each direction,
for every included transcript. Require both known promoter TE classes within
the gene. Count eligibility is independent of upstream NB convergence.

Keep missing parental assignments missing. Round maternal and paternal counts
separately at fitting, audit fractional changes, and omit zero-total rows.
Retain each transcript's original markers. Transcript marker sets may be disjoint.

Declare count provenance and producing pipeline. Known duplicate-fragment
assignment requires upstream correction. With unknown assignment covariance,
the model's shared effects leave read-level allocation dependence unestimated.
MCMC readiness describes computational reliability conditional on these assumptions.

## Priors

| Parameter | Primary prior |
|---|---|
| Gene beta/kappa mean | Normal(0,1) |
| Transcript beta/kappa SD | half-Normal(0,1) |
| Cluster/transcript-cluster/marker SD | half-Normal(0,0.5) |
| Log precision intercept | Normal(log(20),1.5) |
| Transcript log precision SD | half-Normal(0,0.5) |

Five fixed profiles: primary; half effect/random-effect scales; double scales;
primary effects with precision centered at 5; primary effects with precision
centered at 80. Precision widths are fixed. Leave-one-cluster-out fits use the
primary profile. An omission violating any transcript support gate is unavailable.

## Sampling and diagnostics

R implementation: brms with cmdstanr and CmdStan. Default: four chains,
1000 warmup and 1000 retained iterations per chain, one thread per chain,
adapt_delta 0.95 and maximum tree depth 12. Task-specific seeds are recorded.
Size concurrent chains to the CPU allocation and pilot measured peak memory.

Ready fits require Rhat <1.01, bulk and tail ESS ≥400 for model parameters and
derived quantities, zero divergences, zero tree-depth hits and E-BFMI ≥0.3 in
each chain. Failed diagnostics remain in coverage. Longer runs require an
explicit configuration and fresh output directory. Completed private task
objects support interrupted-job continuation with identical settings.

Observed-data residual summaries include fitted parental fractions, absolute
count residuals and Pearson residual dispersion. These checks complement the
sampler diagnostics.

## Outputs and interpretation

`h1_transcript_coverage.tsv` records every reference transcript and count gate.
`h1_gene_coverage.tsv` records mixed-class eligibility. `h1_bayes_jobs.tsv`
assigns stable sequential jobs.

`h1_posterior.tsv.gz` contains posterior means, medians, SDs, 95% equal-tail
intervals, directional probabilities, Monte Carlo errors, Rhat and ESS.
`h1_diagnostics.tsv` retains every prior/omission task and its status.
`h1_residuals.tsv` records observed-data diagnostics.
`h1_coverage.tsv` retains all genes and primary fit status.
`h1_cohort_descriptive.tsv` gives equal-gene descriptive averages. Gene fits
share specimens, so cross-gene posterior dependence remains unestimated.
The cohort table has no posterior combined uncertainty.

No posterior discovery cutoff or frequentist p/q/FDR labels are assigned.
H5/H6/H7 entry does not depend on an H1 significance or readiness gate.

## GUSTO

Use completed transcript POE models for descriptive class contrasts and
recompute count-refit bootstrap diagnostics for the active radius.
`*_bootstrap_centered_tail_prob` names the centered-tail diagnostic;
intervals and SEs require ≥90% successful refits. No H1 BH or primary
frequentist-discovery columns are assigned. A replacement GUSTO inferential
model is outside this implementation.
