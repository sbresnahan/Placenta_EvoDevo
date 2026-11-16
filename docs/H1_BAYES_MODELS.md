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

Count-provenance and producing-pipeline labels are optional descriptive
annotations. Reuse unchanged counts without a new version/provenance audit.
Missing labels and shipped annotation placeholders remain unrecorded; they do
not imply independent fragment assignment or duplicate correction. Valid
existing annotations and cohort-specific settings are retained.
Known duplicate-fragment assignment requires upstream correction, and explicit
invalid assignment labels remain rejected. With unknown assignment covariance,
the model's shared effects leave read-level allocation dependence unestimated.
MCMC readiness describes computational reliability conditional on these assumptions.

## Preparation resources

`h1-inputs` reads the two transcript batch indexes once and each required
count batch once. Features remain complete within their indexed batch.
Workers apply the unchanged count gates to per-gene portions and return
coverage plus compact metadata. The coordinator checks biological clusters
across batches before writing any of the three output tables. This includes
count-ineligible transcripts with informative rows. Completed stages are reused.

`H1_PREP_WORKERS` defaults to 1. Set `H1_PREP_MEMORY_GB` to the total memory
budget when requesting multiple workers. Unix process parallelism is capped
by the requested workers, required batch count, scheduler CPU controls when
present, and the memory estimate. The clean container forwards these controls
and limits native-library threads to one per preparation process.
The coordinator allowance is 1 GiB plus four annotation-table sizes, two
feature-index sizes and one batch-index size. Each worker allowance is four
times the largest required batch's recorded `object_bytes` plus 256 MiB.
The controls use 1024-based GB. This is a conservative planning estimate,
not an enforced RSS ceiling. Missing size metadata select batch-wise serial
loading with an explicit message; no memory assurance is possible in that case.
With size metadata but no budget, a multi-worker request also selects serial
loading. A supplied budget below the estimated one-worker requirement fails
before count loading.

The submission helper derives preparation controls from `NB_CPUS` and
`NB_MEMORY_GB`, rejects controls above those requests, requests one host and
rejects preparation arrays. Direct interactive calls must set a budget and
worker count within the actual allocation. Progress reports show effective
workers, batch start/completion, elapsed times and the final eligible job count.
Do not start concurrent preparations in one output directory. An interrupted
incomplete preparation restarts; it has no per-gene checkpoints.

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

Fitting first loads the RStan/brms/CmdStanR/posterior backend and validates empty
RStan class construction, before new sampling and outside per-task error
handling. This preflight is read-only. The existing R 4.3.1 shared library can be
repaired with `bin/repair_rstan_stack.R` using staged, same-version builds and
explicit apply/rollback maintenance confirmation; see runbook section 0A.
Fitting itself never installs or rebuilds packages.

Normal `h1-fit JOB_INDEX` resumes interrupted genes and reuses compatible saved
tasks. A completed gene with `fit_failure` tasks stops with instructions for
`h1-retry-failed JOB_INDEX`. That explicit retry preserves prior gene-level
failure files/receipt and reuses compatible task objects. Runtime failures,
numerical diagnostic failures and unsupported omissions remain distinct.
Saved profile/mapping comparisons retain column names, types, values and row
order while ignoring process-local table metadata. Model specification,
settings and seeds remain strict. Completed diagnostic failures are not
automatically refitted, and diagnostic thresholds are unchanged.

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
