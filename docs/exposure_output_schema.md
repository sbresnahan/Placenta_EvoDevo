# Exposure module output schema

All outputs initially belong to a private run. TSVs use headers and NA for
unavailable quantities. One planned job is a module/cohort/level/feature/family/
exposure/predictor combination. A failed job is retained, not converted to a
zero effect. Source-path/error columns are removed from shareable exports.

| File | Unit | Principal fields |
|---|---|---|
| preflight.tsv | cohort | crosswalk/quantified/unquantified n, loci/families/pair n, private input bytes, H7 readiness/blocker |
| exposure_coverage.tsv | exposure | n_crosswalk, n_available, n_unavailable, before feature-specific exclusions |
| TE_scaling.tsv | gene/family/predictor | center, scale, n_expressed, predictor_status |
| job_inventory.tsv | module/level/predictor | jobs, distinct batches |
| private/jobs.tsv | planned job | job_id, batch, source batch/input_file, feature/gene/family/exposure/predictor |
| private/TE_predictors.tsv.gz | sample/gene/family/predictor | family_count, cpm, log_cpm, TE_abundance, scaling/status |
| private/prepared.rds | cohort | validated metadata, annotation, manifest, targets, predictor cache, jobs |
| private/run_identity.rds | run | explicit settings/dosage/input size+mtime, readable_jobs_v2 identity schema |
| pilot_results.tsv | pilot job | model fields, elapsed_seconds, input_object_bytes |
| pilot_resources.tsv | module/level | n_planned, n_pilot, n_estimable, median_fit_seconds, max_input_bytes, projected_serial_seconds |
| private/fits/* | batch/stage | every model job row; separate *_timing.tsv contains worker timing |
| model_results.tsv | planned job | model, prediction, mediation and correction fields below |
| model_coverage.tsv | module/level/predictor/status | N, or explicit empty inventory |
| bootstrap_selection.tsv | selected mediation pair | selection_index, job_id, level/feature/family/exposure/predictor |
| bootstrap/*.tsv | selected pair | B, n_valid, valid_fraction, bootstrap_status, *_boot_low/high |
| private/bootstrap_draws/* | selected pair/draw | draw, status, refitted effect metrics |
| frozen_annotation_manifest.tsv | set/species | category/name/source/source_version/retrieved_date |
| private/frozen_annotations.rds | run | frozen memberships, manifest, freeze time |
| gene_evidence.tsv | gene/group | gene_p, n_tested, n_estimable, median counts/units/SNPs, gene length |
| functional_enrichment.tsv | set/group | n_set/background, rank difference, p/q, permutations, retained matching variables, status |
| H5_support_enrichment.tsv | modality/group | assessable n, supported/discovery n, continuous evidence, logistic OR/CI/p/q or sparse/failure status |
| provenance.tsv | run | module/base version, cohort/species/build, covariates, explicit causal assumptions |
| report_summary.tsv | module/level/predictor | planned, estimable, discoveries |
| aggregate_export/MANIFEST.tsv | exported file | relative filename, bytes only |

`job_id` is a sequential `job_00000001`-style label after stable inventory sorting;
it is unique only within its prepared run/cohort. Do not combine job IDs across runs.
TE quantification QC declares `annotation_id` (catalogue ID/release), not verified
content identity. File sizes and times are descriptive metadata, not integrity proof.

Independent reporting uses `bin/export_results.R --exposure-nb --out FILE.tar.gz CONFIG...`.
It requires the aggregate manifest plus model_results, model_coverage, job_inventory,
and provenance tables, not private inputs or an initial-report archive.
`family_id` is retained because it denotes a TE family, not a participant.
Each cohort is under `cohorts/<id>/exposure_nb_v1/`. The archive includes filenames/
sizes, optional-output availability, docs, and descriptive code inventory.

## Model fields

- Identity: job_id, batch, module (`h6`, `h7_total`, `h7_mediation`,
  `h7_modification`), cohort_id, species, level, feature_id, gene_id,
  family_id, exposure, predictor, version, status.
- Support: n_samples, n_biological_units, n_snps, total_counts, mode,
  small_cohort, independence_verified, sample_design_rank, residual_unit_df,
  condition_number, E_min/max, E_p10/p50/p90, min_joint_quadrant when relevant.
- Fit: formula, optimizer_code, pdHess, warnings (private only).
- Interaction: slope/se/p, b/se/p, modification/se/p as applicable.
  `slope` is allele×E; in the mediation model this is the conditional direct
  PFAS path, not the indirect effect.
- Diagnostics: pfas_TE_correlation, pearson_RMS, pearson_abs_max,
  largest_family_residual_fraction, largest_family_count_fraction;
  mediator_residual_RMS, mediator_max_leverage, mediator_max_Cook_distance.
- Predictions: beta_p10/low/high, beta_p50, beta_p90/low/high,
  delta_signed/low/high, delta_magnitude/low/high,
  direction_reversal_fitted/supported, prediction_M, prediction_scale.
  Low/high here mean confidence limits, not expression direction.
- Paths: a/se/p/df, b/se/p, c_prime/se/p, indirect_beta/low/high,
  indirect_path_p, implied_total_beta/low/high, interval_method,
  confounding_assumption, temporal_ordering.
  indirect_beta and implied_total_beta are NA in secondary modification.
- Standardization: scenario_E_low/high, beta_scenario_00/10/11,
  natural_direct/indirect/total_beta,
  natural_direct/indirect/total_magnitude,
  fraction_direction_reversal, mediator_scenario_outside_observed.
- Correction: slope_q, indirect_q, direct_q, modification_q and
  *_q_within_exposure; slope_q_pooled_cohorts is NA until multiple cohorts are
  explicitly combined; primary_p/q and discovery select the module's relevant
  test. FDR discovery concerns signed-slope/indirect-path/modification evidence,
  not every descriptive prediction or magnitude interval.

Column sets are unioned across modules. A quantity is absent/NA when it cannot
be estimated or is not meaningful. Statuses describe insufficient units/counts,
missing/unidentifiable covariates, no within-cross variation, rank/conditioning
failure, nonconvergence/Hessian, nonfinite uncertainty, mediator failure or
input/model exception. Do not drop these rows before coverage/multiplicity checks.

The `quantification_plan` folder contains generated GTF, annotation QC, per-library
shell scripts and a jobs table with `executed=FALSE`. `BAM_probe.tsv` is a bounded
tag probe, not proof that the entire library or tool environment is valid.
