# Conditional SNP-level parental counts, fitted separately for each gene.
source(file.path(root, "R/h1_annotations.R"), local = TRUE)
source(file.path(root, "R/h1_prepare.R"), local = TRUE)
source(file.path(root, "R/h1_backend.R"), local = TRUE)
# Capture the source root now; preparation does not require the tool environment.
h1_cmdstan_runtime <- local({
  project_root <- root
  function(...) {
    runtime <- new.env(parent = baseenv())
    sys.source(file.path(project_root, "R/cmdstan_runtime.R"), envir = runtime)
    runtime$poe_cmdstan_runtime(project_root, ...)
  }
})
H1_MODEL_SPEC <- "h1_joint_conditional_v1"
h1_settings <- function(cfg) {
  x <- h1_normalize_annotations(nb_default(cfg$h1_bayes, list()))
  out <- list(chains = 4L, warmup = 1000L, iterations = 1000L,
    adapt_delta = .99, max_treedepth = 12L, seed = 1147L, parallel_chains = 4L,
    count_provenance = NA_character_, producing_pipeline = NA_character_)
  for (k in intersect(names(x), names(out))) out[[k]] <- x[[k]]
  h1_validate_controls(out)
  out
}
h1_profiles <- function() data.table(profile = c("primary", "half_effect_scales",
  "double_effect_scales", "precision_5", "precision_80"), effect_scale = c(1, .5, 2, 1, 1),
  precision_center = c(20, 20, 20, 5, 80))
h1_formula <- function() brms::bf(
  P | trials(N) ~ 0 + Intercept + cross_c + offset(dosage_offset) +
    (1 + cross_c || transcript) + (1 | cluster) + (1 | tx_cluster) + (1 | marker),
  phi ~ 1 + (1 | transcript))
h1_priors <- function(profile) {
  f <- profile$effect_scale
  c(brms::set_prior(sprintf("normal(0,%g)", f), class = "b", coef = "Intercept"),
    brms::set_prior(sprintf("normal(0,%g)", f), class = "b", coef = "cross_c"),
    brms::set_prior(sprintf("normal(0,%g)", f), class = "sd", group = "transcript"),
    brms::set_prior(sprintf("normal(0,%g)", .5 * f), class = "sd", group = "cluster"),
    brms::set_prior(sprintf("normal(0,%g)", .5 * f), class = "sd", group = "tx_cluster"),
    brms::set_prior(sprintf("normal(0,%g)", .5 * f), class = "sd", group = "marker"),
    brms::set_prior(sprintf("normal(%g,1.5)", log(profile$precision_center)),
      class = "Intercept", dpar = "phi"),
    brms::set_prior("normal(0,0.5)", class = "sd", group = "transcript", dpar = "phi"))
}
h1_counts <- function(z, annotation, cfg, s) {
  if (!nrow(z)) z <- data.table(feature_id = character(), sample_id = character(), SNPpos = character(),
    maternal_count = numeric(), paternal_count = numeric(), cluster_id = character(),
    cross = character(), independence_verified = logical())
  require_columns(z, c("feature_id", "sample_id", "SNPpos", "maternal_count", "paternal_count",
    "cluster_id", "cross", "independence_verified"))
  assert_unique(z, c("feature_id", "sample_id", "SNPpos"))
  d <- copy(z)
  d[, known_parental := is.finite(maternal_count) & is.finite(paternal_count)]
  if (any(d$maternal_count < 0, na.rm = TRUE) || any(d$paternal_count < 0, na.rm = TRUE))
    stop("Negative allele counts")
  audit <- d[, .(n_input_rows = .N, n_missing_parental_rows = sum(!known_parental),
    n_fractional_parental_values = sum(maternal_count != round(maternal_count), na.rm = TRUE) +
      sum(paternal_count != round(paternal_count), na.rm = TRUE),
    absolute_rounding_change = sum(abs(maternal_count - round(maternal_count)), na.rm = TRUE) +
      sum(abs(paternal_count - round(paternal_count)), na.rm = TRUE)), by = feature_id]
  d <- d[known_parental == TRUE]
  d[, `:=`(P = round(paternal_count), M = round(maternal_count))]
  d[, N := M + P]
  zero <- d[, .(n_zero_total_rows = sum(N == 0)), by = feature_id]
  d <- d[N > 0]
  if (any(!d$cross %in% c("focal_mother", "alternate_mother")))
    stop("Bayesian H1 requires established reciprocal direction metadata")
  meta <- unique(d[, .(sample_id, cluster_id, cross, independence_verified)])
  if (anyDuplicated(meta$sample_id) ||
      any(meta[, uniqueN(cross), by = cluster_id]$V1 > 1)) stop("Inconsistent biological clusters")
  stats <- d[, .(n_snps = uniqueN(SNPpos), n_samples = uniqueN(sample_id), total_counts = sum(N),
    n_cross0_clusters = uniqueN(cluster_id[cross == "focal_mother" & independence_verified %in% TRUE]),
    n_cross1_clusters = uniqueN(cluster_id[cross == "alternate_mother" & independence_verified %in% TRUE]),
    all_independent = all(independence_verified %in% TRUE & !is.na(cluster_id) & nzchar(cluster_id))),
    by = feature_id]
  cov <- merge(annotation[, .(feature_id = transcript_id, gene_id, promoter_te_label)],
    audit, by = "feature_id", all.x = TRUE)
  cov <- merge(cov, zero, by = "feature_id", all.x = TRUE)
  cov <- merge(cov, stats, by = "feature_id", all.x = TRUE)
  cov[, eligibility := fcase(!promoter_te_label %in% c("TE", "nonTE"), "unknown_TE_class",
    is.na(n_snps), "no_informative_counts", n_snps < s$min_snps, "insufficient_snps",
    n_samples < s$min_samples, "insufficient_samples", total_counts < s$min_total_counts,
    "insufficient_counts", !all_independent %in% TRUE, "independence_unverified",
    n_cross0_clusters < s$min_per_cross | n_cross1_clusters < s$min_per_cross,
    "insufficient_reciprocal_clusters", default = "eligible")]
  d <- d[feature_id %in% cov[eligibility == "eligible"]$feature_id]
  d[, cross_c := fifelse(cross == "focal_mother", -.5, .5)]
  d[, dosage_offset := log(cfg$dosage$paternal / cfg$dosage$maternal)]
  list(data = d, coverage = cov)
}
h1_prepare <- function(cfg, s, root) {
  if (cfg$profile != "public_cross") stop("Bayesian H1 uses reciprocal small-cohort profiles")
  hs <- h1_settings(cfg); nb_te_require(cfg, s, root); nb_check_prepare(cfg, s, root)
  outputs <- nb_path(s, c("h1_transcript_coverage.tsv", "h1_gene_coverage.tsv", "h1_bayes_jobs.tsv"))
  if (nb_current(s, "h1-inputs", NULL)) return(invisible(TRUE))
  a <- reg_read(nb_path(s, "transcript_TE_annotation.tsv"))
  coverage <- h1_prepare_coverage(a, cfg, s)
  genes <- coverage[, .(n_TE = sum(eligibility == "eligible" & promoter_te_label == "TE"),
    n_nonTE = sum(eligibility == "eligible" & promoter_te_label == "nonTE"),
    n_transcripts = .N, n_count_eligible = sum(eligibility == "eligible")), by = gene_id]
  genes[, eligibility := fifelse(n_TE > 0 & n_nonTE > 0, "eligible", "no_eligible_mixed_classes")]
  jobs <- copy(genes[eligibility == "eligible"]); setorder(jobs, gene_id); jobs[, job_index := .I]
  coverage[, te_tss_window_bp := reg_settings(cfg)$te_tss_window_bp]
  genes[, `:=`(count_provenance = hs$count_provenance,
    assignment_covariance = "shared_random_effects_leave_unknown_read_assignment_covariance_unestimated")]
  for (i in seq_along(outputs)) nb_write(list(coverage, genes, jobs)[[i]], outputs[i])
  nb_finish(s, "h1-inputs", list(settings = hs, profiles = h1_profiles(),
    parent = nb_require(s, "te-inputs")$key), outputs)
  message(sprintf("H1 inputs complete: transcripts=%d; genes=%d; eligible mixed-class jobs=%d",
    nrow(coverage), nrow(genes), nrow(jobs)))
}
h1_encode <- function(d) {
  d <- copy(d)
  ids <- sort(unique(d$feature_id))
  map <- data.table(feature_id = ids, transcript = sprintf("tx%06d", seq_along(ids)))
  d <- merge(d, map, by = "feature_id")
  d[, cluster := factor(sprintf("cl%06d", match(cluster_id, sort(unique(cluster_id)))))]
  d[, marker := factor(sprintf("mk%06d", match(SNPpos, sort(unique(SNPpos)))))]
  d[, transcript := factor(transcript, levels = map$transcript)]
  d[, tx_cluster := interaction(transcript, cluster, drop = TRUE)]
  list(data = d, map = map)
}
h1_derived <- function(draws, map, annotation) {
  # Arrays retain chain identity for Monte Carlo diagnostics of every estimand.
  out <- list()
  beta <- kappa <- list()
  for (i in seq_len(nrow(map))) {
    nm <- map$transcript[i]
    beta[[i]] <- draws[, , "b_Intercept"] + draws[, , paste0("r_transcript[", nm, ",Intercept]")]
    kappa[[i]] <- draws[, , "b_cross_c"] + draws[, , paste0("r_transcript[", nm, ",cross_c]")]
    id <- map$feature_id[i]
    out[[paste0("beta:", id)]] <- beta[[i]]
    out[[paste0("kappa:", id)]] <- kappa[[i]]
    out[[paste0("cross0:", id)]] <- beta[[i]] - .5 * kappa[[i]]
    out[[paste0("cross1:", id)]] <- beta[[i]] + .5 * kappa[[i]]
    out[[paste0("cancellation:", id)]] <- .5 *
      (abs(beta[[i]] - .5 * kappa[[i]]) + abs(beta[[i]] + .5 * kappa[[i]])) - abs(beta[[i]])
  }
  te <- annotation$promoter_te_label[match(map$feature_id, annotation$transcript_id)] == "TE"
  avg <- function(a, use, absolute = FALSE) Reduce("+", lapply(a[use],
    function(v) if (absolute) abs(v) else v)) / sum(use)
  for (class in c("TE", "nonTE")) {
    use <- if (class == "TE") te else !te
    out[[paste0(class, "_signed_mean")]] <- avg(beta, use)
    out[[paste0(class, "_magnitude_mean")]] <- avg(beta, use, TRUE)
  }
  out$Delta_mag <- out$TE_magnitude_mean - out$nonTE_magnitude_mean
  out$Delta_signed <- out$TE_signed_mean - out$nonTE_signed_mean
  ans <- array(NA_real_, c(dim(draws)[1:2], length(out)),
    dimnames = list(NULL, dimnames(draws)[[2]], names(out)))
  for (i in seq_along(out)) ans[, , i] <- out[[i]]
  posterior::as_draws_array(ans)
}
h1_summarize <- function(draws) {
  d <- as.data.table(posterior::summarise_draws(draws, mean, median,
    sd = stats::sd, low = function(x) unname(quantile(x, .025)),
    high = function(x) unname(quantile(x, .975)), rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail,
    mcse_mean = posterior::mcse_mean, mcse_sd = posterior::mcse_sd,
    mcse_median = function(x) unname(posterior::mcse_quantile(x, probs = .5)),
    mcse_low = function(x) unname(posterior::mcse_quantile(x, probs = .025)),
    mcse_high = function(x) unname(posterior::mcse_quantile(x, probs = .975))))
  d[, `:=`(prob_gt_zero = vapply(variable, function(v) mean(draws[, , v] > 0), 0),
    prob_lt_zero = vapply(variable, function(v) mean(draws[, , v] < 0), 0),
    mcse_prob_gt_zero = vapply(variable, function(v)
      posterior::mcse_mean((draws[, , v] > 0) * 1), 0),
    mcse_prob_lt_zero = vapply(variable, function(v)
      posterior::mcse_mean((draws[, , v] < 0) * 1), 0))]
  d
}
h1_diagnostics <- function(fit, derived, hs) {
  pars <- as.data.table(posterior::summarise_draws(posterior::as_draws_array(fit),
    "rhat", "ess_bulk", "ess_tail"))
  der <- h1_summarize(derived)
  all <- rbindlist(list(pars, der[, .(variable, rhat, ess_bulk, ess_tail)]))
  np <- as.data.table(brms::nuts_params(fit))
  energy <- np[Parameter == "energy__", .(ebfmi = mean(diff(Value)^2) / var(Value)), by = Chain]
  bad <- all[!is.finite(rhat) | rhat >= 1.01 | !is.finite(ess_bulk) | ess_bulk < 400 |
    !is.finite(ess_tail) | ess_tail < 400]
  divergent <- sum(np[Parameter == "divergent__"]$Value)
  tree <- sum(np[Parameter == "treedepth__"]$Value >= hs$max_treedepth)
  ready <- !nrow(bad) && divergent == 0 && tree == 0 &&
    nrow(energy) == hs$chains && all(is.finite(energy$ebfmi) & energy$ebfmi >= .3)
  list(status = if (ready) "ready" else "diagnostic_failure",
    summary = data.table(status = if (ready) "ready" else "diagnostic_failure",
      n_bad_parameters = nrow(bad), max_rhat = max(all$rhat, na.rm = TRUE),
      min_bulk_ess = min(all$ess_bulk, na.rm = TRUE), min_tail_ess = min(all$ess_tail, na.rm = TRUE),
      divergences = divergent, treedepth_hits = tree, min_ebfmi = min(energy$ebfmi),
      operator_action = if (ready) "none" else "inspect_fit_and_explicitly_configure_longer_run_in_fresh_output"),
    parameters = all, energy = energy)
}
h1_residuals <- function(fit, d) {
  ids <- unique(round(seq(1, posterior::ndraws(posterior::as_draws_array(fit)), length.out = 200)))
  mu <- brms::posterior_epred(fit, draw_ids = ids)
  phi <- brms::posterior_linpred(fit, dpar = "phi", transform = TRUE, draw_ids = ids)
  p <- sweep(mu, 2, d$N, "/")
  variance <- sweep(p * (1 - p), 2, d$N, "*") *
    (sweep(phi, 2, d$N, "+") / (1 + phi))
  residual <- sweep(-mu, 2, d$P, "+")
  rbindlist(lapply(split(seq_len(nrow(d)), d$feature_id), function(ii) data.table(
    feature_id = d$feature_id[ii[1]], n_rows = length(ii),
    mean_observed_fraction = mean(d$P[ii] / d$N[ii]),
    mean_fitted_fraction = mean(p[, ii, drop = FALSE]),
    mean_Pearson_squared = mean((residual^2 / pmax(variance, 1e-12))[, ii, drop = FALSE]),
    mean_absolute_count_residual = mean(abs(residual[, ii, drop = FALSE])))))
}
h1_fit_once <- function(d, ann, hs, profile, seed, file) {
  encoded <- h1_encode(d)
  poe_h1_backend_check()
  if (file.exists(file)) {
    saved <- readRDS(file)
    if (!identical(saved$model_spec, H1_MODEL_SPEC) ||
        !identical(saved$settings, hs) || !h1_task_table_equal(saved$profile, profile) ||
        !identical(saved$seed, seed) || !h1_task_table_equal(saved$map, encoded$map))
      stop("H1 task settings differ; use a fresh output directory")
    summary <- h1_summarize(saved$derived); summary[, status := saved$diagnostics$status]
    return(list(summary = summary, diagnostics = saved$diagnostics$summary,
      residuals = h1_residuals(saved$fit, encoded$data)))
  }
  h1_cmdstan_runtime()
  fit <- brms::brm(h1_formula(), data = encoded$data, family = brms::beta_binomial(),
    prior = h1_priors(profile), backend = "cmdstanr", chains = hs$chains,
    iter = hs$warmup + hs$iterations, warmup = hs$warmup, cores = hs$parallel_chains,
    seed = seed, threads = brms::threading(1), refresh = 0,
    control = list(adapt_delta = hs$adapt_delta, max_treedepth = hs$max_treedepth),
    save_pars = brms::save_pars(all = TRUE))
  draws <- h1_derived(posterior::as_draws_array(fit), encoded$map, ann)
  diagnostics <- h1_diagnostics(fit, draws, hs)
  summary <- h1_summarize(draws); summary[, status := diagnostics$status]
  nb_save(list(model_spec = H1_MODEL_SPEC, fit = fit, map = encoded$map, derived = draws, settings = hs,
    profile = profile, seed = seed, diagnostics = diagnostics), file)
  list(summary = summary, diagnostics = diagnostics$summary, residuals = h1_residuals(fit, encoded$data))
}
h1_gene <- function(cfg, s, root, job_index, retry_failed = FALSE) {
  hs <- h1_settings(cfg); nb_require(s, "h1-inputs")
  id <- as.integer(job_index)
  job <- reg_read(nb_path(s, "h1_bayes_jobs.tsv"))[job_index == id]
  if (nrow(job) != 1) stop("Unknown H1 job index")
  ann <- reg_read(nb_path(s, "transcript_TE_annotation.tsv"))[gene_id == job$gene_id]
  counts <- h1_counts(nb_gene_data(s, ann$transcript_id), ann, cfg, s)
  d <- counts$data; ann <- ann[transcript_id %in% d$feature_id]
  stage <- sprintf("h1_%06d", job_index)
  paths <- nb_path(s, "h1_bayes", sprintf("gene_%06d_%s.tsv", job_index,
    c("posterior", "diagnostics", "residuals")))
  if (nb_current(s, stage, NULL)) {
    prior <- reg_read(paths[2L])
    require_columns(prior, "status")
    if (!any(prior$status %in% "fit_failure")) {
      message("H1 gene already completed; retaining saved fits and diagnostic status")
      return(invisible(TRUE))
    }
    if (!isTRUE(retry_failed)) stop("H1 gene has recorded fit failures. After backend validation use ",
      "'h1-retry-failed ", job_index, "'; completed diagnostic failures are not automatically refitted.")
  }
  poe_h1_backend_check()
  h1_cmdstan_runtime()
  if (isTRUE(retry_failed) && any(file.exists(c(paths, nb_receipt(s, stage)))))
    h1_archive_gene_failure(s, job_index, paths)
  profiles <- h1_profiles()
  tasks <- rbindlist(list(profiles[, .(profile, omitted_cluster = "")],
    data.table(profile = "primary", omitted_cluster = sort(unique(d$cluster_id)))))
  summaries <- diagnostics <- residuals <- list()
  for (i in seq_len(nrow(tasks))) {
    message(sprintf("H1 gene %d task %d/%d: profile=%s; %s", job_index, i, nrow(tasks),
      tasks$profile[i], if (nzchar(tasks$omitted_cluster[i])) "cluster omission" else "full data"))
    profile <- profiles[profile == tasks$profile[i]]
    omit <- tasks$omitted_cluster[i]
    dd <- if (nzchar(omit)) d[cluster_id != omit] else d
    support <- if (nrow(dd)) h1_counts(dd, ann, cfg, s) else NULL
    eligible <- !is.null(support) && all(support$coverage$eligibility == "eligible")
    seed <- as.integer((hs$seed + job_index * 1009 + i * 100003) %% .Machine$integer.max)
    result <- if (!eligible) list(diagnostics = data.table(status = "unavailable_LOCO_support",
      operator_action = "retain_eligibility_thresholds")) else tryCatch(
      h1_fit_once(dd, ann, hs, profile, seed,
        nb_path(s, "private", "h1_bayes", sprintf("gene_%06d_task_%04d.rds", job_index, i))),
      error = function(e) {
        message("H1 task ", i, " failed: ", conditionMessage(e))
        list(diagnostics = data.table(status = "fit_failure",
          error = conditionMessage(e), operator_action = "inspect_error_and_rerun_explicitly"))
      })
    label <- function(x) {if (is.null(x)) return(NULL)
      x[, `:=`(gene_id = job$gene_id, profile = tasks$profile[i], task_index = i,
        omitted_cluster = if (nzchar(omit)) paste0("omission_", match(omit, sort(unique(d$cluster_id)))) else "",
        seed = seed, te_tss_window_bp = reg_settings(cfg)$te_tss_window_bp)]; x}
    summaries[[i]] <- label(result$summary)
    diagnostics[[i]] <- label(result$diagnostics)
    residuals[[i]] <- label(result$residuals)
  }
  empty <- data.table(gene_id = character(), status = character())
  for (i in seq_along(paths)) {
    x <- rbindlist(list(summaries, diagnostics, residuals)[[i]], fill = TRUE)
    nb_write(if (ncol(x)) x else empty, paths[i])
  }
  nb_finish(s, stage, list(settings = hs, profiles = profiles, gene_id = job$gene_id,
    parent = nb_require(s, "h1-inputs")$key), paths)
  message("H1 gene tasks recorded. Review diagnostics; completion does not imply numerical readiness.")
}
h1_combine <- function(cfg, s, root) {
  nb_require(s, "h1-inputs")
  jobs <- reg_read(nb_path(s, "h1_bayes_jobs.tsv"))
  coverage <- reg_read(nb_path(s, "h1_gene_coverage.tsv"))
  for (j in jobs$job_index) nb_require(s, sprintf("h1_%06d", j))
  # Full-data rows are labelled omitted_cluster = "". A gene table with no
  # cluster-omission rows (every omission lacked support) has an all-blank
  # column, which the reader returns as logical NA; restore the label so the
  # full-data filters below keep those genes.
  read_gene <- function(path) {
    x <- reg_read(path)
    if ("omitted_cluster" %in% names(x)) {
      x[, omitted_cluster := as.character(omitted_cluster)]
      x[is.na(omitted_cluster), omitted_cluster := ""]
    }
    x
  }
  objects <- lapply(c("posterior", "diagnostics", "residuals"), function(type) {
    if (!nrow(jobs)) return(data.table(gene_id = character(), status = character()))
    rbindlist(lapply(nb_path(s, "h1_bayes",
      sprintf("gene_%06d_%s.tsv", jobs$job_index, type)), read_gene), fill = TRUE)
  })
  diag <- objects[[2]]
  coverage[, fit_status := "ineligible"]
  if (nrow(diag)) {
    primary <- diag[profile == "primary" & omitted_cluster == ""]
    coverage[primary, on = "gene_id", fit_status := i.status]
  }
  coverage[eligibility == "eligible" & fit_status == "ineligible", fit_status := "missing_job"]
  objects[[4]] <- coverage
  p <- objects[[1]]
  objects[[5]] <- if ("variable" %in% names(p)) p[
    variable %in% c("Delta_mag", "Delta_signed") & omitted_cluster == "",
    .(n_genes = .N, n_ready_genes = sum(status == "ready"),
      descriptive_equal_gene_mean = mean(mean), descriptive_equal_gene_median = median(median),
      interpretation = "separate_gene_fits_cross_gene_dependence_unestimated"),
    by = .(profile, variable)] else data.table(status = "no_posterior_available")
  files <- nb_path(s, paste0("h1_", c("posterior.tsv.gz", "diagnostics.tsv", "residuals.tsv",
    "coverage.tsv", "cohort_descriptive.tsv")))
  for (i in seq_along(files)) nb_write(objects[[i]], files[i])
  nb_finish(s, "te-combine", list(method = "joint_gene_beta_binomial", settings = h1_settings(cfg)), files)
}
