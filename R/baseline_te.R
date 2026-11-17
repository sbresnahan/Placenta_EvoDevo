# Comparisons of individually fitted transcript/gene POE. No shared-marker join.
nb_estimable <- function(d) d[status %in% c("eligible_NB", "exploratory_NB") &
  is.finite(beta_adjusted) & is.finite(se) & se > 0]
nb_repeat_path <- function(cfg) {
  z <- reg_settings(cfg)
  if (cfg$species == "fish") reg_fish_paths(z)$intervals else
    file.path(z$output_dir, "repeat_intervals.tsv.gz")
}
nb_te_resource_key <- function(cfg, s, root) {
  z <- reg_settings(cfg); path <- nb_repeat_path(cfg)
  candidates <- c(path, file.path(z$output_dir, "resource_manifest.tsv"),
    if (!is.null(z$contig_aliases)) z$contig_aliases,
    if (cfg$species == "fish") vapply(c("regulatory-resources", "fish-discovery", "fish-masking"),
      function(stage) reg_receipt_path(z, stage), ""))
  state <- setNames(vapply(candidates, function(p)
    nb_file_key(p)[[1]], ""), candidates)
  list(nb_require(s, "combine")$key, resources = state,
    reference = nb_file_key(cfg$paths$tx_map)[[1]], code = nb_code_version(root))
}
nb_te_require <- function(cfg, s, root) {
  nb_require(s, "te-inputs", nb_te_resource_key(cfg, s, root))
}
nb_te_annotation <- function(tx, repeats = NULL, full_reference = FALSE, radius = 3000L, contigs = NULL) {
  x <- copy(tx)
  require_columns(x, c("transcript_id", "gene_id", "chrom", "tss"))
  x[, `:=`(tss_te_label = "not_assessed", te_class = NA_character_,
    te_family = NA_character_, te_ids = NA_character_)]
  if (is.null(repeats) || !nrow(repeats))
    return(reg_promoter_annotation(x, repeats, radius, contigs, full_reference))
  x[, tss_te_label := if (full_reference) "nonTE" else "not_assessed"]
  h <- reg_hits(x[, .(chrom, start = tss, end = tss)], repeats)
  if (nrow(h)) {
    a <- cbind(h, repeats[h$subject, .(is_te, is_unknown, repeat_class, repeat_family, repeat_id)])
    a <- a[, .(tss_te_label = if (any(is_unknown)) "unknown_repeat" else
      if (any(is_te)) "TE" else if (full_reference) "nonTE" else "not_assessed",
      te_class = paste(sort(unique(repeat_class[is_te])), collapse = ";"),
      te_family = paste(sort(unique(repeat_family[is_te])), collapse = ";"),
      te_ids = paste(sort(unique(repeat_id[is_te])), collapse = ";")), by = query]
    x[a$query, `:=`(tss_te_label = a$tss_te_label, te_class = a$te_class,
      te_family = a$te_family, te_ids = a$te_ids)]
  }
  reg_promoter_annotation(x, repeats, radius, contigs, full_reference)
}
nb_te_prepare <- function(cfg, s, root) {
  nb_check_prepare(cfg, s, root); nb_require(s, "combine")
  path <- nb_repeat_path(cfg); repeats <- NULL; full <- FALSE
  z <- reg_settings(cfg)
  if (file.exists(path)) {
    repeats <- reg_read(path)
    require_columns(repeats, c("chrom", "start", "end", "is_te", "is_unknown",
      "repeat_class", "repeat_family", "repeat_id"))
    manifest <- file.path(z$output_dir, "resource_manifest.tsv")
    full <- if (cfg$species == "fish") reg_valid(z, "fish-masking") else
      file.exists(manifest) && any(reg_read(manifest)[resource_id == "repeats"]$scope == "full_reference")
    target <- reg_read(cfg$paths$tx_map)$chrom
    aliases <- if (!is.null(z$contig_aliases)) reg_read(z$contig_aliases) else NULL
    repeats[, chrom := reg_alias(chrom, unique(target), aliases)]
    repeats <- repeats[!is.na(chrom)]
    if (!nrow(repeats)) stop("No repeat contigs match transcript reference")
  }
  key <- nb_te_resource_key(cfg, s, root)
  if (nb_current(s, "te-inputs", key)) return(invisible(TRUE))
  tx <- nb_te_annotation(reg_read(cfg$paths$tx_map), repeats, full,
    z$te_tss_window_bp, reg_native_contigs(cfg))
  nb_write(tx, nb_path(s, "transcript_TE_annotation.tsv"))
  fits <- reg_read(nb_path(s, "transcript_POE.tsv"))
  d <- merge(tx, fits, by.x = "transcript_id", by.y = "feature_id", all.x = TRUE)
  d[is.na(status), status := "no_informative_counts"]
  nb_write(d, nb_path(s, "transcript_TE_POE.tsv"))
  eligible <- nb_estimable(d)
  mixed <- eligible[promoter_te_label %in% c("TE", "nonTE"), .(
    n_TE = sum(promoter_te_label == "TE"), n_nonTE = sum(promoter_te_label == "nonTE")), by = gene_id]
  mixed <- mixed[n_TE > 0 & n_nonTE > 0]
  setorder(mixed, gene_id); mixed[, job_index := .I]
  nb_write(mixed, nb_path(s, "within_gene_jobs.tsv"))
  availability <- data.table(cohort_id = cfg$cohort_id,
    repeats = if (is.null(repeats)) "not_available" else if (full) "full_reference" else "partial_reference",
    n_transcripts = nrow(tx), n_estimable_transcripts = nrow(eligible),
    te_tss_window_bp = z$te_tss_window_bp,
    n_TE = sum(tx$promoter_te_label == "TE"), n_nonTE = sum(tx$promoter_te_label == "nonTE"),
    n_unknown_or_unassessed = sum(!tx$promoter_te_label %in% c("TE", "nonTE")), n_mixed_genes = nrow(mixed),
    status = if (is.null(repeats)) "TE_comparisons_unavailable_missing_repeats" else
      if (!full) "TE_comparisons_limited_partial_repeat_reference" else "ready")
  nb_write(availability, nb_path(s, "TE_availability.tsv"))
  if (cfg$profile == "human_trio") {
    coverage <- d[, .(transcript_id, gene_id, promoter_te_label, status,
      eligibility = fifelse(!promoter_te_label %in% c("TE", "nonTE"), "unknown_TE_class",
        fifelse(is.finite(beta_adjusted) & is.finite(se) & se > 0 &
          status %in% c("eligible_NB", "exploratory_NB"), "eligible", status)),
      te_tss_window_bp)]
    nb_write(coverage, nb_path(s, "h1_transcript_coverage.tsv"))
  }
  nb_finish(s, "te-inputs", key, nb_path(s, c("transcript_TE_annotation.tsv",
    "transcript_TE_POE.tsv", "within_gene_jobs.tsv", "TE_availability.tsv")))
}
# Within-gene class contrasts of individually fitted transcript POE estimates.
# Raw |b| is inflated by estimation noise (E|b| = SE*sqrt(2/pi) at zero effect),
# so classes that differ in precision differ in raw magnitude even without a
# true difference. b^2 - SE^2 is unbiased for beta^2 given an accurate SE, as in
# the between-gene analysis; its class-mean difference is the primary magnitude
# field (squared-ln units). The root-scale version is a descriptive companion;
# raw |b| is retained as a secondary field.
NB_TE_METRICS <- c("signed", "magnitude", "magnitude_sq_corrected", "magnitude_rms_corrected")
NB_TE_STATS_VERSION <- "within_gene_stats_v2"
# Versioned stage names: genes or bootstrap chunks recorded by an earlier
# statistic definition are recomputed rather than mixed with new ones.
nb_te_stage <- function(job_id, chunk = NULL)
  paste0(sprintf("te_%06d_%s", as.integer(job_id), NB_TE_STATS_VERSION),
    if (!is.null(chunk)) sprintf("_chunk_%05d", as.integer(chunk)) else "")
nb_te_rms <- function(q) sqrt(pmax(q, 0))
nb_te_statistics <- function(estimates, annotations) {
  annotations <- copy(annotations)
  if (!"promoter_te_label" %in% names(annotations))
    annotations[, promoter_te_label := tss_te_label]
  require_columns(estimates, c("feature_id", "beta_adjusted", "se"), "transcript estimates")
  d <- merge(estimates[, .(feature_id, beta_adjusted, se)],
    annotations[, .(transcript_id, tss_te_label = promoter_te_label)], by.x = "feature_id",
    by.y = "transcript_id")
  d[, q := beta_adjusted^2 - se^2]
  a <- d[tss_te_label == "TE"]; b <- d[tss_te_label == "nonTE"]
  if (!nrow(a) || !nrow(b)) return(as.data.table(c(list(contrast = character()),
    setNames(rep(list(numeric()), length(NB_TE_METRICS)), NB_TE_METRICS))))
  pair <- CJ(te = a$feature_id, nonte = b$feature_id)
  ia <- match(pair$te, a$feature_id); ib <- match(pair$nonte, b$feature_id)
  pair[, `:=`(signed = a$beta_adjusted[ia] - b$beta_adjusted[ib],
    magnitude = abs(a$beta_adjusted[ia]) - abs(b$beta_adjusted[ib]),
    magnitude_sq_corrected = a$q[ia] - b$q[ib],
    magnitude_rms_corrected = nb_te_rms(a$q[ia]) - nb_te_rms(b$q[ib]))]
  pair[, contrast := paste(te, nonte, sep = "__vs__")]
  out <- pair[, c("contrast", NB_TE_METRICS), with = FALSE]
  qa <- mean(a$q); qb <- mean(b$q)
  out <- rbind(out, data.table(contrast = "class_mean",
    signed = mean(a$beta_adjusted) - mean(b$beta_adjusted),
    magnitude = mean(abs(a$beta_adjusted)) - mean(abs(b$beta_adjusted)),
    magnitude_sq_corrected = qa - qb,
    magnitude_rms_corrected = nb_te_rms(qa) - nb_te_rms(qb)))
  out[, `:=`(TE_signed_mean = mean(a$beta_adjusted), nonTE_signed_mean = mean(b$beta_adjusted),
    TE_magnitude_mean = mean(abs(a$beta_adjusted)), nonTE_magnitude_mean = mean(abs(b$beta_adjusted)),
    TE_sq_corrected_mean = qa, nonTE_sq_corrected_mean = qb,
    TE_mean_se = mean(a$se), nonTE_mean_se = mean(b$se),
    raw_magnitude_zero_effect_inflation_difference = sqrt(2 / pi) * (mean(a$se) - mean(b$se)))]
  out
}
nb_bootstrap_units <- function(z, public) {
  units <- unique(z[, .(sample_id, cluster_id, cross)])
  if (anyDuplicated(units$sample_id)) stop("Inconsistent sample bootstrap metadata")
  cl <- unique(units[, .(cluster_id, cross)])
  if (anyDuplicated(cl$cluster_id)) stop("Biological cluster spans reciprocal directions")
  strata <- if (public) split(cl$cluster_id, cl$cross) else list(all = cl$cluster_id)
  chosen <- unlist(lapply(strata, function(ids) ids[sample.int(length(ids), length(ids), replace = TRUE)]),
    use.names = FALSE)
  rbindlist(lapply(seq_along(chosen), function(j) {
    target <- chosen[j]
    part <- copy(z[cluster_id == target])
    part[, `:=`(sample_id = paste0(sample_id, "_boot", j), cluster_id = paste0("boot", j))]
    part
  }))
}
nb_gene_data <- function(s, transcript_ids) {
  map <- reg_read(nb_path(s, "transcript_feature_batches.tsv"))
  manifest <- reg_read(nb_path(s, "transcript_batches.tsv"))
  batches <- unique(map[feature_id %in% transcript_ids]$batch)
  rbindlist(lapply(batches, function(b) {
    row <- manifest[batch == b]
    if (nrow(row) != 1L || !file.exists(row$input_file)) stop("Within-gene input batch missing")
    inp <- readRDS(row$input_file)
    inp$data[feature_id %in% transcript_ids]
  }))
}
nb_te_gene <- function(cfg, s, root, job_id) {
  if (cfg$profile != "human_trio") stop("Use h1-inputs and h1-fit for reciprocal Bayesian H1")
  nb_check_prepare(cfg, s, root); r <- nb_te_require(cfg, s, root)
  jobs <- reg_read(nb_path(s, "within_gene_jobs.tsv"))
  job <- jobs[job_index == as.integer(job_id)]
  if (nrow(job) != 1L) stop("Unknown within-gene job index")
  gene <- job$gene_id
  all <- reg_read(nb_path(s, "transcript_TE_POE.tsv"))
  a <- nb_estimable(all[gene_id == gene & promoter_te_label %in% c("TE", "nonTE")])
  estimates <- copy(a); setnames(estimates, "transcript_id", "feature_id")
  observed <- nb_te_statistics(estimates, a)
  observed[, no_detected_cross_interaction := !any(a$not_strictly_POE %in% TRUE)]
  # Each bootstrap must estimate the same set of transcripts as the observed statistic.
  data <- nb_gene_data(s, a$transcript_id)
  key <- list(r$key, gene, s$bootstrap, s$seed, code = nb_code_version(root),
    statistics = NB_TE_STATS_VERSION)
  stage <- nb_te_stage(job_id)
  output <- nb_path(s, "within_gene", sprintf("gene_%06d.tsv", job_id))
  if (nb_current(s, stage, key)) return(invisible(TRUE))
  chunks <- split(seq_len(s$bootstrap), ceiling(seq_len(s$bootstrap) / s$bootstrap_chunk))
  chunk_files <- character()
  for (k in seq_along(chunks)) {
    chunk_stage <- nb_te_stage(job_id, k)
    file <- nb_path(s, "private", "bootstraps", NB_TE_STATS_VERSION,
      sprintf("gene_%06d_chunk_%05d.rds", job_id, k))
    chunk_files <- c(chunk_files, file)
    if (nb_current(s, chunk_stage, key)) next
    results <- lapply(chunks[[k]], function(draw_id) {
      seed <- (as.numeric(s$seed) + sum(utf8ToInt(gene) * seq_along(utf8ToInt(gene))) * 1009 + draw_id) %% .Machine$integer.max
      set.seed(as.integer(seed))
      sample <- nb_bootstrap_units(data, cfg$profile == "public_cross")
      fits <- rbindlist(lapply(split(sample, sample$feature_id), nb_fit_one, cfg = cfg, s = s), fill = TRUE)
      good <- nb_estimable(fits)
      # Same-transcript-set refits only; each refit uses its own estimates and SEs.
      # Keep only the contrast and metric columns so failed draws bind consistently.
      if (!setequal(good$feature_id, a$transcript_id)) {
        out <- observed[, .(contrast)]
        out[, (NB_TE_METRICS) := NA_real_]
      } else out <- nb_te_statistics(good, a)[, c("contrast", NB_TE_METRICS), with = FALSE]
      out[, draw := draw_id]; out
    })
    nb_save(rbindlist(results), file)
    nb_finish(s, chunk_stage, key, file)
  }
  boot <- rbindlist(lapply(chunk_files, readRDS))
  result <- copy(observed)
  for (j in seq_len(nrow(result))) for (metric in NB_TE_METRICS) {
    v <- boot[contrast == result$contrast[j]][[metric]]
    ok <- is.finite(v); frac <- mean(ok); obs <- result[[metric]][j]
    ci <- if (frac >= .9) quantile(v[ok], c(.025, .975), names = FALSE) else c(NA_real_, NA_real_)
    result[j, (paste0(metric, c("_bootstrap_se", "_bootstrap_ci_low", "_bootstrap_ci_high",
        "_bootstrap_centered_tail_prob", "_bootstrap_valid_fraction"))) :=
      list(if (frac >= .9) sd(v[ok]) else NA_real_, ci[1], ci[2],
        if (frac >= .9) (1 + sum(abs(v[ok] - obs) >= abs(obs))) / (1 + sum(ok)) else NA_real_, frac)]
  }
  result[, `:=`(gene_id = gene, cohort_id = cfg$cohort_id, version = NB_VERSION,
    n_TE = nrow(a[promoter_te_label == "TE"]), n_nonTE = nrow(a[promoter_te_label == "nonTE"]),
    any_not_strictly_POE = any(a$not_strictly_POE %in% TRUE),
    independence_verified = all(a$independence_verified),
    shared_marker_required = FALSE, same_first_donor_required = FALSE,
    interpretation = "descriptive_class_contrast_with_count_refit_bootstrap_diagnostic",
    magnitude_primary_field = "magnitude_sq_corrected", statistics_version = NB_TE_STATS_VERSION,
    te_tss_window_bp = reg_settings(cfg)$te_tss_window_bp)]
  # A separate sensitivity subset never replaces the primary transcript set.
  aa <- a[!not_strictly_POE %in% TRUE]
  if (any(aa$promoter_te_label == "TE") && any(aa$promoter_te_label == "nonTE")) {
    ee <- copy(aa); setnames(ee, "transcript_id", "feature_id")
    sensitivity <- nb_te_statistics(ee, aa)[contrast == "class_mean"]
    sensitivity[, `:=`(gene_id = gene, sensitivity = "no_detected_cross_interaction_descriptive")]
  } else sensitivity <- data.table(gene_id = gene, sensitivity = "no_eligible_TE_nonTE_sensitivity_pair")
  sf <- sub("\\.tsv$", "_sensitivity.tsv", output)
  nb_write(result, output); nb_write(sensitivity, sf)
  nb_finish(s, stage, key, c(output, sf))
}
nb_te_combine <- function(cfg, s, root) {
  nb_check_prepare(cfg, s, root); nb_te_require(cfg, s, root)
  jobs <- reg_read(nb_path(s, "within_gene_jobs.tsv"))
  for (j in jobs$job_index) nb_require(s, nb_te_stage(j))
  files <- nb_path(s, "within_gene", sprintf("gene_%06d.tsv", jobs$job_index))
  result <- if (nrow(jobs)) rbindlist(lapply(files, reg_read), fill = TRUE) else {
    metric_cols <- c(NB_TE_METRICS, paste0(NB_TE_METRICS, "_bootstrap_centered_tail_prob"))
    as.data.table(c(list(gene_id = character(), contrast = character()),
      setNames(rep(list(numeric()), length(metric_cols)), metric_cols),
      list(any_not_strictly_POE = logical())))
  }
  nb_write(result, nb_path(s, "within_gene_contrasts.tsv.gz"))
  class <- result[contrast == "class_mean"]
  summarize <- function(d, name) {
    out <- data.table(subset = name, n_genes = nrow(d),
      status = if (nrow(d) >= 2) "exploratory_gene_bootstrap" else "insufficient_genes")
    for (metric in NB_TE_METRICS) {
      v <- d[[metric]]
      if (length(v)) out[, (paste0(metric, "_mean")) := mean(v)]
      if (length(v) >= 2) {
        set.seed(s$seed)
        boot <- replicate(s$bootstrap, mean(v[sample.int(length(v), length(v), replace = TRUE)]))
        ci <- quantile(boot, c(.025, .975), names = FALSE)
        out[, (paste0(metric, c("_descriptive_ci_low", "_descriptive_ci_high"))) := list(ci[1], ci[2])]
      }
    }
    out
  }
  summary <- rbindlist(list(summarize(class, "all_estimable_mixed_genes"),
    summarize(class[any_not_strictly_POE %in% FALSE], "no_detected_cross_interaction")), fill = TRUE)
  summary[, `:=`(cohort_id = cfg$cohort_id, unit = "gene_equal_weight",
    magnitude_primary_field = "magnitude_sq_corrected", statistics_version = NB_TE_STATS_VERSION,
    caveat = "gene_bootstrap_conditions_on_fitted_effects_not_independent_replication")]
  nb_write(summary, nb_path(s, "within_gene_cohort_summary.tsv"))
  all <- reg_read(nb_path(s, "transcript_TE_POE.tsv"))
  coverage <- all[, .(n_transcripts = .N,
    n_TE = sum(promoter_te_label == "TE"), n_nonTE = sum(promoter_te_label == "nonTE")),
    by = gene_id]
  coverage[, diagnostic_status := "no_estimable_mixed_classes"]
  if (nrow(class)) {
    valid <- as.matrix(class[, paste0(NB_TE_METRICS, "_bootstrap_valid_fraction"), with = FALSE])
    ok <- data.table(gene_id = class$gene_id, ok = rowSums(is.na(valid) | valid < .9) == 0)
    coverage[ok, on = "gene_id", diagnostic_status := fifelse(i.ok,
      "bootstrap_diagnostic_available", "bootstrap_diagnostic_insufficient_success")]
  }
  nb_write(coverage, nb_path(s, "h1_coverage.tsv"))
  nb_finish(s, "te-combine", list(nb_require(s, "te-inputs")$key,
    results = if (length(files)) nb_file_key(files) else "no_mixed_genes"),
    nb_path(s, c("within_gene_contrasts.tsv.gz", "within_gene_cohort_summary.tsv", "h1_coverage.tsv")))
}
nb_between_genes <- function(cfg, s, root) {
  nb_check_prepare(cfg, s, root); nb_te_require(cfg, s, root)
  z <- reg_settings(cfg); z$output_dir <- nb_path(s, "between_genes")
  availability <- reg_read(nb_path(s, "TE_availability.tsv"))
  if (availability$repeats != "full_reference") {
    out <- nb_path(s, "between_genes", "te_genes_result.tsv")
    nb_write(data.table(cohort_id = cfg$cohort_id, status = "unavailable_full_repeat_reference_required",
      version = NB_VERSION), out)
    nb_finish(s, "between-genes", nb_require(s, "te-inputs")$key, out)
    return(invisible(TRUE))
  }
  # The adjusted comparison uses all-SNP gene POE and its count inventory.
  inp <- reg_te_genes_inputs(cfg, reg_settings(cfg), root, reg_te_genes_settings(z),
    baseline = list(s = s))
  outputs <- file.path(z$output_dir, TE_GENES_FILES)
  reg_te_genes_run(cfg, z, inp, reg_te_genes_settings(z), outputs)
  if (cfg$profile == "public_cross") {
    sensitivity <- inp
    sensitivity$genes <- sensitivity$genes[!not_strictly_POE %in% TRUE]
    sensitivity_outputs <- file.path(z$output_dir, "no_detected_cross_interaction", TE_GENES_FILES)
    reg_te_genes_run(cfg, z, sensitivity, reg_te_genes_settings(z), sensitivity_outputs)
    outputs <- c(outputs, sensitivity_outputs)
  }
  nb_finish(s, "between-genes", nb_require(s, "te-inputs")$key, outputs)
}
nb_transcript_between <- function(cfg, s) {
  x <- nb_estimable(reg_read(nb_path(s, "transcript_TE_POE.tsv")))
  x <- x[promoter_te_label %in% c("TE", "nonTE")]
  x[, tss_te_label := promoter_te_label]
  # One class mean per gene, so isoform-rich genes do not supply extra replicates.
  # magnitude_sq_corrected removes noise inflation (b^2 - SE^2); see nb_te_statistics.
  g <- x[, .(signed = mean(beta_adjusted),
    magnitude = mean(abs(beta_adjusted)), magnitude_sq_corrected = mean(beta_adjusted^2 - se^2),
    n_transcripts = .N,
    any_not_strictly_POE = any(not_strictly_POE %in% TRUE)), by = .(gene_id, tss_te_label)]
  out <- data.table(cohort_id = cfg$cohort_id, n_genes = uniqueN(g$gene_id),
    status = "insufficient_groups", estimand = "transcript_class_gene_equal_weight_not_gene_aggregate_POE")
  if (all(c("TE", "nonTE") %in% g$tss_te_label)) {
    difference <- function(a, metric) mean(a[tss_te_label == "TE"][[metric]]) -
      mean(a[tss_te_label == "nonTE"][[metric]])
    blocks <- split(g, g$gene_id)
    set.seed(s$seed)
    # The added metric runs last, so the existing signed/magnitude draws are unchanged.
    for (metric in c("signed", "magnitude", "magnitude_sq_corrected")) {
      obs <- difference(g, metric)
      boot <- replicate(s$bootstrap, {
        a <- rbindlist(blocks[sample.int(length(blocks), length(blocks), replace = TRUE)])
        if (all(c("TE", "nonTE") %in% a$tss_te_label)) difference(a, metric) else NA_real_
      })
      ok <- is.finite(boot)
      ci <- if (mean(ok) >= .9) quantile(boot[ok], c(.025, .975), names = FALSE) else c(NA, NA)
      out[, (paste0(metric, c("_difference", "_ci_low", "_ci_high"))) := list(obs, ci[1], ci[2])]
    }
    out$status <- "descriptive_gene_block_bootstrap"
  }
  nb_write(g, nb_path(s, "transcript_between_gene_table.tsv"))
  nb_write(out, nb_path(s, "transcript_between_gene_summary.tsv"))
}
