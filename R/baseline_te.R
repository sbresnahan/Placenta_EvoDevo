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
    if (file.exists(p)) reg_sha(p) else "not_present", ""), candidates)
  nb_hash(list(nb_require(s, "combine")$key, resources = state,
    reference = reg_sha(cfg$paths$tx_map), code = nb_code_hash(root)))
}
nb_te_require <- function(cfg, s, root) {
  nb_require(s, "te-inputs", nb_te_resource_key(cfg, s, root))
}
nb_te_annotation <- function(tx, repeats = NULL, full_reference = FALSE) {
  x <- copy(tx)
  require_columns(x, c("transcript_id", "gene_id", "chrom", "tss"))
  x[, `:=`(tss_te_label = "not_assessed", te_class = NA_character_,
    te_family = NA_character_, te_ids = NA_character_)]
  if (is.null(repeats) || !nrow(repeats)) return(x)
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
  x
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
  tx <- nb_te_annotation(reg_read(cfg$paths$tx_map), repeats, full)
  nb_write(tx, nb_path(s, "transcript_TE_annotation.tsv"))
  fits <- reg_read(nb_path(s, "transcript_POE.tsv"))
  d <- merge(tx, fits, by.x = "transcript_id", by.y = "feature_id", all.x = TRUE)
  d[is.na(status), status := "no_informative_counts"]
  nb_write(d, nb_path(s, "transcript_TE_POE.tsv"))
  eligible <- nb_estimable(d)
  mixed <- eligible[tss_te_label %in% c("TE", "nonTE"), .(
    n_TE = sum(tss_te_label == "TE"), n_nonTE = sum(tss_te_label == "nonTE")), by = gene_id]
  mixed <- mixed[n_TE > 0 & n_nonTE > 0]
  setorder(mixed, gene_id); mixed[, job_index := .I]
  nb_write(mixed, nb_path(s, "within_gene_jobs.tsv"))
  availability <- data.table(cohort_id = cfg$cohort_id,
    repeats = if (is.null(repeats)) "not_available" else if (full) "full_reference" else "partial_reference",
    n_transcripts = nrow(tx), n_estimable_transcripts = nrow(eligible),
    n_TE = sum(tx$tss_te_label == "TE"), n_nonTE = sum(tx$tss_te_label == "nonTE"),
    n_unknown_or_unassessed = sum(!tx$tss_te_label %in% c("TE", "nonTE")), n_mixed_genes = nrow(mixed),
    status = if (is.null(repeats)) "TE_comparisons_unavailable_missing_repeats" else
      if (!full) "TE_comparisons_limited_partial_repeat_reference" else "ready")
  nb_write(availability, nb_path(s, "TE_availability.tsv"))
  nb_finish(s, "te-inputs", key, nb_path(s, c("transcript_TE_annotation.tsv",
    "transcript_TE_POE.tsv", "within_gene_jobs.tsv", "TE_availability.tsv")))
}
nb_te_statistics <- function(estimates, annotations) {
  d <- merge(estimates[, .(feature_id, beta_adjusted)],
    annotations[, .(transcript_id, tss_te_label)], by.x = "feature_id",
    by.y = "transcript_id")
  a <- d[tss_te_label == "TE"]; b <- d[tss_te_label == "nonTE"]
  if (!nrow(a) || !nrow(b)) return(data.table(contrast = character(), signed = numeric(), magnitude = numeric()))
  pair <- CJ(te = a$feature_id, nonte = b$feature_id)
  pair[, `:=`(signed = a$beta_adjusted[match(te, a$feature_id)] -
    b$beta_adjusted[match(nonte, b$feature_id)],
    magnitude = abs(a$beta_adjusted[match(te, a$feature_id)]) -
      abs(b$beta_adjusted[match(nonte, b$feature_id)]))]
  pair[, contrast := paste(te, nonte, sep = "__vs__")]
  out <- pair[, .(contrast, signed, magnitude)]
  rbind(out, data.table(contrast = "class_mean",
    signed = mean(a$beta_adjusted) - mean(b$beta_adjusted),
    magnitude = mean(abs(a$beta_adjusted)) - mean(abs(b$beta_adjusted))))
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
    if (nrow(row) != 1L || !identical(reg_sha(row$input_file), row$sha256))
      stop("Within-gene input batch missing or changed")
    inp <- readRDS(row$input_file)
    if (!identical(inp$key, nb_require(s, "prepare")$key)) stop("Mixed within-gene inputs")
    inp$data[feature_id %in% transcript_ids]
  }))
}
nb_te_gene <- function(cfg, s, root, job_id) {
  nb_check_prepare(cfg, s, root); r <- nb_te_require(cfg, s, root)
  jobs <- reg_read(nb_path(s, "within_gene_jobs.tsv"))
  job <- jobs[job_index == as.integer(job_id)]
  if (nrow(job) != 1L) stop("Unknown within-gene job index")
  gene <- job$gene_id
  all <- reg_read(nb_path(s, "transcript_TE_POE.tsv"))
  a <- nb_estimable(all[gene_id == gene & tss_te_label %in% c("TE", "nonTE")])
  estimates <- copy(a); setnames(estimates, "transcript_id", "feature_id")
  observed <- nb_te_statistics(estimates, a)
  observed[, no_detected_cross_interaction := !any(a$not_strictly_POE %in% TRUE)]
  # Each bootstrap must estimate the same set of transcripts as the observed statistic.
  data <- nb_gene_data(s, a$transcript_id)
  key <- nb_hash(list(r$key, gene, s$bootstrap, s$seed, code = nb_code_hash(root)))
  stage <- sprintf("te_%06d", job_id)
  output <- nb_path(s, "within_gene", sprintf("gene_%06d.tsv", job_id))
  if (nb_current(s, stage, key)) return(invisible(TRUE))
  chunks <- split(seq_len(s$bootstrap), ceiling(seq_len(s$bootstrap) / s$bootstrap_chunk))
  chunk_files <- character()
  for (k in seq_along(chunks)) {
    chunk_stage <- sprintf("te_%06d_chunk_%05d", job_id, k)
    file <- nb_path(s, "private", "bootstraps", sprintf("gene_%06d_chunk_%05d.rds", job_id, k))
    chunk_files <- c(chunk_files, file)
    if (nb_current(s, chunk_stage, key)) next
    results <- lapply(chunks[[k]], function(draw_id) {
      seed <- (as.numeric(s$seed) + sum(utf8ToInt(gene) * seq_along(utf8ToInt(gene))) * 1009 + draw_id) %% .Machine$integer.max
      set.seed(as.integer(seed))
      sample <- nb_bootstrap_units(data, cfg$profile == "public_cross")
      fits <- rbindlist(lapply(split(sample, sample$feature_id), nb_fit_one, cfg = cfg, s = s), fill = TRUE)
      good <- nb_estimable(fits)
      if (!setequal(good$feature_id, a$transcript_id)) {
        out <- observed[, .(contrast, signed = NA_real_, magnitude = NA_real_)]
      } else out <- nb_te_statistics(good, a)
      out[, draw := draw_id]; out
    })
    nb_save(rbindlist(results), file)
    nb_finish(s, chunk_stage, key, file)
  }
  boot <- rbindlist(lapply(chunk_files, readRDS))
  result <- copy(observed)
  for (j in seq_len(nrow(result))) for (metric in c("signed", "magnitude")) {
    v <- boot[contrast == result$contrast[j]][[metric]]
    ok <- is.finite(v); frac <- mean(ok); obs <- result[[metric]][j]
    ci <- if (frac >= .9) quantile(v[ok], c(.025, .975), names = FALSE) else c(NA_real_, NA_real_)
    result[j, (paste0(metric, c("_se", "_ci_low", "_ci_high", "_p", "_valid_fraction"))) :=
      list(if (frac >= .9) sd(v[ok]) else NA_real_, ci[1], ci[2],
        if (frac >= .9) (1 + sum(abs(v[ok] - obs) >= abs(obs))) / (1 + sum(ok)) else NA_real_, frac)]
  }
  result[, `:=`(gene_id = gene, cohort_id = cfg$cohort_id, version = NB_VERSION,
    n_TE = nrow(a[tss_te_label == "TE"]), n_nonTE = nrow(a[tss_te_label == "nonTE"]),
    any_not_strictly_POE = any(a$not_strictly_POE %in% TRUE),
    independence_verified = all(a$independence_verified),
    shared_marker_required = FALSE, same_first_donor_required = FALSE,
    interpretation = "within_gene_baseline_NB_bootstrap_magnitude_exploratory")]
  # A separate sensitivity subset never replaces the primary transcript set.
  aa <- a[!not_strictly_POE %in% TRUE]
  if (any(aa$tss_te_label == "TE") && any(aa$tss_te_label == "nonTE")) {
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
  for (j in jobs$job_index) nb_require(s, sprintf("te_%06d", j))
  files <- nb_path(s, "within_gene", sprintf("gene_%06d.tsv", jobs$job_index))
  result <- if (nrow(jobs)) rbindlist(lapply(files, reg_read), fill = TRUE) else
    data.table(gene_id = character(), contrast = character(), signed = numeric(), magnitude = numeric(),
      signed_p = numeric(), magnitude_p = numeric(), any_not_strictly_POE = logical())
  result[, `:=`(signed_q = p.adjust(signed_p, "BH"), magnitude_q = p.adjust(magnitude_p, "BH"))]
  nb_write(result, nb_path(s, "within_gene_contrasts.tsv.gz"))
  class <- result[contrast == "class_mean"]
  summarize <- function(d, name) {
    out <- data.table(subset = name, n_genes = nrow(d),
      status = if (nrow(d) >= 2) "exploratory_gene_bootstrap" else "insufficient_genes")
    for (metric in c("signed", "magnitude")) {
      v <- d[[metric]]
      if (length(v)) out[, (paste0(metric, "_mean")) := mean(v)]
      if (length(v) >= 2) {
        set.seed(s$seed)
        boot <- replicate(s$bootstrap, mean(v[sample.int(length(v), length(v), replace = TRUE)]))
        ci <- quantile(boot, c(.025, .975), names = FALSE)
        out[, (paste0(metric, c("_ci_low", "_ci_high", "_p"))) :=
          list(ci[1], ci[2], (1 + sum(abs(boot - mean(v)) >= abs(mean(v)))) / (1 + length(boot)))]
      }
    }
    out
  }
  summary <- rbindlist(list(summarize(class, "all_estimable_mixed_genes"),
    summarize(class[any_not_strictly_POE %in% FALSE], "no_detected_cross_interaction")), fill = TRUE)
  summary[, `:=`(cohort_id = cfg$cohort_id, unit = "gene_equal_weight",
    caveat = "gene_bootstrap_conditions_on_fitted_effects_not_independent_replication")]
  nb_write(summary, nb_path(s, "within_gene_cohort_summary.tsv"))
  nb_finish(s, "te-combine", nb_hash(list(nb_require(s, "te-inputs")$key,
    results = if (length(files)) nb_file_key(files) else "no_mixed_genes")),
    nb_path(s, c("within_gene_contrasts.tsv.gz", "within_gene_cohort_summary.tsv")))
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
  # Existing adjusted comparison is reused, but its outcome and depth inventory
  # must come from all-SNP baseline outputs, never legacy completed_counts().
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
  x <- x[tss_te_label %in% c("TE", "nonTE")]
  # One class mean per gene, so isoform-rich genes do not supply extra replicates.
  g <- x[, .(signed = mean(beta_adjusted),
    magnitude = mean(abs(beta_adjusted)), n_transcripts = .N,
    any_not_strictly_POE = any(not_strictly_POE %in% TRUE)), by = .(gene_id, tss_te_label)]
  out <- data.table(cohort_id = cfg$cohort_id, n_genes = uniqueN(g$gene_id),
    status = "insufficient_groups", estimand = "transcript_class_gene_equal_weight_not_gene_aggregate_POE")
  if (all(c("TE", "nonTE") %in% g$tss_te_label)) {
    difference <- function(a, metric) mean(a[tss_te_label == "TE"][[metric]]) -
      mean(a[tss_te_label == "nonTE"][[metric]])
    blocks <- split(g, g$gene_id)
    set.seed(s$seed)
    for (metric in c("signed", "magnitude")) {
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
