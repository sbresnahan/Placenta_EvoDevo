reg_annotate_events <- function(tx, event_map, repeats, seqinfo, aliases = NULL,
    full_reference = TRUE) {
  require_columns(tx, c("transcript_id", "gene_id", "chrom", "strand", "tss",
    "first_start", "first_end"))
  require_columns(event_map, c("transcript_id", "event_id", "window"))
  require_columns(seqinfo, c("chrom", "length"))
  assert_unique(tx, "transcript_id", "Transcript map")
  assert_unique(seqinfo, "chrom", "Reference contigs")
  if (any(!tx$chrom %in% seqinfo$chrom)) stop("Transcript contigs absent from reference .fai")
  if (anyNA(tx[, .(tss, first_start, first_end)]) ||
      any(!tx$strand %in% c("+", "-"))) stop("Invalid transcript TSS/exon structure")
  limits <- seqinfo$length[match(tx$chrom, seqinfo$chrom)]
  if (any(tx$tss < 1 | tx$tss > limits | tx$first_end > limits))
    stop("Transcript coordinates exceed reference contigs")
  r <- copy(repeats)
  r[, chrom := reg_alias(chrom, seqinfo$chrom, aliases)]
  unmatched <- sum(is.na(r$chrom))
  r <- r[!is.na(chrom)]
  if (!nrow(r)) stop("No repeats map to this reference; check assembly/chromosome names")
  if (any(r$end > seqinfo$length[match(r$chrom, seqinfo$chrom)]))
    stop("Repeat coordinates exceed reference contigs")
  tx <- copy(tx)
  tx[, `:=`(tss_te = FALSE, tss_unknown_repeat = FALSE,
    tss_repeat_overlap = FALSE, te_class = NA_character_, te_family = NA_character_,
    te_ids = NA_character_, te_orientation = NA_character_,
    first_exon_te = FALSE, window_te = FALSE)]
  h <- reg_hits(tx[, .(chrom, start = tss, end = tss)], r)
  if (nrow(h)) {
    lab <- cbind(h, r[h$subject, .(is_te, is_unknown, repeat_class,
      repeat_family, repeat_id, repeat_strand)])
    lab[, orientation := ifelse(repeat_strand == tx$strand[query], "sense", "antisense")]
    a <- lab[, .(tss_te = any(is_te), tss_unknown_repeat = any(is_unknown),
      te_class = paste(sort(unique(repeat_class[is_te])), collapse = ";"),
      te_family = paste(sort(unique(repeat_family[is_te])), collapse = ";"),
      te_ids = paste(sort(unique(repeat_id[is_te])), collapse = ";"),
      te_orientation = paste(sort(unique(orientation[is_te])), collapse = ";")), by = query]
    tx[a$query, `:=`(tss_te = a$tss_te, tss_unknown_repeat = a$tss_unknown_repeat,
      tss_repeat_overlap = TRUE, te_class = a$te_class, te_family = a$te_family,
      te_ids = a$te_ids, te_orientation = a$te_orientation)]
  }
  te <- r[is_te %in% TRUE]
  h <- reg_hits(tx[, .(chrom, start = first_start, end = first_end)], te)
  tx[unique(h$query), first_exon_te := TRUE]
  win <- tx[, .(chrom, start = pmax(1, ifelse(strand == "+", tss - 500, tss - 100)),
    end = pmin(limits, ifelse(strand == "+", tss + 100, tss + 500)))]
  h <- reg_hits(win, te); tx[unique(h$query), window_te := TRUE]
  tx[, tss_te_label := fifelse(tss_unknown_repeat, "unknown_repeat",
    fifelse(tss_te, "TE", if (full_reference) "nonTE" else "not_assessed"))]
  m <- unique(event_map[window == 50L, .(transcript_id, event_id)])
  if (any(!m$transcript_id %in% tx$transcript_id)) stop("Event map refers to absent transcripts")
  t <- merge(m, tx, by = "transcript_id")
  events <- t[, .(gene_id = gene_id[1], chrom = chrom[1], strand = strand[1],
    tss_te_label = if (any(tss_te_label == "unknown_repeat")) "unknown_repeat" else
      if (any(tss_te_label == "not_assessed")) "not_assessed" else
      if (all(tss_te_label == "TE")) "TE" else
      if (all(tss_te_label == "nonTE")) "nonTE" else "mixed",
    n_transcripts = .N, te_ids = paste(sort(unique(te_ids[!is.na(te_ids)])), collapse = ";"),
    te_class = paste(sort(unique(te_class[!is.na(te_class)])), collapse = ";"),
    first_exon_te = any(first_exon_te), window_te = any(window_te)), by = event_id]
  if (any(t[, uniqueN(gene_id), by = event_id]$V1 != 1L))
    stop("Event spans multiple genes")
  if (any(t[, uniqueN(paste(chrom, strand)), by = event_id]$V1 != 1L))
    stop("Event spans chromosomes/strands")
  list(transcripts = t, events = events, repeats = r,
    coverage = data.table(n_reference_contigs = nrow(seqinfo),
      n_transcripts = nrow(tx), unmatched_repeat_records = unmatched,
      n_events = nrow(events), n_TE = sum(events$tss_te_label == "TE"),
      n_nonTE = sum(events$tss_te_label == "nonTE"),
      n_unknown = sum(events$tss_te_label == "unknown_repeat"),
      n_mixed = sum(events$tss_te_label == "mixed"), full_reference = full_reference))
}
reg_te_pairs <- function(events, labels) {
  require_columns(events, c("window", "event_id", "gene_id", "first_donor",
    "tss_min", "distinguishing_bp"))
  # Same-first-donor TSS events have nested first exons: the downstream event's
  # first exon lies inside the upstream one, so it never has distinguishing
  # sequence. Requiring distinguishing_bp > 0 for both events therefore excluded
  # every pair. Distinguishing sequence is kept as support_* and
  # direct_coordinate_both; counts are Ornaments allocations (model_scope).
  e <- copy(events[window == 50L])
  assert_unique(e, "event_id", "50bp event catalog")
  e <- merge(e, labels[, .(event_id, tss_te_label, chrom, strand)],
    by = "event_id", all.x = TRUE, suffixes = c("_catalog", ""))
  a <- e[tss_te_label == "TE" & !is.na(first_donor),
    .(gene_id, chrom, strand, first_donor, event_te = event_id, tss_te = tss_min,
      support_te = distinguishing_bp)]
  b <- e[tss_te_label == "nonTE" & !is.na(first_donor),
    .(gene_id, chrom, strand, first_donor, event_nonTE = event_id, tss_nonTE = tss_min,
      support_nonTE = distinguishing_bp)]
  p <- merge(a, b, by = c("gene_id", "chrom", "strand", "first_donor"), allow.cartesian = TRUE)
  p[, tss_distance_bp := abs(tss_te - tss_nonTE)]
  p <- p[tss_distance_bp > 50]
  setorder(p, gene_id, -support_te, -support_nonTE, event_te, event_nonTE)
  p
}
reg_fit_two <- function(Y, X, cluster) {
  if (any(!is.finite(Y)) || anyNA(X) || anyNA(cluster)) return(NULL)
  n <- nrow(X); p <- ncol(X); G <- uniqueN(cluster)
  if (qr(X)$rank != p || n <= p || G < 2L) return(NULL)
  fit <- lm.fit(X, Y); bread <- solve(crossprod(X))
  lev <- rowSums((X %*% bread) * X)
  if (any(lev >= 1 - 1e-10)) return(NULL)
  diff_res <- fit$residuals[, 1] - fit$residuals[, 2]
  if (G == n) {
    score <- X * (diff_res / (1 - lev))
    variance <- bread %*% crossprod(score) %*% bread
  } else {
    score <- rowsum(X * diff_res, cluster, reorder = FALSE)
    variance <- (G / (G - 1)) * ((n - 1) / (n - p)) *
      (bread %*% crossprod(score) %*% bread)
  }
  se <- sqrt(max(0, variance[1, 1]))
  list(beta = fit$coefficients[1, ], se = se, df = min(n - p, G - 1L),
    n_clusters = G)
}
reg_pair_inference <- function(Y, X, cluster, strata, B, seed) {
  fit <- reg_fit_two(Y, X, cluster)
  if (is.null(fit)) return(list(status = "unidentifiable_rank_clusters_or_leverage"))
  beta <- fit$beta
  signed <- beta[1] - beta[2]; magnitude <- abs(beta[1]) - abs(beta[2])
  # Resample whole biological clusters jointly for both events. Reciprocal
  # orientation strata are preserved; human families have one common stratum.
  cs <- unique(data.table(cluster = cluster, stratum = strata))
  if (anyDuplicated(cs$cluster)) return(list(status = "cluster_crosses_reciprocal_strata"))
  by_cluster <- split(seq_len(nrow(Y)), cluster)
  groups <- split(cs$cluster, cs$stratum)
  set.seed(seed)
  boot <- replicate(B, {
    picked <- unlist(lapply(groups, function(ids) ids[sample.int(length(ids),
      length(ids), replace = TRUE)]), use.names = FALSE)
    ii <- unlist(by_cluster[picked], use.names = FALSE)
    xb <- X[ii, , drop = FALSE]
    if (qr(xb)$rank < ncol(xb)) c(NA_real_, NA_real_) else {
      bb <- lm.fit(xb, Y[ii, , drop = FALSE])$coefficients[1, ]
      c(bb[1] - bb[2], abs(bb[1]) - abs(bb[2]))
    }
  })
  ok <- apply(boot, 2, function(v) all(is.finite(v)))
  valid <- mean(ok)
  ci <- if (valid >= .9) t(apply(boot[, ok, drop = FALSE], 1, quantile,
    probs = c(.025, .975))) else matrix(NA_real_, 2, 2)
  pmag <- if (valid >= .9)
    (1 + sum(abs(boot[2, ok] - magnitude) >= abs(magnitude))) / (1 + sum(ok)) else NA_real_
  list(status = if (valid < .9) "unstable_bootstrap" else
      if (fit$se <= sqrt(.Machine$double.eps)) "degenerate_signed_variance" else
      "eligible_conditional_inference",
    estimate_te_log2 = unname(beta[1]), estimate_nonTE_log2 = unname(beta[2]),
    signed_difference_log2 = unname(signed), signed_se = fit$se,
    signed_ci_low = ci[1, 1], signed_ci_high = ci[1, 2],
    signed_p_value = if (fit$se > sqrt(.Machine$double.eps))
      2 * pt(-abs(signed / fit$se), fit$df) else NA_real_,
    magnitude_difference_log2 = unname(magnitude),
    magnitude_ci_low = ci[2, 1], magnitude_ci_high = ci[2, 2],
    magnitude_p_value = pmag, bootstrap_valid_fraction = valid,
    residual_df = fit$df, n_clusters = fit$n_clusters,
    magnitude_test = "centered_paired_bootstrap_exploratory_nonregular_near_zero")
}
reg_one_te_pair <- function(pair, event_counts, gene_counts, units, cov, cfg, z) {
  human <- cfg$profile == "human_trio"
  id <- if (human) "sample_id" else "unit_id"
  keys <- c(id, "locus_id", if (human) "family_id")
  require_columns(event_counts, c(keys, "event_id", "maternal_count", "paternal_count"))
  e <- copy(event_counts[event_id %in% c(pair$event_te, pair$event_nonTE)])
  if (any(e$maternal_count < 0 | e$paternal_count < 0, na.rm = TRUE))
    stop("Negative expected counts")
  assert_unique(e, c(id, "event_id", "locus_id"), "Event allocation keys")
  e[, depth := maternal_count + paternal_count]
  e <- e[is.finite(depth) & depth >= z$min_depth &
    is.finite(maternal_count) & is.finite(paternal_count)]
  a <- e[event_id == pair$event_te]; b <- e[event_id == pair$event_nonTE]
  d <- merge(a, b, by = keys, suffixes = c("_te", "_nonTE"))
  out <- as.list(pair)
  out <- c(out, list(status = "no_shared_callable_locus", n_units = 0L,
    min_pair_depth = NA_real_, locus_id = NA_character_))
  if (!nrow(d)) return(as.data.table(out))
  if (human) require_columns(d, "family_id")
  d[, biological_id := as.character(if (human) family_id else get(id))]
  scored <- d[, .(n = uniqueN(biological_id),
    weakest_depth = min(depth_te, depth_nonTE)), by = locus_id]
  setorder(scored, -n, -weakest_depth, locus_id)
  d <- d[locus_id == scored$locus_id[1]]
  setorderv(d, c("biological_id", "depth_te", "depth_nonTE", id), c(1, -1, -1, 1))
  d <- d[!duplicated(biological_id)]
  out$n_units <- nrow(d); out$locus_id <- scored$locus_id[1]
  out$min_pair_depth <- min(d$depth_te, d$depth_nonTE)
  if (human) {
    ori <- gene_counts[gene_id == pair$gene_id,
      .(ref_maternal = {
        v <- unique(na.omit(ref_maternal)); if (length(v) == 1L) v[1] else NA
      }), by = .(sample_id, locus_id)]
    d <- merge(d, ori, by = c("sample_id", "locus_id"), all.x = TRUE)
    if (anyNA(d$ref_maternal)) {
      out$status <- "missing_or_discordant_cis_orientation"; return(as.data.table(out))
    }
    o <- ifelse(d$ref_maternal, 1, -1)
    out$n_ref_maternal <- sum(o == 1); out$n_alt_maternal <- sum(o == -1)
    if (nrow(d) < cfg$models$min_families ||
        min(out$n_ref_maternal, out$n_alt_maternal) < cfg$models$min_orientation) {
      out$status <- "insufficient_families_or_cis_orientation"; return(as.data.table(out))
    }
    require_columns(cov, c("sample_id", unlist(cfg$models$covariates)))
    assert_unique(cov, "sample_id", "Covariates")
    cv <- copy(cov)[match(d$sample_id, sample_id)]
    vars <- unlist(cfg$models$covariates)
    dat <- data.frame(cis = o)
    for (v in vars) {
      dat[[v]] <- cv[[v]]
      if (is.numeric(dat[[v]])) dat[[v]] <- dat[[v]] - mean(dat[[v]])
    }
    if (anyNA(dat)) {out$status <- "missing_covariates"; return(as.data.table(out))}
    X <- tryCatch(model.matrix(~., dat), error = function(e) NULL)
    if (is.null(X)) {out$status <- "invalid_covariates"; return(as.data.table(out))}
    other <- setdiff(colnames(X), c("(Intercept)", "cis"))
    if (length(other)) X[, other] <- scale(X[, other, drop = FALSE], center = TRUE, scale = FALSE)
    cluster <- as.character(d$family_id); strata <- rep("human", nrow(d))
  } else {
    u <- units[match(d$unit_id, unit_id)]
    if (anyNA(u$maternal_parent) || any(!u$maternal_parent %in%
        c(cfg$parents$focal, cfg$parents$alternate))) stop("Unmatched reciprocal parental labels")
    o <- ifelse(u$maternal_parent == cfg$parents$focal, 1, -1)
    out$n_focal_mother <- sum(o == 1); out$n_alternate_mother <- sum(o == -1)
    if (length(unique(o)) < 2L) {
      out$status <- "single_direction_cis_POE_confounded"; return(as.data.table(out))
    }
    if (nrow(d) < reg_default(cfg$models$min_callable_units, 3L)) {
      out$status <- "insufficient_residual_units"; return(as.data.table(out))
    }
    if (anyNA(u$cluster_id) || any(!nzchar(u$cluster_id))) {
      out$status <- "missing_biological_cluster_ids"; return(as.data.table(out))
    }
    X <- cbind("(Intercept)" = 1, cis = o)
    cluster <- as.character(u$cluster_id); strata <- as.character(o)
  }
  offset <- log2(cfg$dosage$maternal / cfg$dosage$paternal)
  if (!is.finite(offset)) stop("Invalid parental dosage")
  Y <- cbind(log2((d$maternal_count_te + .5) / (d$paternal_count_te + .5)) - offset,
    log2((d$maternal_count_nonTE + .5) / (d$paternal_count_nonTE + .5)) - offset)
  support_cols <- c("distinguishing_coordinate_te", "distinguishing_coordinate_nonTE")
  out$direct_coordinate_both <- all(support_cols %in% names(d)) &&
    all(d$distinguishing_coordinate_te %in% TRUE & d$distinguishing_coordinate_nonTE %in% TRUE)
  seed <- as.integer((as.numeric(z$seed) + sum(utf8ToInt(
    paste(pair$event_te, pair$event_nonTE)))) %% .Machine$integer.max)
  ans <- reg_pair_inference(Y, X, cluster, strata, z$bootstrap, seed)
  out[names(ans)] <- ans
  as.data.table(out)
}
reg_analyze_pairs <- function(pairs, event_counts, gene_counts, units, cov, cfg, z) {
  template <- copy(pairs)[0, status := character()]
  for (v in c("estimate_te_log2", "estimate_nonTE_log2", "signed_difference_log2",
      "magnitude_difference_log2", "signed_p_value", "magnitude_p_value",
      "signed_q_value", "magnitude_q_value"))
    template[, (v) := numeric()]
  if (!nrow(pairs)) return(template)
  by_event <- split(event_counts, event_counts$event_id)
  ans <- lapply(seq_len(nrow(pairs)), function(i) {
    p <- pairs[i]
    # Prepend the empty schema so pairs whose events have no allelic counts
    # still carry the columns (a pair with no counts is scored, not an error).
    counts <- rbindlist(c(list(event_counts[0]), by_event[c(p$event_te, p$event_nonTE)]),
      fill = TRUE)
    reg_one_te_pair(p, counts, gene_counts, units, cov, cfg, z)
  })
  x <- rbindlist(ans, fill = TRUE)
  for (col in c("signed_p_value", "magnitude_p_value"))
    if (!col %in% names(x)) x[, (col) := NA_real_]
  setorder(x, gene_id, -n_units, -min_pair_depth, event_te, event_nonTE)
  x[, pair_rank_by_callability := seq_len(.N), by = gene_id]
  x[, selected_for_inference := pair_rank_by_callability <= z$max_pairs_per_gene]
  x[, `:=`(signed_q_value = NA_real_, magnitude_q_value = NA_real_)]
  use <- which(x$selected_for_inference & x$status == "eligible_conditional_inference")
  x[use, `:=`(signed_q_value = p.adjust(signed_p_value, "BH"),
    magnitude_q_value = p.adjust(magnitude_p_value, "BH"))]
  x[selected_for_inference == FALSE, status := paste0("outside_callability_cap:", status)]
  x[, `:=`(effect_scale = "log2_maternal_paternal_dosage_adjusted",
    evidence = "expected_counts_same_locus_joint_unit_bootstrap_not_independent_molecules",
    bh_family = "cohort_selected_pairs_separate_signed_and_magnitude",
    effect_type = "TSS_conditional_on_first_donor")]
  x
}
reg_cohort_summary <- function(genes, cfg, z) {
  out <- data.table(cohort_id = cfg$cohort_id, species = cfg$species,
    n_genes = nrow(genes), n_estimable = 0L, status = "no_estimable_TE_nonTE_pairs")
  if (!"magnitude_difference_log2" %in% names(genes)) return(out)
  d <- genes[status == "eligible_conditional_inference" &
    is.finite(magnitude_difference_log2)]
  out$n_estimable <- nrow(d)
  if (!nrow(d)) return(out)
  out[, `:=`(mean_magnitude_difference_log2 = mean(d$magnitude_difference_log2),
    mean_signed_difference_log2 = mean(d$signed_difference_log2),
    status = if (nrow(d) < 2L) "descriptive_single_gene" else "exploratory_gene_bootstrap")]
  if (nrow(d) >= 2L) {
    set.seed(z$seed)
    boot <- replicate(z$bootstrap, mean(d$magnitude_difference_log2[
      sample.int(nrow(d), nrow(d), replace = TRUE)]))
    ci <- quantile(boot, c(.025, .975))
    out[, `:=`(magnitude_ci_low = unname(ci[1]), magnitude_ci_high = unname(ci[2]))]
  }
  out[, interpretation := "within_cohort_gene_summary_not_independent_cross_species_replication"]
  out
}
reg_te_worker <- function(cfg, z, root) {
  reference_gate(cfg, root)
  if (!stage_valid(cfg$paths$model_receipt,
      unlist(cfg$paths[c("model_gene", "model_promoter", "model_pairs")])))
    stop("Complete frozen primary model outputs required")
  m <- load_manifests(cfg); counts <- completed_counts(cfg, root, m)
  if (!reg_valid(z, "regulatory-resources")) stop("Run regulatory-resources first")
  repeat_path <- if (cfg$species == "fish") reg_fish_paths(z)$intervals else
    file.path(z$output_dir, "repeat_intervals.tsv.gz")
  if (cfg$species == "fish") {
    if (!reg_valid(z, "fish-masking")) stop("Run fish-repeats discovery and masking first")
  }
  names <- c("transcript_TE_annotation.tsv", "event_TE_annotation.tsv",
    "TE_annotation_coverage.tsv", "te_tss_event_pairs.tsv", "te_tss_gene_summary.tsv",
    "te_tss_cohort_summary.tsv", "te_tss_provenance.tsv")
  outputs <- file.path(z$output_dir, names)
  if (!reg_begin(z, "te-tss", outputs)) return(invisible(TRUE))
  seq <- fread(file = cfg$paths$genome_fai, header = FALSE, select = 1:2,
    col.names = c("chrom", "length"))
  tx <- reg_read(cfg$paths$tx_map); map <- reg_read(cfg$paths$event_map)
  r <- reg_read(repeat_path)
  aliases <- if (!is.null(z$contig_aliases)) reg_read(z$contig_aliases) else NULL
  resources <- reg_read(file.path(z$output_dir, "resource_manifest.tsv"))
  full <- cfg$species == "fish" ||
    identical(resources[resource_id == "repeats"]$scope, "full_reference")
  anno <- reg_annotate_events(tx, map, r, seq, aliases, full)
  pairs <- reg_te_pairs(reg_read(cfg$paths$events), anno$events)
  result <- reg_analyze_pairs(pairs, counts$event, counts$gene, m$units,
    if (cfg$profile == "human_trio") reg_read(cfg$paths$covariates) else NULL, cfg, z)
  genes <- if (nrow(result)) result[!duplicated(gene_id)] else copy(result)
  genes[, `:=`(signed_gene_q_value = NA_real_, magnitude_gene_q_value = NA_real_)]
  genes[status == "eligible_conditional_inference",
    `:=`(signed_gene_q_value = p.adjust(signed_p_value, "BH"),
      magnitude_gene_q_value = p.adjust(magnitude_p_value, "BH"))]
  summaries <- reg_cohort_summary(genes, cfg, z)
  objects <- list(anno$transcripts, anno$events, anno$coverage, result, genes, summaries,
    data.table(cohort_id = cfg$cohort_id, build = z$build,
      receipt_policy = "immutable_inputs_output_presence", version = z$version, bootstrap = z$bootstrap, seed = z$seed,
      primary_results_modified = FALSE, same_first_donor = TRUE,
      magnitude_caveat = "absolute_effects_and_bootstrap_nonregular_near_zero",
      model_scope = "conditional_on_annotation_quantifier_and_frozen_locus"))
  for (i in seq_along(outputs)) reg_write(objects[[i]], outputs[i])
  reg_finish(z, "te-tss", outputs)
}
