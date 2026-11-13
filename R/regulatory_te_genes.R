# Hypothesis #1, between-gene complement (te-genes): do genes with at least one
# TE-derived TSS show a different gene-level parent-of-origin effect (POE) than
# comparable genes whose TSSs are all non-TE?
#
# Exposure  : gene has >=1 multi-exon transcript TSS inside a TE ("TE"), or no
#             TSS inside a TE or unknown repeat ("nonTE"); others are excluded.
# Outcome   : gene POE estimate b and its SE. Magnitude uses b^2 - SE^2, which is
#             unbiased for the true squared effect (noise inflates |b|); signed
#             uses b. Public crosses: primary model_gene (log2, dosage-adjusted);
#             GUSTO: additional_v1 gene_POE (ln, adjusted for cis and covariates).
# Confounders: depth, informative loci, units, SE, number of TSS events and
#             transcripts, gene length, local TE density, promoter GC.
# Design    : propensity-score subclassification (ATT weights) on common support,
#             balance table, 2 Mb genomic block bootstrap, within-stratum label
#             permutation, and a placebo that circularly shifts the TE track per
#             chromosome (keeps TE clustering, breaks the TSS-TE link).
# This is an exploratory association within one cohort, not a causal estimate.
TE_GENES_FILES <- c("te_genes_gene_table.tsv", "te_genes_balance.tsv",
  "te_genes_result.tsv", "te_genes_placebo.tsv", "te_genes_provenance.tsv")
TE_GENES_COVARIATES <- c("log_depth", "log_informative_loci", "n_units", "log_se",
  "log_tss_events", "log_transcripts", "log_gene_length", "te_density", "promoter_gc")

reg_te_genes_settings <- function(z) {
  s <- reg_default(z$te_genes, list())
  out <- list(min_genes = reg_default(s$min_genes_per_group, 5),
    strata = reg_default(s$propensity_strata, 10),
    block_bp = reg_default(s$block_bp, 2e6),
    flank_bp = reg_default(s$te_density_flank_bp, 1e5),
    gc_flank_bp = reg_default(s$gc_flank_bp, 1000),
    permutations = reg_default(s$permutations, 2000),
    placebo = reg_default(s$placebo_shifts, 200),
    min_shift_bp = reg_default(s$min_shift_bp, 1e6))
  for (k in names(out)) positive(out[[k]], paste0("regulatory.te_genes.", k))
  for (k in c("min_genes", "strata", "permutations", "placebo"))
    out[[k]] <- as.integer(out[[k]])
  out
}

# Merge overlapping/adjacent intervals per contig.
reg_union <- function(x) {
  if (!nrow(x)) return(data.table(chrom = character(), start = numeric(), end = numeric()))
  x <- x[, .(chrom = as.character(chrom), start = as.numeric(start), end = as.numeric(end))]
  setorder(x, chrom, start, end)
  x[, grp := cumsum(c(TRUE, start[-1L] > cummax(end)[-.N] + 1)), by = chrom]
  x <- x[, .(start = min(start), end = max(end)), by = .(chrom, grp)]
  x[, grp := NULL]
  x[]
}

# TRUE where a point lies inside the (merged) intervals of u.
reg_points_in <- function(chrom, pos, u) {
  out <- logical(length(pos))
  for (ch in intersect(unique(chrom), unique(u$chrom))) {
    i <- which(chrom == ch)
    v <- u[u$chrom == ch]
    k <- findInterval(pos[i], v$start)
    out[i] <- k > 0L & pos[i] <= v$end[pmax(k, 1L)]
  }
  out
}

# Base pairs of each window covered by the merged intervals u.
reg_window_coverage <- function(w, u) {
  cov <- numeric(nrow(w))
  if (!nrow(u) || !nrow(w)) return(cov)
  x <- data.table(chrom = as.character(w$chrom), start = as.numeric(w$start),
    end = as.numeric(w$end), row = seq_len(nrow(w)))
  y <- copy(u); setkey(y, chrom, start, end)
  h <- foverlaps(x, y, nomatch = 0L)
  if (nrow(h)) {
    s <- h[, .(bp = sum(pmin(end, i.end) - pmax(start, i.start) + 1)), by = row]
    cov[s$row] <- s$bp
  }
  cov
}

reg_gene_exposure <- function(tss, te_u, rep_u) {
  d <- data.table(gene_id = tss$gene_id,
    in_te = reg_points_in(tss$chrom, tss$tss, te_u),
    in_rep = reg_points_in(tss$chrom, tss$tss, rep_u))
  d[, .(exposure = if (any(in_te)) "TE" else if (any(in_rep)) "unknown_repeat" else "nonTE"),
    by = gene_id]
}

reg_promoter_gc <- function(genome, genes, contigs, flank) {
  len <- contigs$length[match(genes$chrom, contigs$chrom)]
  st <- pmax(1, genes$tss5 - flank)
  en <- pmin(len, genes$tss5 + flank)
  ok <- is.finite(len) & is.finite(st) & is.finite(en) & en > st
  gc <- rep(NA_real_, nrow(genes))
  if (any(ok)) {
    gr <- GenomicRanges::GRanges(genes$chrom[ok], IRanges::IRanges(st[ok], en[ok]))
    s <- Rsamtools::scanFa(Rsamtools::FaFile(genome), gr)
    f <- Biostrings::letterFrequency(s, c("A", "C", "G", "T"))
    gc[ok] <- (f[, "G"] + f[, "C"]) / pmax(1, rowSums(f))
  }
  gc
}

# Circular shift of merged intervals by one offset per contig. A rotation keeps
# merged intervals disjoint, so shifting the merged track equals merging the
# shifted track; intervals that would wrap past the contig end are dropped.
reg_shift_offsets <- function(contigs, min_shift) {
  stats::setNames(floor(vapply(contigs$length, function(L)
    if (L > 2 * min_shift) stats::runif(1, min_shift, L - min_shift) else stats::runif(1, 0, L), 0)),
    contigs$chrom)
}
reg_shift_union <- function(u, contigs, offs) {
  if (!nrow(u)) return(u)
  len <- contigs$length[match(u$chrom, contigs$chrom)]
  start <- ((u$start - 1 + offs[u$chrom]) %% len) + 1
  out <- data.table(chrom = u$chrom, start = start, end = start + (u$end - u$start))
  out <- out[out$end <= len]
  setorder(out, chrom, start)
  out[]
}

# Propensity-score subclassification on one analysis set (TE/nonTE genes with
# complete covariates). Returns NULL when either group is too small.
reg_te_genes_estimate <- function(d, covs, k, min_genes) {
  y <- d$exposure == "TE"
  if (min(sum(y), sum(!y)) < min_genes) return(NULL)
  X <- as.matrix(d[, covs, with = FALSE])
  sdv <- apply(X, 2, stats::sd)
  X <- X[, is.finite(sdv) & sdv > 0, drop = FALSE]
  if (!ncol(X)) return(NULL)
  X <- scale(X)
  fit <- tryCatch(suppressWarnings(stats::glm.fit(cbind(1, X), as.numeric(y),
    family = stats::binomial())), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ps <- fit$fitted.values
  lo <- max(min(ps[y]), min(ps[!y])); hi <- min(max(ps[y]), max(ps[!y]))
  inside <- ps >= lo & ps <= hi
  if (min(sum(y & inside), sum(!y & inside)) < min_genes) return(NULL)
  br <- unique(stats::quantile(ps[y & inside], probs = seq(0, 1, length.out = k + 1L), names = FALSE))
  stratum <- rep(NA_integer_, length(ps))
  stratum[inside] <- if (length(br) < 2L) 1L else
    findInterval(ps[inside], br, rightmost.closed = TRUE, all.inside = TRUE)
  n1 <- tapply(y[inside], stratum[inside], sum)
  n0 <- tapply(!y[inside], stratum[inside], sum)
  valid <- names(n1)[n1 > 0 & n0 > 0]
  stratum[!(as.character(stratum) %in% valid)] <- NA_integer_
  ratio <- stats::setNames(as.numeric(n1[valid] / n0[valid]), valid)
  w <- reg_te_genes_weights(y, stratum, ratio)
  if (min(sum(y & w > 0), sum(!y & w > 0)) < min_genes) return(NULL)
  diffs <- reg_te_genes_diff(d, y, w)
  list(magnitude = diffs[["magnitude"]], signed = diffs[["signed"]], y = y, ps = ps,
    stratum = stratum, ratio = ratio, weight = w, n_te = sum(y & w > 0),
    n_nonte = sum(!y & w > 0), n_trimmed = sum(!inside), n_strata = length(valid))
}
reg_te_genes_weights <- function(y, stratum, ratio) {
  w <- ifelse(is.na(stratum), 0, ifelse(y, 1, ratio[as.character(stratum)]))
  w[is.na(w)] <- 0
  w
}
reg_te_genes_diff <- function(d, y, w) {
  t <- y & w > 0; c0 <- !y & w > 0
  wm <- function(v, i) sum(w[i] * v[i]) / sum(w[i])
  c(magnitude = wm(d$magnitude_unbiased, t) - wm(d$magnitude_unbiased, c0),
    signed = wm(d$estimate, t) - wm(d$estimate, c0))
}

reg_te_genes_balance <- function(d, covs, w) {
  y <- d$exposure == "TE"
  rbindlist(lapply(covs, function(v) {
    x <- d[[v]]
    sdp <- sqrt((stats::var(x[y]) + stats::var(x[!y])) / 2)
    t <- y & w > 0; c0 <- !y & w > 0
    m1w <- sum(w[t] * x[t]) / sum(w[t]); m0w <- sum(w[c0] * x[c0]) / sum(w[c0])
    data.table(covariate = v, mean_TE = mean(x[y]), mean_nonTE = mean(x[!y]),
      smd_before = (mean(x[y]) - mean(x[!y])) / sdp,
      mean_TE_weighted = m1w, mean_nonTE_weighted = m0w, smd_after = (m1w - m0w) / sdp)
  }))
}

reg_te_genes_inputs <- function(cfg, z, root, s, baseline = NULL) {
  contigs <- fread(file = cfg$paths$genome_fai, header = FALSE, select = 1:2,
    col.names = c("chrom", "length"))
  contigs[, chrom := as.character(chrom)]
  tx <- reg_read(cfg$paths$tx_map)
  require_columns(tx, c("transcript_id", "gene_id", "chrom", "strand", "start", "end", "tss"),
    "Transcript map")
  tx[, chrom := as.character(chrom)]
  em <- reg_read(cfg$paths$event_map)
  require_columns(em, c("transcript_id", "gene_id", "event_id", "window"), "Event map")
  em50 <- em[em$window == 50L]
  tss <- tx[tx$transcript_id %in% em50$transcript_id, .(gene_id, chrom, tss)]
  genes <- tx[, .(n_chrom = uniqueN(chrom), chrom = chrom[1], strand = strand[1],
    start = min(start), end = max(end), n_transcripts = .N,
    tss5 = if (strand[1] == "+") min(tss) else max(tss)), by = gene_id]
  genes <- genes[n_chrom == 1L][, n_chrom := NULL]
  ev <- unique(em50[, .(gene_id, event_id)])[, .(n_tss_events = .N), by = gene_id]
  genes <- merge(genes, ev, by = "gene_id", all.x = TRUE)
  genes[is.na(n_tss_events), n_tss_events := 0L]

  repeat_path <- if (identical(cfg$species, "fish")) reg_fish_paths(z)$intervals else
    file.path(z$output_dir, "repeat_intervals.tsv.gz")
  r <- reg_read(repeat_path)
  aliases <- if (!is.null(z$contig_aliases)) reg_read(z$contig_aliases) else NULL
  r[, chrom := reg_alias(chrom, contigs$chrom, aliases)]
  r <- r[!is.na(r$chrom), .(chrom, start = as.numeric(start), end = as.numeric(end),
    is_te = is_te %in% TRUE, is_unknown = is_unknown %in% TRUE)]
  te_u <- reg_union(r[r$is_te])

  if (!is.null(baseline)) {
    bs <- baseline$s
    nb_require(bs, "combine")
    poe <- nb_estimable(reg_read(nb_path(bs, "gene_POE.tsv")))
    poe[, `:=`(gene_id = feature_id, estimate = beta_adjusted, n_units = n_samples)]
    scale <- "ln_paternal_maternal_dosage_adjusted_all_SNP_NB"
    source_table <- "baseline_nb_v1/gene_POE.tsv"
  } else if (cfg$profile == "human_trio") {
    if (!stage_valid(additional_receipt(cfg), additional_paths(cfg)))
      stop("te-genes for human cohorts needs additional_v1/gene_POE.tsv: run the additional stage first")
    poe <- reg_read(additional_paths(cfg)[["genes"]])
    require_columns(poe, c("gene_id", "estimate", "se", "status", "n_families"), "Gene POE")
    poe[, n_units := n_families]
    scale <- "ln_maternal_paternal_adjusted_for_cis_and_covariates"
    source_table <- "additional_v1/gene_POE.tsv"
  } else {
    poe <- reg_read(cfg$paths$model_gene)
    require_columns(poe, c("entity_id", "estimate", "se", "status", "n_units"), "Gene POE")
    setnames(poe, "entity_id", "gene_id")
    scale <- "log2_maternal_paternal_dosage_adjusted"
    source_table <- "model_gene.tsv"
  }
  poe_flags <- if (!is.null(baseline)) poe[, .(gene_id, not_strictly_POE, independence_verified)] else NULL
  poe <- poe[(if (is.null(baseline)) poe$status == "eligible_conditional_inference" else
      poe$status %in% c("eligible_NB", "exploratory_NB")) & is.finite(poe$estimate) &
    is.finite(poe$se) & poe$se > 0, .(gene_id, estimate, se, n_units = as.numeric(n_units))]
  assert_unique(poe, "gene_id", "Gene POE")
  if (!is.null(poe_flags)) poe <- merge(poe, poe_flags, by = "gene_id", all.x = TRUE)

  if (is.null(baseline)) {
    counts <- completed_counts(cfg, root, load_manifests(cfg))$gene
    require_columns(counts, c("gene_id", "locus_id", "maternal_count", "paternal_count"), "Gene counts")
    counts[, depth := maternal_count + paternal_count]
    dep <- counts[is.finite(depth), .(total_depth = sum(depth),
      n_informative_loci = uniqueN(locus_id[depth >= z$min_depth])), by = gene_id]
  } else {
    dep <- reg_read(nb_path(baseline$s, "gene_inventory.tsv.gz"))[, .(
      total_depth = sum(total_counts), n_informative_loci = sum(n_samples_observed > 0)),
      by = .(gene_id = feature_id)]
  }

  len <- contigs$length[match(genes$chrom, contigs$chrom)]
  win <- data.table(chrom = genes$chrom, start = pmax(1, genes$start - s$flank_bp),
    end = pmin(len, genes$end + s$flank_bp))
  genes[, te_density := reg_window_coverage(win, te_u) / pmax(1, win$end - win$start + 1)]
  genes[, promoter_gc := reg_promoter_gc(cfg$paths$genome, genes, contigs, s$gc_flank_bp)]

  g <- merge(genes, poe, by = "gene_id")
  g <- merge(g, dep, by = "gene_id", all.x = TRUE)
  g[is.na(total_depth), `:=`(total_depth = 0, n_informative_loci = 0L)]
  g[, `:=`(log_depth = log10(total_depth + 1), log_informative_loci = log1p(n_informative_loci),
    log_se = log(se), log_tss_events = log1p(n_tss_events), log_transcripts = log1p(n_transcripts),
    log_gene_length = log10(end - start + 1), magnitude_unbiased = estimate^2 - se^2,
    block = paste(chrom, floor(((start + end) / 2) / s$block_bp), sep = ":"))]
  complete <- stats::complete.cases(g[, TE_GENES_COVARIATES, with = FALSE]) &
    is.finite(rowSums(as.matrix(g[, TE_GENES_COVARIATES, with = FALSE])))
  list(genes = g[complete], dropped_incomplete = sum(!complete), tss = tss, repeats = r,
    te_u = te_u, contigs = contigs, scale = scale, source_table = source_table)
}

reg_te_genes_worker <- function(cfg, z, root) {
  if (!reg_valid(z, "te-tss")) stop("Complete te-tss first (repeat track and TE labels)")
  s <- reg_te_genes_settings(z)
  outputs <- file.path(z$output_dir, TE_GENES_FILES)
  if (!reg_begin(z, "te-genes", outputs)) return(invisible(TRUE))
  set.seed(z$seed)
  inp <- reg_te_genes_inputs(cfg, z, root, s)
  reg_te_genes_run(cfg, z, inp, s, outputs)
  reg_finish(z, "te-genes", outputs)
}
reg_te_genes_run <- function(cfg, z, inp, s, outputs) {
  set.seed(z$seed)
  covs <- TE_GENES_COVARIATES
  rep_u <- reg_union(inp$repeats[inp$repeats$is_te | inp$repeats$is_unknown])
  exposure <- reg_gene_exposure(inp$tss, inp$te_u, rep_u)
  g <- merge(inp$genes, exposure, by = "gene_id", all.x = TRUE)
  g[is.na(exposure), exposure := "no_multi_exon_TSS_event"]
  A <- g[exposure %in% c("TE", "nonTE")]
  est <- reg_te_genes_estimate(A, covs, s$strata, s$min_genes)

  result <- data.table(cohort_id = cfg$cohort_id, species = cfg$species,
    effect_scale = inp$scale, gene_POE_source = inp$source_table,
    n_genes_with_POE = nrow(g), n_TE_genes = sum(A$exposure == "TE"),
    n_nonTE_genes = sum(A$exposure == "nonTE"),
    n_unknown_repeat_genes_excluded = sum(g$exposure == "unknown_repeat"),
    n_incomplete_covariates_excluded = inp$dropped_incomplete,
    status = "insufficient_genes_per_group")
  if ("not_strictly_POE" %in% names(g))
    result[, `:=`(n_not_strictly_POE = sum(g$not_strictly_POE %in% TRUE),
      n_independence_unverified = sum(!g$independence_verified),
      cross_flag_policy = "retained_main_analysis_separate_sensitivity")]
  balance <- data.table(covariate = character(), smd_before = numeric(), smd_after = numeric())
  placebo <- data.table(draw = integer(), n_TE_genes = integer(),
    magnitude_difference = numeric(), signed_difference = numeric())
  A[, `:=`(propensity = NA_real_, stratum = NA_integer_, weight = 0)]

  if (!is.null(est)) {
    A[, `:=`(propensity = est$ps, stratum = est$stratum, weight = est$weight)]
    balance <- reg_te_genes_balance(A, covs, est$weight)
    # Block bootstrap: resample 2 Mb genomic blocks, refit the propensity model.
    blocks <- split(seq_len(nrow(A)), A$block)
    boot <- t(vapply(seq_len(z$bootstrap), function(b) {
      pick <- unlist(blocks[sample.int(length(blocks), length(blocks), replace = TRUE)],
        use.names = FALSE)
      e <- reg_te_genes_estimate(A[pick], covs, s$strata, s$min_genes)
      if (is.null(e)) c(NA_real_, NA_real_) else c(e$magnitude, e$signed)
    }, numeric(2)))
    # Permutation: shuffle TE labels within propensity strata (fixed weights rule).
    keep <- which(!is.na(est$stratum))
    perm <- t(vapply(seq_len(s$permutations), function(p) {
      y <- est$y
      for (idx in split(keep, est$stratum[keep]))
        y[idx] <- if (length(idx) > 1L) y[sample(idx)] else y[idx]
      reg_te_genes_diff(A, y, reg_te_genes_weights(y, est$stratum, est$ratio))
    }, numeric(2)))
    # Placebo: circularly shift the TE track per chromosome; recompute exposure.
    base <- g[, !"exposure"]
    placebo <- rbindlist(lapply(seq_len(s$placebo), function(i) {
      offs <- reg_shift_offsets(inp$contigs, s$min_shift_bp)
      ex <- reg_gene_exposure(inp$tss, reg_shift_union(inp$te_u, inp$contigs, offs),
        reg_shift_union(rep_u, inp$contigs, offs))
      Ai <- merge(base, ex, by = "gene_id")[exposure %in% c("TE", "nonTE")]
      e <- reg_te_genes_estimate(Ai, covs, s$strata, s$min_genes)
      data.table(draw = i, n_TE_genes = sum(Ai$exposure == "TE"),
        magnitude_difference = if (is.null(e)) NA_real_ else e$magnitude,
        signed_difference = if (is.null(e)) NA_real_ else e$signed)
    }))
    observed <- c(magnitude = est$magnitude, signed = est$signed)
    for (j in 1:2) {
      nm <- names(observed)[j]; b <- boot[, j]; ok <- is.finite(b)
      se_b <- if (sum(ok) > 1L) stats::sd(b[ok]) else NA_real_
      ci <- if (sum(ok) > 1L) stats::quantile(b[ok], c(.025, .975), names = FALSE) else c(NA_real_, NA_real_)
      pl <- placebo[[paste0(nm, "_difference")]]; pl <- pl[is.finite(pl)]
      result[, (paste0(nm, c("_difference", "_bootstrap_se", "_ci_low", "_ci_high",
        "_p_bootstrap", "_p_permutation", "_p_placebo"))) := list(observed[[nm]], se_b,
        ci[1], ci[2], if (is.finite(se_b) && se_b > 0) 2 * stats::pnorm(-abs(observed[[nm]] / se_b)) else NA_real_,
        (1 + sum(abs(perm[, j]) >= abs(observed[[nm]]))) / (1 + nrow(perm)),
        if (length(pl)) (1 + sum(abs(pl) >= abs(observed[[nm]]))) / (1 + length(pl)) else NA_real_)]
    }
    max_smd <- max(abs(balance$smd_after), na.rm = TRUE)
    result[, `:=`(n_TE_genes_matched = est$n_te, n_nonTE_genes_matched = est$n_nonte,
      n_trimmed_off_support = est$n_trimmed, n_propensity_strata = est$n_strata,
      n_genomic_blocks = length(blocks),
      bootstrap_valid_fraction = mean(is.finite(boot[, 1])),
      n_placebo_valid = sum(is.finite(placebo$magnitude_difference)),
      max_abs_smd_after = max_smd,
      status = if (max_smd > 0.1) "exploratory_residual_imbalance_gt_0.1" else
        "exploratory_propensity_subclassification",
      interpretation = paste("within_cohort_between_gene_association_not_causal;",
        "magnitude_is_difference_in_unbiased_squared_POE"))]
  }

  provenance <- data.table(setting = c("version", "seed", "bootstrap", "min_depth",
      names(s), "covariates", "repeat_source", "gene_POE_source"),
    value = c(z$version, z$seed, z$bootstrap, z$min_depth, unlist(lapply(s, as.character)),
      paste(covs, collapse = ";"),
      if (identical(cfg$species, "fish")) "fish_de_novo_RepeatMasker" else "resource_manifest repeats",
      inp$source_table))
  out_cols <- c("gene_id", "chrom", "start", "end", "exposure", "estimate", "se",
    "magnitude_unbiased", covs, "propensity", "stratum", "weight", "block")
  out_cols <- c(out_cols, intersect(c("not_strictly_POE", "independence_verified"), names(A)))
  reg_write(A[, out_cols, with = FALSE], outputs[1])
  reg_write(balance, outputs[2])
  reg_write(result, outputs[3])
  reg_write(placebo, outputs[4])
  reg_write(provenance, outputs[5])
}
