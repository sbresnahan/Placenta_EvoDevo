ex_failure <- function(status, details = list()) list(row = data.table::as.data.table(c(
  list(status = status, version = EX_VERSION), details)), fit = NULL)
ex_long <- function(z, metadata) {
  ex_need(z, c("feature_id", "SNPpos", "sample_id", "maternal_count", "paternal_count", "offset"))
  ex_unique(z, c("feature_id", "SNPpos", "sample_id"), "SNP count input")
  if (data.table::uniqueN(z$feature_id) != 1L) stop("Fit one feature at a time")
  if (any(!grepl("^[^:]+:[1-9][0-9]*:[ACGT]:[ACGT]$", z$SNPpos))) stop("Full genomic SNP keys required")
  for (nm in c("maternal_count", "paternal_count"))
    if (!ex_numeric(z[[nm]], nonnegative = TRUE, allow_na = TRUE)) stop("Invalid allele counts")
  if (!ex_numeric(z$offset)) stop("Invalid count offset")
  # Metadata must agree with the baseline preparation where both declare values.
  for (nm in intersect(c("cluster_id", "sex", "ancestry", "cross"), names(z))) {
    if (!nm %in% names(metadata)) next
    ref <- metadata[[nm]][match(z$sample_id, metadata$sample_id)]
    bad <- !is.na(ref) & !is.na(z[[nm]]) & as.character(ref) != as.character(z[[nm]])
    if (any(bad)) stop("Crosswalk disagrees with baseline ", nm)
  }
  meta <- data.table::copy(metadata); ex_unique(meta, "sample_id", "fit metadata")
  x <- data.table::copy(z[, .(feature_id, SNPpos, sample_id, maternal_count, paternal_count, offset)])
  x <- merge(x, meta, by = "sample_id", all = FALSE, sort = FALSE)
  x <- x[is.finite(maternal_count) & is.finite(paternal_count)]
  x <- data.table::rbindlist(list(
    x[, `:=`(allele = "maternal", count = round(maternal_count))][],
    data.table::copy(x)[, `:=`(allele = "paternal", count = round(paternal_count))][]))
  x[, c("maternal_count", "paternal_count") := NULL]
  x[, `:=`(allele = factor(allele, c("maternal", "paternal")),
    SNPpos = factor(SNPpos), sample_id = factor(sample_id))]
  x
}
ex_sample_design <- function(x, public, mediator = FALSE) {
  vars <- c(if (public) "cross", "E", if (mediator) "M",
    intersect(c("sex", "ancestry"), names(x)))
  one <- unique(x[, c("sample_id", "cluster_id", vars), with = FALSE])
  ex_unique(one, "sample_id", "sample covariates")
  fixed <- stats::as.formula(paste("~", paste(vars, collapse = " + ")))
  X <- stats::model.matrix(fixed, one)
  list(one = one, X = X, rank = qr(X)$rank, kappa = kappa(X, exact = TRUE))
}
ex_nb <- function(z, metadata, cfg, s, mode = c("h6", "total", "mediation", "modification")) {
  mode <- match.arg(mode); public <- identical(cfg$profile, "public_cross")
  if (mode != "h6") ex_require_human(cfg)
  ex_need(metadata, c("sample_id", "cluster_id", "independence_verified", "E"))
  if (mode %in% c("mediation", "modification")) ex_need(metadata, "M")
  keep <- is.finite(metadata$E)
  if ("M" %in% names(metadata) && mode %in% c("mediation", "modification")) keep <- keep & is.finite(metadata$M)
  meta <- data.table::copy(metadata[keep])
  x <- ex_long(z, meta)
  if (!nrow(x)) return(ex_failure("no_complete_informative_samples"))
  positive <- x[, .(depth = sum(count)), by = sample_id][depth > 0]$sample_id
  x <- x[sample_id %in% positive]
  nu <- data.table::uniqueN(x$cluster_id); ns <- data.table::uniqueN(x$sample_id)
  n_snp <- data.table::uniqueN(as.character(x[count > 0]$SNPpos))
  info <- list(n_samples = ns, n_biological_units = nu, n_snps = n_snp, total_counts = sum(x$count),
    mode = mode, small_cohort = public, independence_verified = all(x$independence_verified))
  fail <- function(status) ex_failure(status, info)
  min_units <- if (public) 4L else 15L
  if (ns < min_units || nu < min_units) return(fail("insufficient_biological_units"))
  if (n_snp < 3L || sum(x$count) < 50) return(fail("insufficient_SNPs_or_counts"))
  if (!all(x$independence_verified) && !isTRUE(s$allow_unverified_independence))
    return(fail("unverified_biological_independence"))
  covars <- character()
  for (nm in c("sex", "ancestry")) {
    if (!nm %in% names(x)) {if (!public) return(fail(paste0("missing_", nm))); next}
    if (anyNA(x[[nm]])) {if (!public || any(!is.na(x[[nm]]))) return(fail(paste0("missing_", nm))); next}
    if (data.table::uniqueN(x[[nm]]) > 1L) {
      x[, (nm) := factor(get(nm))]; covars <- c(covars, nm)
    } else if (!public) return(fail(paste0("unidentifiable_", nm)))
  }
  drop <- setdiff(c("sex", "ancestry"), covars)
  if (length(intersect(drop, names(x)))) x[, (intersect(drop, names(x))) := NULL]
  within_centers <- NULL
  if (public) {
    ex_need(x, "cross")
    units <- unique(x[, .(cluster_id, cross)])
    if (anyDuplicated(units$cluster_id)) return(fail("cluster_spans_cross_directions"))
    n_cross <- units[, .N, by = cross]
    if (nrow(n_cross) != 2L || anyNA(n_cross$cross) || min(n_cross$N) < 2L)
      return(fail("insufficient_reciprocal_units"))
    x[, cross := factor(cross, c("focal_mother", "alternate_mother"))]
    sample_e <- unique(x[, .(sample_id, cross, E)])
    v <- sample_e[, .(sd = stats::sd(E), center = mean(E)), by = cross]
    if (any(!is.finite(v$sd) | v$sd <= 0)) return(fail("no_within_cross_exposure_variation"))
    within_centers <- v
    x[, E_original := E]
    x[, E := E - v$center[match(cross, v$cross)]]
  }
  sample_e <- unique(x[, .(sample_id, E)])
  if (data.table::uniqueN(sample_e$E) < 2L || !is.finite(stats::sd(sample_e$E)) ||
      stats::sd(sample_e$E) <= 0) return(fail("constant_or_insufficient_exposure_values"))
  has_m <- mode %in% c("mediation", "modification")
  if (has_m && (data.table::uniqueN(x$M) < 3L || stats::sd(unique(x[, .(sample_id, M)])$M) <= 0))
    return(fail("constant_or_insufficient_mediator_values"))
  sdsgn <- ex_sample_design(x, public, has_m)
  if (mode == "modification") {
    Xs <- stats::model.matrix(stats::as.formula(paste("~ E * M",
      if (length(covars)) paste("+", paste(covars, collapse = " + ")) else "")), sdsgn$one)
    sdsgn$X <- Xs; sdsgn$rank <- qr(Xs)$rank; sdsgn$kappa <- kappa(Xs, exact = TRUE)
  }
  info$sample_design_rank <- sdsgn$rank; info$residual_unit_df <- nu - sdsgn$rank
  info$condition_number <- sdsgn$kappa
  if (sdsgn$rank < ncol(sdsgn$X)) return(fail("rank_deficient_sample_design"))
  if (nu <= sdsgn$rank) return(fail("no_residual_biological_unit_support"))
  if (!is.finite(sdsgn$kappa) || sdsgn$kappa > s$max_condition) return(fail("ill_conditioned_design"))
  if (mode == "modification") {
    one <- sdsgn$one
    quad <- table(factor(one$E > stats::median(one$E), levels = c(FALSE, TRUE)),
                  factor(one$M > stats::median(one$M), levels = c(FALSE, TRUE)))
    info$min_joint_quadrant <- min(quad)
    if (nu < 30L || min(quad) < 5L) return(fail("insufficient_joint_exposure_support"))
  }
  terms <- if (mode == "modification") "allele * E * M" else
    if (has_m) "allele * (E + M)" else
      if (public) "allele * cross + allele * E" else "allele * E"
  rhs <- paste(c(terms, covars), collapse = " + ")
  form <- stats::as.formula(paste("count ~", rhs, "+ (1|SNPpos) + (1|sample_id)"))
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), x)
  if (qr(X)$rank < ncol(X)) return(fail("rank_deficient_count_design"))
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(glmmTMB::glmmTMB(form, data = x,
    family = glmmTMB::nbinom2(), offset = offset,
    control = glmmTMB::glmmTMBControl(optimizer = nlminb, parallel = 1L,
      optCtrl = list(rel.tol = 1e-5, eval.max = 300L, iter.max = 300L))),
    warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
    error = function(e) e)
  info$formula <- paste(deparse(form), collapse = "")
  info$warnings <- paste(unique(warnings), collapse = " | ")
  if (inherits(fit, "error")) {info$warnings <- conditionMessage(fit); return(fail("fit_error"))}
  info$optimizer_code <- fit$fit$convergence; info$pdHess <- isTRUE(fit$sdr$pdHess)
  if (fit$fit$convergence != 0L || !isTRUE(fit$sdr$pdHess)) return(fail("failed_convergence_or_Hessian"))
  b <- glmmTMB::fixef(fit)$cond; V <- as.matrix(stats::vcov(fit)$cond)
  if (any(!is.finite(b)) || any(!is.finite(V))) return(fail("nonfinite_uncertainty"))
  required <- c("allelepaternal", "allelepaternal:E", if (has_m) "allelepaternal:M",
    if (mode == "modification") "allelepaternal:E:M")
  if (!all(required %in% names(b))) return(fail("missing_required_coefficient"))
  info$E_min <- min(x$E); info$E_max <- max(x$E)
  info$E_p10 <- as.numeric(stats::quantile(sample_e$E, .1)); info$E_p50 <- stats::median(sample_e$E)
  info$E_p90 <- as.numeric(stats::quantile(sample_e$E, .9))
  for (pair in list(c("slope", "allelepaternal:E"), c("b", "allelepaternal:M"),
                    c("modification", "allelepaternal:E:M"))) {
    if (!pair[2] %in% names(b)) next
    info[[pair[1]]] <- b[[pair[2]]]
    info[[paste0(pair[1], "_se")]] <- sqrt(V[pair[2], pair[2]])
    info[[paste0(pair[1], "_p")]] <- 2 * stats::pnorm(-abs(info[[pair[1]]] / info[[paste0(pair[1], "_se")]]))
  }
  info$pfas_TE_correlation <- if (has_m) stats::cor(sdsgn$one$E, sdsgn$one$M) else NA_real_
  rp <- stats::residuals(fit, type = "pearson")
  influence <- data.table::data.table(cluster_id = x$cluster_id, squared = rp^2,
    count = x$count)[, .(residual_mass = sum(squared), count_mass = sum(count)), by = cluster_id]
  info$pearson_RMS <- sqrt(mean(rp^2))
  info$pearson_abs_max <- max(abs(rp))
  info$largest_family_residual_fraction <- max(influence$residual_mass) / sum(influence$residual_mass)
  info$largest_family_count_fraction <- max(influence$count_mass) / sum(influence$count_mass)
  row <- data.table::as.data.table(c(list(status = "eligible", version = EX_VERSION), info))
  result <- list(row = row, fit = fit, coefficients = b, covariance = V, data = x,
    samples = sdsgn$one, public = public, within_centers = within_centers,
    dosage_null = log(cfg$dosage$paternal / cfg$dosage$maternal))
  pred <- ex_prediction_summary(result)
  result$row <- cbind(row, data.table::as.data.table(pred))
  result
}
ex_weights <- function(result, E, M = 0) {
  b <- result$coefficients
  w <- stats::setNames(rep(0, length(b)), names(b)); w["allelepaternal"] <- 1
  values <- c("allelepaternal:E" = E, "allelepaternal:M" = M, "allelepaternal:E:M" = E * M,
    "allelepaternal:crossalternate_mother" = .5)
  nms <- intersect(names(values), names(w)); w[nms] <- values[nms]; w
}
ex_contrast <- function(result, w, null = 0, alpha = .05) {
  est <- sum(w * result$coefficients) - null
  se <- sqrt(max(0, as.numeric(t(w) %*% result$covariance %*% w)))
  c(beta = est, se = se, low = est - stats::qnorm(1 - alpha / 2) * se,
    high = est + stats::qnorm(1 - alpha / 2) * se,
    p = if (se > 0) 2 * stats::pnorm(-abs(est / se)) else NA_real_)
}
ex_abs_range <- function(lo, hi) c(min = if (lo <= 0 && hi >= 0) 0 else min(abs(c(lo, hi))),
                                  max = max(abs(c(lo, hi))))
ex_prediction_summary <- function(result) {
  e <- as.numeric(stats::quantile(result$samples$E, c(.1, .5, .9)))
  # With a mediator, marginal predictions here hold M at its observed median;
  # natural-effect predictions are separately generated by mediation standardization.
  m <- if ("M" %in% names(result$samples)) stats::median(result$samples$M) else 0
  wl <- ex_weights(result, e[1], m); wh <- ex_weights(result, e[3], m)
  low <- ex_contrast(result, wl, result$dosage_null)
  mid <- ex_contrast(result, ex_weights(result, e[2], m), result$dosage_null)
  high <- ex_contrast(result, wh, result$dosage_null)
  delta <- ex_contrast(result, wh - wl)
  # Simultaneous endpoint rectangle induced by a two-contrast confidence ellipse.
  radius <- sqrt(stats::qchisq(.95, df = 2))
  arl <- ex_abs_range(low["beta"] - radius * low["se"], low["beta"] + radius * low["se"])
  arh <- ex_abs_range(high["beta"] - radius * high["se"], high["beta"] + radius * high["se"])
  list(beta_p10 = unname(low["beta"]), beta_p10_low = unname(low["low"]), beta_p10_high = unname(low["high"]),
    beta_p50 = unname(mid["beta"]), beta_p90 = unname(high["beta"]),
    beta_p90_low = unname(high["low"]), beta_p90_high = unname(high["high"]),
    delta_signed = unname(delta["beta"]), delta_signed_low = unname(delta["low"]),
    delta_signed_high = unname(delta["high"]),
    delta_magnitude = unname(abs(high["beta"]) - abs(low["beta"])),
    delta_magnitude_low = unname(arh["min"] - arl["max"]),
    delta_magnitude_high = unname(arh["max"] - arl["min"]),
    direction_reversal_fitted = unname(low["beta"] * high["beta"] < 0),
    direction_reversal_supported = unname((low["high"] < 0 && high["low"] > 0) ||
                                          (high["high"] < 0 && low["low"] > 0)),
    prediction_M = m, prediction_scale = if (result$public) "within_cross_TE_z_units" else "supplied_exposure_units")
}
