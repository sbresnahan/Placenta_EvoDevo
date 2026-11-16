ex_product_interval <- function(a_low, a_high, b_low, b_high) {
  if (any(!is.finite(c(a_low, a_high, b_low, b_high)))) return(c(low = NA_real_, high = NA_real_))
  if (a_low > a_high || b_low > b_high) stop("Reversed coefficient interval")
  r <- range(outer(c(a_low, a_high), c(b_low, b_high), "*"))
  c(low = r[1], high = r[2])
}
ex_mediator <- function(samples) {
  ex_need(samples, c("sample_id", "cluster_id", "E", "M", "sex", "ancestry"))
  ex_unique(samples, "sample_id", "mediator sample table")
  fit <- stats::lm(M ~ E + sex + ancestry, data = samples)
  X <- stats::model.matrix(fit); nu <- data.table::uniqueN(samples$cluster_id)
  if (qr(X)$rank < ncol(X) || nu <= ncol(X)) stop("Unidentifiable mediator model")
  V <- sandwich::vcovCL(fit, cluster = samples$cluster_id, type = "HC3")
  a <- stats::coef(fit)[["E"]]; se <- sqrt(V["E", "E"]); df <- nu - ncol(X)
  if (!all(is.finite(c(a, se))) || se <= 0) stop("Nonfinite mediator uncertainty")
  list(fit = fit, covariance = V, a = a, a_se = se, a_df = df,
    a_p = 2 * stats::pt(-abs(a / se), df = df), residuals = stats::residuals(fit))
}
ex_mediation_scenarios <- function(outcome, mediator, samples, exposure = stats::quantile(samples$E, c(.25, .75))) {
  e0 <- unname(exposure[1]); e1 <- unname(exposure[2])
  if (!is.finite(e0) || !is.finite(e1) || e1 <= e0) stop("Unsupported exposure contrast")
  d0 <- data.table::copy(samples); d1 <- data.table::copy(samples)
  d0[, E := e0]; d1[, E := e1]
  m0 <- as.numeric(stats::predict(mediator$fit, newdata = d0)) + mediator$residuals
  m1 <- as.numeric(stats::predict(mediator$fit, newdata = d1)) + mediator$residuals
  beta <- function(e, m) vapply(m, function(mm)
    sum(ex_weights(outcome, e, mm) * outcome$coefficients) - outcome$dosage_null, numeric(1))
  b00 <- beta(e0, m0); b10 <- beta(e1, m0); b11 <- beta(e1, m1)
  list(scenario_E_low = e0, scenario_E_high = e1,
    beta_scenario_00 = mean(b00), beta_scenario_10 = mean(b10), beta_scenario_11 = mean(b11),
    natural_direct_beta = mean(b10 - b00), natural_indirect_beta = mean(b11 - b10),
    natural_total_beta = mean(b11 - b00),
    natural_direct_magnitude = mean(abs(b10) - abs(b00)),
    natural_indirect_magnitude = mean(abs(b11) - abs(b10)),
    natural_total_magnitude = mean(abs(b11) - abs(b00)),
    fraction_direction_reversal = mean(b00 * b11 < 0),
    mediator_scenario_outside_observed = mean(c(m0, m1) < min(samples$M) | c(m0, m1) > max(samples$M)))
}
ex_mediation <- function(z, metadata, cfg, s, modification = FALSE) {
  ex_require_human(cfg)
  outcome <- ex_nb(z, metadata, cfg, s, if (modification) "modification" else "mediation")
  if (outcome$row$status != "eligible") return(outcome)
  samples <- data.table::copy(outcome$samples)
  med <- tryCatch(ex_mediator(samples), error = function(e) e)
  if (inherits(med, "error")) {
    outcome$row[, `:=`(status = "mediator_model_failed", mediator_warning = conditionMessage(med))]
    return(outcome)
  }
  a <- med$a; b <- outcome$coefficients[["allelepaternal:M"]]
  cp <- outcome$coefficients[["allelepaternal:E"]]
  b_se <- sqrt(outcome$covariance["allelepaternal:M", "allelepaternal:M"])
  # Each coefficient interval has 97.5% coverage; their rectangle has >=95%.
  za <- stats::qt(1 - .05 / 4, med$a_df); zb <- stats::qnorm(1 - .05 / 4)
  indirect_ci <- ex_product_interval(a - za * med$a_se, a + za * med$a_se, b - zb * b_se, b + zb * b_se)
  # Three simultaneous intervals give a conservative interval for cp + a*b.
  a3 <- stats::qt(1 - .05 / 6, med$a_df); z3 <- stats::qnorm(1 - .05 / 6)
  prod3 <- ex_product_interval(a - a3 * med$a_se, a + a3 * med$a_se, b - z3 * b_se, b + z3 * b_se)
  cp_se <- sqrt(outcome$covariance["allelepaternal:E", "allelepaternal:E"])
  vals <- list(a = a, a_se = med$a_se, a_p = med$a_p, a_df = med$a_df,
    mediator_residual_RMS = sqrt(mean(med$residuals^2)),
    mediator_max_leverage = max(stats::hatvalues(med$fit)),
    mediator_max_Cook_distance = max(stats::cooks.distance(med$fit)),
    c_prime = cp, c_prime_se = cp_se, c_prime_p = outcome$row$slope_p,
    indirect_beta = if (!modification) a * b else NA_real_,
    indirect_low = if (!modification) unname(indirect_ci["low"]) else NA_real_,
    indirect_high = if (!modification) unname(indirect_ci["high"]) else NA_real_,
    indirect_path_p = if (!modification) max(med$a_p, outcome$row$b_p) else NA_real_,
    implied_total_beta = if (!modification) cp + a * b else NA_real_,
    implied_total_low = if (!modification) cp - z3 * cp_se + unname(prod3["low"]) else NA_real_,
    implied_total_high = if (!modification) cp + z3 * cp_se + unname(prod3["high"]) else NA_real_,
    interval_method = "Bonferroni_simultaneous_coefficient_rectangle",
    confounding_assumption = "none_remaining_given_sex_ancestry",
    temporal_ordering = "assumed_PFAS_then_TE_then_POE")
  outcome$row <- cbind(outcome$row, data.table::as.data.table(vals),
    data.table::as.data.table(ex_mediation_scenarios(outcome, med, samples)))
  outcome$mediator <- med
  outcome
}
ex_bootstrap_pair <- function(z, metadata, cfg, s, B = s$bootstrap) {
  # B can be smaller only for explicitly labeled software fixtures.
  original <- ex_mediation(z, metadata, cfg, s)
  if (original$row$status != "eligible") return(list(summary = original$row, draws = data.table::data.table()))
  ids <- unique(original$samples$cluster_id)
  active <- as.character(original$samples$sample_id)
  z <- z[sample_id %in% active]; metadata <- metadata[sample_id %in% active]
  metrics <- c("indirect_beta", "c_prime", "implied_total_beta", "natural_indirect_beta",
    "natural_direct_beta", "natural_total_beta", "natural_indirect_magnitude",
    "natural_direct_magnitude", "natural_total_magnitude")
  draws <- data.table::rbindlist(lapply(seq_len(B), function(i) {
    set.seed((s$seed + i) %% .Machine$integer.max)
    chosen <- sample(ids, length(ids), replace = TRUE)
    blocks <- lapply(seq_along(chosen), function(j) {
      m <- data.table::copy(metadata[cluster_id == chosen[j]])
      old <- m$sample_id; new <- paste0("draw", j, "_", seq_along(old))
      zz <- data.table::copy(z[sample_id %in% old])
      zz[, sample_id := new[match(sample_id, old)]]
      # Authoritative metadata are rejoined by ex_long.
      zz[, (intersect(c("cluster_id", "sex", "ancestry", "cross"), names(zz))) := NULL]
      m[, `:=`(sample_id = new, cluster_id = paste0("draw_family_", j))]
      list(z = zz, m = m)
    })
    ans <- tryCatch(ex_mediation(data.table::rbindlist(lapply(blocks, `[[`, "z")),
      data.table::rbindlist(lapply(blocks, `[[`, "m")), cfg, s), error = function(e) ex_failure("bootstrap_error"))
    if (ans$row$status == "eligible") cbind(data.table::data.table(draw = i, status = "eligible"),
      ans$row[, ..metrics]) else data.table::data.table(draw = i, status = ans$row$status)
  }), fill = TRUE)
  nvalid <- sum(draws$status == "eligible")
  summary <- data.table::data.table(B = B, n_valid = nvalid, valid_fraction = nvalid / B,
    bootstrap_status = if (nvalid / B >= s$bootstrap_min_valid) "available" else "insufficient_valid_refits")
  for (nm in metrics) {
    v <- if (nm %in% names(draws)) draws[[nm]][draws$status == "eligible"] else numeric()
    ci <- if (nvalid / B >= s$bootstrap_min_valid && length(v) && all(is.finite(v)))
      as.numeric(stats::quantile(v, c(.025, .975))) else c(NA_real_, NA_real_)
    summary[, (paste0(nm, c("_boot_low", "_boot_high"))) := as.list(ci)]
  }
  list(summary = summary, draws = draws)
}
