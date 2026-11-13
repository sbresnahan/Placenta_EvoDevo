choose_sentinels <- function(counts, pairs, min_depth = 20) {
  require_columns(counts, c("sample_id","family_id","event_id","locus_id",
                            "maternal_count","paternal_count","ref_maternal"))
  c <- data.table::copy(counts)
  assert_unique(c, c("sample_id","event_id","locus_id"), "Event locus counts")
  c[, depth := maternal_count + paternal_count]
  c <- c[!is.na(ref_maternal) & is.finite(depth) & depth >= min_depth]
  # Choose one offspring per family globally by callability, before imbalance.
  max_event <- c[,.(depth=max(depth)),by=.(sample_id,family_id,event_id)]
  samples <- max_event[, .(callable_events = uniqueN(event_id), total_depth = sum(depth)),
               by = .(sample_id, family_id)]
  data.table::setorder(samples, family_id, -callable_events, -total_depth, sample_id)
  keep <- samples[, .SD[1], by = family_id]$sample_id
  c <- c[sample_id %in% keep]
  results <- list()
  for (i in seq_len(nrow(pairs))) {
    p <- pairs[i]
    aa <- c[event_id == p$event_a]
    bb <- c[event_id == p$event_b]
    if (!nrow(aa) || !nrow(bb)) next
    ab <- merge(aa, bb, by = c("sample_id","family_id"), suffixes = c("_a","_b"),
                allow.cartesian = TRUE)
    if (!"distinguishing_coordinate_a" %in% names(ab))
      ab[, `:=`(distinguishing_coordinate_a=FALSE,distinguishing_coordinate_b=FALSE)]
    score <- ab[, .(direct_coordinates = as.integer(all(distinguishing_coordinate_a) &&
                       all(distinguishing_coordinate_b)), n_donors = uniqueN(sample_id),
                     min_depth = min(pmin(depth_a, depth_b))),
                 by = .(locus_id_a, locus_id_b)]
    if (!nrow(score)) next
    # Do not prefer a direct-coordinate pair below the primary donor floor.
    # Sample size gating is still applied in the model, not guessed here.
    data.table::setorder(score, -n_donors, -direct_coordinates, -min_depth, locus_id_a, locus_id_b)
    sel <- score[1]
    d <- ab[locus_id_a == sel$locus_id_a & locus_id_b == sel$locus_id_b]
    d[, `:=`(gene_id = p$gene_id, gene_name = p$gene_name,
              pair_type = p$pair_type, event_a = p$event_a, event_b = p$event_b)]
    results[[length(results) + 1L]] <- d
  }
  ans<-data.table::rbindlist(results, fill = TRUE)
  if(!ncol(ans)) ans<-data.table::data.table(gene_id=character(),sample_id=character(),
    family_id=character(),event_a=character(),event_b=character())
  ans
}
hc3_intercept <- function(y, X) {
  if(any(!is.finite(y))||any(!is.finite(X)))return(NULL)
  if (nrow(X) - ncol(X) < 5L || qr(X)$rank != ncol(X)) return(NULL)
  fit <- lm.fit(X, y)
  b <- solve(crossprod(X))
  h <- rowSums((X %*% b) * X)
  if (any(h >= 1 - 1e-8)) return(NULL)
  v <- b %*% crossprod(X, X * (fit$residuals / (1 - h))^2) %*% b
  se <- sqrt(v[1,1])
  if (!is.finite(se) || se <= 0) return(NULL)
  df <- nrow(X) - ncol(X)
  list(estimate = unname(fit$coefficients[1]), se = se, df = df,
       p = 2 * pt(-abs(fit$coefficients[1] / se), df))
}
fit_promoter_pair <- function(d, covariates = character(), n_boot = 2000L,
                              min_donors = 10L, min_orientation = 3L, seed = 1147L) {
  assert_unique(d,"family_id","One offspring per family at model entry")
  assert_unique(d,"sample_id","One row per offspring at model entry")
  require_columns(d,c("maternal_count_a","paternal_count_a","maternal_count_b","paternal_count_b"))
  values<-unlist(d[,.(maternal_count_a,paternal_count_a,maternal_count_b,paternal_count_b)])
  if(any(!is.finite(values))||any(values<0))stop("Invalid paired counts")
  base <- data.table::data.table(gene_id = d$gene_id[1], gene_name = d$gene_name[1],
    event_a = d$event_a[1], event_b = d$event_b[1], pair_type = d$pair_type[1],
    n_families = data.table::uniqueN(d$family_id), status = "insufficient_donors")
  if (nrow(d) < min_donors) return(base)
  for (nm in c("ref_maternal_a","ref_maternal_b")) {
    if (anyNA(d[[nm]]) || min(table(factor(d[[nm]], levels = c(FALSE,TRUE)))) < min_orientation) {
      base$status <- "insufficient_cis_orientation"; return(base)
    }
  }
  require_columns(d, covariates, "Model covariates")
  z <- if (length(covariates)) as.data.frame(d[, ..covariates]) else
    data.frame(row.names = seq_len(nrow(d)))
  if (length(covariates) && anyNA(z)) {
    base$status <- "missing_covariates"; return(base)
  }
  d <- data.table::copy(d)
  d[, contrast := log((maternal_count_a + 0.5) / (paternal_count_a + 0.5)) -
                    log((maternal_count_b + 0.5) / (paternal_count_b + 0.5))]
  z$cis_a <- ifelse(d$ref_maternal_a, 1, -1)
  z$cis_b <- ifelse(d$ref_maternal_b, 1, -1)
  # Perfect LD gives redundant cis terms; use the same cis span once, explicitly.
  if (identical(z$cis_a, z$cis_b) || identical(z$cis_a, -z$cis_b)) {
    z$cis_b <- NULL
    cis_basis <- "one_redundant_cis_term_removed"
  } else cis_basis <- "two_cis_terms"
  for (nm in covariates) if (is.numeric(z[[nm]])) z[[nm]] <- z[[nm]] - mean(z[[nm]])
  X <- tryCatch(model.matrix(~ ., z),error=function(e)NULL)
  if(is.null(X)) {base$status<-"unidentifiable_covariates";return(base)}
  # Categorical covariates are also centered at the analyzed cohort distribution.
  cn <- setdiff(colnames(X), c("(Intercept)", "cis_a", "cis_b"))
  if (length(cn)) X[,cn] <- scale(X[,cn,drop=FALSE], center = TRUE, scale = FALSE)
  m <- hc3_intercept(d$contrast, X)
  if (is.null(m)) {base$status <- "unidentifiable_or_leverage"; return(base)}
  set.seed(seed)
  boot <- replicate(n_boot, {
    i <- sample.int(nrow(d), nrow(d), replace = TRUE)
    if (qr(X[i,,drop=FALSE])$rank != ncol(X)) NA_real_ else
      unname(lm.fit(X[i,,drop=FALSE], d$contrast[i])$coefficients[1])
  })
  frac <- mean(is.finite(boot))
  ci <- if (frac >= 0.9) quantile(boot, c(.025,.975), na.rm = TRUE) else c(NA,NA)
  base[, `:=`(estimate = m$estimate, se = m$se, p_value = if (frac >= .9) m$p else NA_real_,
    ci_low = unname(ci[1]), ci_high = unname(ci[2]), residual_df = m$df,
    bootstrap_valid_fraction = frac, cis_basis = cis_basis,
    status = ifelse(frac >= .9, "eligible_conditional_inference", "unstable_bootstrap"),
    uncertainty = "between_families_conditional_on_Ornaments")]
  base
}
fit_all_pairs <- function(paired, covariates = character(), n_boot = 2000L,
                          min_donors = 10L, min_orientation = 3L, seed = 1147L) {
  if (!nrow(paired)) return(data.table::data.table(
    gene_id = character(), status = character(), p_value = numeric(),
    q_value = numeric(), estimate = numeric(), discordant = logical()))
  out <- data.table::rbindlist(lapply(split(paired, paired$gene_id), function(d)
    fit_promoter_pair(d, covariates, n_boot, min_donors, min_orientation, seed)), fill = TRUE)
  if (!"p_value" %in% names(out)) out[, p_value := NA_real_]
  out[, q_value := p.adjust(p_value, "BH")]
  if (!"estimate" %in% names(out)) out[, estimate := NA_real_]
  out[, discordant := !is.na(q_value) & q_value <= .05 & abs(estimate) >= log(2)]
  out
}
equivalence_call <- function(estimate, se, df, odds_margin = 1.25) {
  lo <- estimate - qt(.95, df) * se; hi <- estimate + qt(.95, df) * se
  data.table::data.table(ci90_low = lo, ci90_high = hi,
    equivalent = lo > -log(odds_margin) & hi < log(odds_margin))
}
matched_te_model <- function(x, match_variables, min_genes = 30L, max_distance = Inf) {
  require_columns(x, c("gene_id","pair_type","discordant","estimate","status", match_variables))
  assert_unique(x, "gene_id", "TE model genes")
  x <- x[status == "eligible_conditional_inference" &
    complete.cases(x[, c("estimate", match_variables), with = FALSE])]
  if (!length(match_variables)) stop("Explicit matching variables required")
  treated <- which(x$pair_type == "TE_nonTE")
  controls <- which(x$pair_type == "nonTE_nonTE")
  if (min(length(treated), length(controls)) < min_genes)
    return(list(summary=data.table::data.table(status = "below_gene_feasibility_floor"),
                matches=data.table::data.table()))
  Z <- scale(as.matrix(x[, ..match_variables]))
  if (any(!is.finite(Z))) stop("Constant or nonnumeric matching variable")
  used <- integer(); sets <- list()
  for (j in treated[order(x$gene_id[treated])]) {
    available <- setdiff(controls, used)
    if (!length(available)) break
    dist <- rowSums((Z[available,,drop=FALSE] -
                      matrix(Z[j,], nrow = length(available), ncol = ncol(Z), byrow = TRUE))^2)
    k <- available[order(dist, x$gene_id[available])[1]]
    if (sqrt(min(dist)) > max_distance) next
    used <- c(used, k)
    tmp <- data.table::copy(x[c(j,k)])
    tmp[, `:=`(matched_set = length(sets) + 1L, exposed = c(1L,0L))]
    sets[[length(sets) + 1L]] <- tmp
  }
  matched <- data.table::rbindlist(sets)
  if (length(sets) < min_genes)
    return(list(summary = data.table::data.table(status = "below_matched_floor"),
                matches = matched))
  delta <- matched[, abs(estimate[exposed == 1L]) - abs(estimate[exposed == 0L]),
                   by = matched_set]$V1
  continuous <- if (sd(delta) > 0) t.test(delta) else NULL
  # Binding strata locally avoids depending on the user's attached packages.
  strata <- survival::strata
  fit <- tryCatch(survival::clogit(discordant ~ exposed + strata(matched_set),
    data = matched, control = survival::coxph.control(iter.max = 50)),
    error = function(e) NULL, warning = function(w) NULL)
  valid <- !is.null(fit) && all(is.finite(coef(fit))) &&
    all(is.finite(vcov(fit))) && all(diag(vcov(fit)) > 0)
  ci <- if (valid) suppressMessages(confint(fit)) else c(NA_real_,NA_real_)
  list(summary = data.table::data.table(
    status = if (valid) "estimated_review_balance" else "clogit_unestimable",
    n_sets = length(sets),
    log_odds_ratio = if (valid) unname(coef(fit)[1]) else NA_real_,
    ci_low = ci[1], ci_high = ci[2],
    p_value = if (valid) summary(fit)$coefficients[1,"Pr(>|z|)"] else NA_real_,
    continuous_abs_effect_difference = mean(delta),
    continuous_ci_low = if (is.null(continuous)) NA_real_ else continuous$conf.int[1],
    continuous_ci_high = if (is.null(continuous)) NA_real_ else continuous$conf.int[2]),
    matches = matched, model = fit)
}
