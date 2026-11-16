ex_finish <- function(d) {
  if (!nrow(d)) return(d)
  d <- data.table::copy(d)
  ex_need(d, c("module", "cohort_id", "level", "predictor", "exposure"))
  families <- c("module", "cohort_id", "level", "predictor")
  for (pair in list(c("slope_p", "slope_q"), c("indirect_path_p", "indirect_q"),
                    c("c_prime_p", "direct_q"), c("modification_p", "modification_q"))) {
    if (!pair[1] %in% names(d)) next
    p <- pair[1]; q <- pair[2]
    d[, (q) := stats::p.adjust(get(p), method = "BH", n = .N), by = families]
    d[, (paste0(q, "_within_exposure")) := stats::p.adjust(get(p), method = "BH", n = .N),
      by = c(families, "exposure")]
    if (p == "slope_p") d[module == "h6", slope_q_pooled_cohorts :=
      if (data.table::uniqueN(cohort_id) > 1L) stats::p.adjust(slope_p, "BH", n = .N) else NA_real_,
      by = .(level, predictor)]
  }
  d[, primary_p := if ("indirect_path_p" %in% names(d))
    data.table::fifelse(module == "h7_mediation", indirect_path_p,
      if ("modification_p" %in% names(d)) data.table::fifelse(module == "h7_modification", modification_p, slope_p)
      else slope_p) else slope_p]
  d[, primary_q := NA_real_]
  if ("slope_q" %in% names(d)) d[module %in% c("h6", "h7_total"), primary_q := slope_q]
  if ("indirect_q" %in% names(d)) d[module == "h7_mediation", primary_q := indirect_q]
  if ("modification_q" %in% names(d)) d[module == "h7_modification", primary_q := modification_q]
  d[, discovery := status == "eligible" & is.finite(primary_q) & primary_q < .05]
  d
}
ex_gene_scores <- function(d, spans) {
  groups <- intersect(c("module", "cohort_id", "level", "predictor", "exposure", "gene_id"), names(d))
  d[, gene_length := spans$end[match(gene_id, spans$gene_id)] - spans$start[match(gene_id, spans$gene_id)] + 1]
  d[, {
    p <- primary_p[is.finite(primary_p) & status == "eligible"]
    list(gene_p = if (length(p)) min(1, .N * min(p)) else NA_real_,
      n_tested = .N, n_estimable = length(p),
      total_counts = if (all(is.na(total_counts))) NA_real_ else stats::median(total_counts, na.rm = TRUE),
      n_biological_units = if (all(is.na(n_biological_units))) NA_real_ else stats::median(n_biological_units, na.rm = TRUE),
      n_snps = if (all(is.na(n_snps))) NA_real_ else stats::median(n_snps, na.rm = TRUE),
      gene_length = gene_length[1])
  }, by = groups]
}
ex_annotation_check <- function(manifest, memberships, cfg) {
  ex_need(manifest, c("set_id", "category", "name", "species", "source", "source_version", "retrieved_date"))
  ex_unique(manifest, c("set_id", "species"), "annotation manifest")
  for (nm in c("set_id", "category", "name", "species", "source", "source_version", "retrieved_date"))
    if (any(!ex_text(manifest[[nm]]))) stop("Undeclared annotation metadata: ", nm)
  ex_need(memberships, c("set_id", "gene_id", "species"))
  ex_unique(memberships, c("set_id", "gene_id", "species"), "gene-set membership")
  if (any(!paste(memberships$species, memberships$set_id) %in%
      paste(manifest$species, manifest$set_id))) stop("Unregistered gene set/species")
  if (any(memberships$species != cfg$species)) stop("Gene sets must use species-native IDs or a documented orthology map")
  if (cfg$species != "human" && "orthology_source" %in% names(manifest) &&
      any(manifest$species == cfg$species & grepl("human", manifest$source) & !ex_text(manifest$orthology_source)))
    stop("Document cross-species orthology source")
  invisible(TRUE)
}
ex_match_strata <- function(d) {
  vars <- c("total_counts", "n_biological_units", "n_snps", "gene_length", "n_tested")
  ex_need(d, vars)
  if (anyNA(d[, ..vars]) || any(!is.finite(as.matrix(d[, ..vars]))))
    stop("Missing detectability covariates; do not impute annotation-matching inputs")
  bins <- lapply(d[, ..vars], function(v) {
    if (data.table::uniqueN(v) < 2L) rep(1L, length(v)) else
      1L + as.integer(rank(v, ties.method = "average") > length(v) / 2)
  })
  # Coarsen sparse strata deterministically, without using gene-set membership.
  k <- length(bins)
  repeat {
    group <- if (k) do.call(paste, c(bins[seq_len(k)], sep = "_")) else rep("all", nrow(d))
    if (min(table(group)) >= 5L || k == 0L) break
    k <- k - 1L
  }
  list(stratum = group, matching_variables_retained = paste(vars[seq_len(k)], collapse = ","))
}
ex_enrich_group <- function(d, memberships, manifest, s) {
  d <- d[is.finite(gene_p)]
  if (nrow(d) < 10L) return(data.table::data.table(status = "insufficient_tested_genes"))
  ex_unique(d, "gene_id", "enrichment gene universe")
  strat <- ex_match_strata(d); d[, stratum := strat$stratum]
  # Continuous, unsigned gene evidence; one score per gene.
  score <- rank(-log10(pmax(d$gene_p, 1e-300)), ties.method = "average") / nrow(d)
  sets <- intersect(manifest$set_id, memberships$set_id)
  data.table::rbindlist(lapply(sets, function(set) {
    in_set <- d$gene_id %in% memberships[set_id == set]$gene_id
    n <- sum(in_set); n_bg <- nrow(d) - n
    if (n < 5L || n_bg < 5L || n > 2000L) return(data.table::data.table(set_id = set,
      n_set = n, n_background = n_bg, status = "insufficient_set_or_background"))
    obs <- mean(score[in_set]) - mean(score[!in_set])
    strata <- split(seq_len(nrow(d)), d$stratum)
    set.seed(s$seed)
    perm <- replicate(s$permutations, {
      ii <- unlist(lapply(strata, function(idx) {
        k <- sum(in_set[idx]); if (k) sample(idx, k, replace = FALSE) else integer()
      }), use.names = FALSE)
      mean(score[ii]) - mean(score[-ii])
    })
    data.table::data.table(set_id = set, n_set = n, n_background = n_bg,
      evidence_rank_difference = obs, p = (1 + sum(perm >= obs)) / (1 + length(perm)),
      permutations = s$permutations, matching_variables = strat$matching_variables_retained,
      status = "available")
  }), fill = TRUE)
}
ex_enrich <- function(scores, memberships, manifest, cfg, s) {
  ex_annotation_check(manifest, memberships, cfg)
  memberships <- memberships[species == cfg$species]; manifest <- manifest[species == cfg$species]
  groups <- c("module", "cohort_id", "level", "predictor", "exposure")
  out <- scores[, ex_enrich_group(.SD, memberships, manifest, s), by = groups]
  if ("p" %in% names(out)) out[, q := stats::p.adjust(p, "BH"), by = .(module, cohort_id, level, predictor)]
  out
}
ex_h5_adapter <- function(context) {
  ex_need(context, c("transcript_id", "gene_id", "n_exact_isoform_qtl_links", "n_event_qtl_links",
    "n_qtl_variants_in_TSS_TE", "n_peak_overlaps", "qtl_assessed", "chromatin_assessed",
    "TE_insertion_assessed"))
  ex_unique(context, "transcript_id", "H5 context")
  out <- context[, .(feature_id = transcript_id, gene_id,
    exact_QTL = data.table::fifelse(qtl_assessed, n_exact_isoform_qtl_links > 0, NA),
    event_QTL = data.table::fifelse(qtl_assessed, n_event_qtl_links > 0, NA),
    same_TE_lead = data.table::fifelse(qtl_assessed & TE_insertion_assessed, n_qtl_variants_in_TSS_TE > 0, NA),
    cCRE = data.table::fifelse(chromatin_assessed, n_peak_overlaps > 0, NA))]
  if ("n_qtl_variants_in_proximal_TE" %in% names(context)) out[, proximal_same_TE_lead :=
    data.table::fifelse(context$qtl_assessed & context$TE_insertion_assessed,
      context$n_qtl_variants_in_proximal_TE > 0, NA)]
  out
}
ex_h5_followup <- function(results, annotations, cfg, supports = intersect(
    c("exact_QTL", "event_QTL", "same_TE_lead", "proximal_same_TE_lead", "cCRE"), names(annotations))) {
  ex_require_human(cfg)
  ex_unique(annotations, "feature_id", "H5 annotations")
  groups <- c("module", "cohort_id", "level", "predictor", "exposure")
  d <- merge(results[level == "transcript" & status == "eligible"], annotations, by = c("feature_id", "gene_id"))
  if (!nrow(d)) return(data.table::data.table(status = "no_assessable_transcripts"))
  out <- data.table::rbindlist(lapply(supports, function(nm) {
    d[!is.na(get(nm)), {
      z <- data.table::copy(.SD); z[, support := as.integer(get(nm))]
      # Multiple family tests remain clustered within their host gene.
      base <- list(support = nm, n_features = data.table::uniqueN(z$feature_id),
        n_genes = data.table::uniqueN(z$gene_id), n_supported = sum(z$support),
        n_discoveries = sum(z$discovery),
        mean_evidence_supported = if (any(z$support == 1)) mean(-log10(pmax(z$primary_p[z$support == 1], 1e-300))) else NA_real_,
        mean_evidence_unsupported = if (any(z$support == 0)) mean(-log10(pmax(z$primary_p[z$support == 0], 1e-300))) else NA_real_)
      if (length(unique(z$support)) < 2L || sum(z$discovery) < 5L ||
          sum(!z$discovery) < 5L || data.table::uniqueN(z$gene_id) < 10L) {
        as.list(c(base, list(status = "sparse_or_constant_annotation")))
      } else {
        warning_text <- character()
        fit <- tryCatch(withCallingHandlers(stats::glm(discovery ~ support + log1p(total_counts) +
          log1p(n_snps) + n_biological_units, data = z, family = stats::binomial()),
          warning = function(w) {warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning")}),
          error = function(e) e)
        if (inherits(fit, "error") || length(warning_text) || !isTRUE(fit$converged) ||
            any(abs(stats::coef(fit)) > 20) || any(!is.finite(stats::coef(fit)))) {
          as.list(c(base, list(status = "separated_or_failed_annotation_model")))
        } else {
          V <- tryCatch(sandwich::vcovCL(fit, cluster = z$gene_id, type = "HC1"), error = function(e) NULL)
          if (is.null(V) || !is.finite(V["support", "support"]) || V["support", "support"] <= 0)
            as.list(c(base, list(status = "unavailable_cluster_uncertainty"))) else {
              b <- stats::coef(fit)[["support"]]; se <- sqrt(V["support", "support"])
              as.list(c(base, list(status = "available", odds_ratio = exp(b),
                low = exp(b - 1.96 * se), high = exp(b + 1.96 * se), p = 2 * stats::pnorm(-abs(b / se)))))
            }
        }
      }
    }, by = groups]
  }), fill = TRUE)
  if ("p" %in% names(out)) out[, q := stats::p.adjust(p, "BH"), by = .(module, cohort_id, level, predictor)]
  out
}
ex_h5_gxe_followup <- function(results, context, cfg) {
  ex_need(context, c("transcript_id", "gene_id", "input_id", "tier", "ancestry",
    "modality", "exposure", "context_type", "exact_QTL", "event_QTL", "gene_QTL",
    "same_TE_lead", "proximal_same_TE_lead"))
  context <- data.table::copy(context)
  data.table::setnames(context, "exposure", "context_exposure")
  groups <- c("input_id", "tier", "ancestry", "modality", "context_exposure", "context_type")
  # Each ancestry/modality/tier has its own assessment and multiplicity family.
  context[, {
    ann <- data.table::copy(.SD)
    data.table::setnames(ann, "transcript_id", "feature_id")
    ex_h5_followup(results, ann, cfg, supports = c("exact_QTL", "event_QTL", "gene_QTL",
      "same_TE_lead", "proximal_same_TE_lead"))
  }, by = groups]
}
