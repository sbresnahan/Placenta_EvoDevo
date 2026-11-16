# Summary-table adapter for the ancestry-first PANTRY GxE aggregate contract.
gxe_specs <- function(cfg, z) {
  specs <- reg_context_manifest(cfg, z, "gxe_tables")
  if (is.null(specs)) return(NULL)
  require_columns(specs, c("tier", "source_format", "coordinate_system"))
  expected <- c(tier1 = "gxe_tier1_results.tsv.gz", tier1_meta = "gxe_tier1_meta.tsv.gz",
    tier2 = "gxe_tier2_results.tsv.gz")
  if (any(!specs$tier %in% names(expected)) || anyDuplicated(specs$tier) ||
      any(basename(specs$path) != expected[specs$tier]) ||
      any(specs$source_format != "CHROM:POS:REF:ALT") ||
      any(specs$coordinate_system != "one_based"))
    stop("GxE input requires unique authoritative tier aggregates and explicit one-based variant format")
  if (is.null(z$human_context$exposure_manifest)) stop("Supply ppBMI exposure_manifest")
  e <- reg_read(reg_path(z$human_context$exposure_manifest, cfg))
  require_columns(e, c("exposure", "exposure_type", "transformation", "scale_definition",
    "range_definition", "producing_config", "producing_provenance", "enabled"))
  e <- e[exposure == "ppBMI"]
  if (nrow(e) != 1L || !isTRUE(e$enabled) ||
      !e$exposure_type %in% c("continuous", "binary"))
    stop("Exposure manifest requires one enabled ppBMI row and explicit continuous/binary definition")
  for (k in c("transformation", "scale_definition", "range_definition"))
    reg_text(e[[k]], paste("ppBMI", k))
  for (k in c("producing_config", "producing_provenance")) {
    e[, (k) := reg_path(e[[k]], cfg)]; reg_require_file(e[[k]])
  }
  list(tables = specs, exposure = e)
}
gxe_read <- function(spec, exposure, target, aliases = NULL) {
  raw <- reg_read(spec$path)
  require_columns(raw, c("modality", "phenotype_id", "exposure"))
  raw <- raw[exposure == "ppBMI"]
  if (!nrow(raw)) stop("Configured GxE aggregate has no ppBMI rows: ", spec$input_id)
  tier <- spec$tier
  if (tier == "tier1") {
    require_columns(raw, c("ancestry", "variant_id", "b_gi", "b_gi_se", "b_g", "b_e",
      "pval_nominal", "pval_beta", "qval"))
    ancestry <- raw$ancestry; variant <- raw$variant_id
    beta <- raw$b_gi; se <- raw$b_gi_se
    support <- ifelse(is.na(raw$qval), NA, raw$qval <= .05)
    tested <- is.finite(raw$pval_beta)
  } else {
    require_columns(raw, c("beta_meta", "se_meta"))
    ancestry <- rep("cross_ancestry", nrow(raw))
    if (tier == "tier1_meta") {
      require_columns(raw, c("meta_variant_id", "scope", "lead_concordant",
        "p_feature_acat", "q_feature_acat", "n_ancestries"))
      if (any(!raw$scope %in% c("ancestry_only", "shared_concordant_lead", "shared_discordant_lead")))
        stop("Unknown tier1 synthesis scope")
      discordant <- !raw$lead_concordant %in% TRUE
      if (any(discordant & (!is.na(raw$meta_variant_id) | !is.na(raw$beta_meta) | !is.na(raw$se_meta))))
        stop("Discordant or ancestry-only synthesis cannot declare a common variant effect")
      variant <- raw$meta_variant_id
      tested <- is.finite(raw$p_feature_acat)
    } else {
      require_columns(raw, c("variant_id", "p_meta", "q_meta"))
      variant <- raw$variant_id; tested <- is.finite(raw$p_meta)
    }
    beta <- raw$beta_meta; se <- raw$se_meta; support <- rep(NA, nrow(raw))
  }
  for (k in intersect(c("pval_nominal", "pval_beta", "qval", "p_feature_acat",
      "q_feature_acat", "p_meta", "q_meta"), names(raw))) {
    if ((!is.numeric(raw[[k]]) && !all(is.na(raw[[k]]))) || any(!is.na(raw[[k]]) &
        (!is.finite(raw[[k]]) | raw[[k]] < 0 | raw[[k]] > 1)))
      stop("Invalid GxE probability field: ", k)
  }
  if (any(!is.na(se) & (!is.finite(se) | se <= 0))) stop("Invalid interaction SE")
  valid <- !is.na(variant) & nzchar(variant)
  if (any(!grepl("^[^:]+:[1-9][0-9]*:[ACGT]+:[ACGT]+$", variant[valid])))
    stop("Variant does not match declared CHROM:POS:REF:ALT format")
  parts <- tstrsplit(variant, ":", fixed = TRUE)
  getpart <- function(i) if (length(parts) >= i) parts[[i]] else rep(NA_character_, nrow(raw))
  chrom <- reg_alias(getpart(1), target, aliases)
  pos <- as.numeric(getpart(2)); ref <- getpart(3); alt <- getpart(4)
  if (any(valid & is.na(chrom))) stop("GxE contigs disagree with native reference")
  for (k in intersect(c("chrom", "pos", "ref", "alt"), names(raw))) {
    expected <- switch(k, chrom = chrom, pos = pos, ref = ref, alt = alt)
    supplied <- if (k == "chrom") reg_alias(as.character(raw[[k]]), target, aliases) else raw[[k]]
    if (any(valid & (is.na(supplied) | supplied != expected))) stop("GxE supplied coordinate/allele disagreement")
  }
  x <- data.table(input_id = spec$input_id, source_row = seq_len(nrow(raw)),
    context_type = "gxe_ppBMI", tier = tier, ancestry = ancestry, modality = raw$modality,
    phenotype_id = raw$phenotype_id, exposure = "ppBMI", kind = "lead",
    variant_id = variant, chrom = chrom, pos = pos, ref = ref, alt = alt,
    interaction_beta = beta, interaction_se = se, primary_source_support = support,
    source_tested = tested, synthesis_support = if (tier == "tier1_meta")
      ifelse(is.na(raw$q_feature_acat), NA, raw$q_feature_acat <= .05) else NA,
    exploratory_source_flag = if (tier == "tier2")
      ifelse(is.na(raw$q_meta), NA, raw$q_meta <= .05) else NA,
    statistical_scope = switch(tier, tier1 = "feature_cis_search_calibrated_ancestry_modality_exposure",
      tier1_meta = "feature_ACAT_cross_ancestry_synthesis",
      tier2 = "prioritized_min_cis_meta_p_then_BH_within_cis_search_uncalibrated"))
  if (anyNA(x[, .(modality, phenotype_id, ancestry)])) stop("Missing GxE feature keys")
  assert_unique(x, c("modality", "phenotype_id", "ancestry"), "GxE feature rows")
  original <- copy(raw); setnames(original, paste0("source_", names(original)))
  cbind(x, original)
}
gxe_context <- function(cfg, z, tx, map, repeats, target, aliases, s) {
  spec <- gxe_specs(cfg, z)
  if (is.null(spec)) {
    p <- nb_path(s, "context", "gxe_resource_coverage.tsv")
    nb_write(data.table(context_type = "gxe_ppBMI", status = "not_configured"), p)
    return(p)
  }
  q <- rbindlist(lapply(seq_len(nrow(spec$tables)), function(i)
    gxe_read(spec$tables[i], spec$exposure, target, aliases)), fill = TRUE)
  if (is.null(z$human_context$phenotype_map)) stop("ppBMI context requires verified phenotype_map")
  mapped <- reg_context_map(q, reg_read(reg_path(z$human_context$phenotype_map, cfg)),
    unique(map[, .(event_id, gene_id)]))
  links <- merge(mapped$links, map[, .(event_id, transcript_id)], by = "event_id", allow.cartesian = TRUE)
  links[, exact_transcript_phenotype := nb_pantry_transcript(modality, phenotype_id) ==
    sub("\\.[0-9]+$", "", transcript_id)]
  links[is.na(exact_transcript_phenotype), exact_transcript_phenotype := FALSE]
  if (!is.null(repeats) && nrow(links)) {
    setnames(links, "event_id", "context_event_id"); links[, event_id := transcript_id]
    unit <- copy(tx); unit[, event_id := transcript_id]
    links <- reg_context_variant_links(links, unit, repeats, target, aliases)
    links[, event_id := context_event_id][, context_event_id := NULL]
  }
  for (k in c("variant_overlaps_TSS_TE", "variant_overlaps_promoter_TE",
      "variant_overlaps_proximal_only_TE"))
    if (!k %in% names(links)) links[, (k) := NA]
  strata <- unique(q[, .(input_id, context_type, tier, ancestry, modality, exposure, statistical_scope)])
  summary <- cbind(tx[rep(seq_len(nrow(tx)), each = nrow(strata)),
    .(transcript_id, gene_id, promoter_te_label, te_tss_window_bp, tss_te_label)],
    strata[rep(seq_len(nrow(strata)), times = nrow(tx))])
  summary[, `:=`(assessment_status = "untested_or_unmapped", n_tested_features = 0L,
    exact_QTL = NA, event_QTL = NA, gene_QTL = NA, same_TE_lead = NA,
    proximal_same_TE_lead = NA, cCRE = NA)]
  keys <- c("transcript_id", "input_id", "tier", "ancestry", "modality", "exposure")
  if (nrow(links)) {
    # Only tier1 support populates the primary follow-up supports.
    a <- links[, {
      tested <- source_tested %in% TRUE
      ok <- primary_source_support %in% TRUE
      assessed <- any(tested) && tier[1] == "tier1" && any(!is.na(primary_source_support[tested]))
      flag <- function(v) if (assessed) any(ok & v) else NA
      list(assessment_status = if (any(tested)) "assessed_mapped" else "missing_statistic",
        n_tested_features = uniqueN(phenotype_id[tested]),
        exact_QTL = if (any(tested & exact_transcript_phenotype)) flag(exact_transcript_phenotype) else NA,
        event_QTL = if (any(tested & link_level == "verified_phenotype_event"))
          flag(link_level == "verified_phenotype_event") else NA,
        gene_QTL = if (any(tested & link_level != "verified_phenotype_event"))
          flag(link_level != "verified_phenotype_event") else NA,
        same_TE_lead = if (is.null(repeats)) NA else flag(variant_overlaps_TSS_TE %in% TRUE),
        proximal_same_TE_lead = if (is.null(repeats)) NA else flag(variant_overlaps_proximal_only_TE %in% TRUE))
    }, by = keys]
    summary[a, on = keys, `:=`(assessment_status = i.assessment_status,
      n_tested_features = i.n_tested_features, exact_QTL = i.exact_QTL,
      event_QTL = i.event_QTL, gene_QTL = i.gene_QTL, same_TE_lead = i.same_TE_lead,
      proximal_same_TE_lead = i.proximal_same_TE_lead)]
  }
  summary[!promoter_te_label %in% c("TE", "nonTE"), proximal_same_TE_lead := NA]
  summary[!tss_te_label %in% c("TE", "nonTE"), same_TE_lead := NA]
  resource <- q[, .(n_source_features = .N, n_tested = sum(source_tested),
    n_missing_statistics = sum(!source_tested)), by = .(input_id, tier, ancestry, modality, exposure)]
  files <- nb_path(s, "context", c("gxe_ppBMI_context.tsv.gz", "gxe_ppBMI_links.tsv.gz",
    "gxe_mapping_coverage.tsv", "gxe_resource_coverage.tsv", "gxe_exposure_provenance.tsv"))
  for (i in seq_along(files)) nb_write(list(summary, links, mapped$coverage, resource,
    spec$exposure)[[i]], files[i])
  files
}
