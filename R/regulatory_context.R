reg_context_manifest <- function(cfg, z, key) {
  value <- z$human_context[[key]]
  if (is.null(value)) return(NULL)
  d <- reg_read(reg_path(value, cfg))
  require_columns(d, c("input_id", "path", "build"))
  assert_unique(d, "input_id", paste("human_context", key))
  d[, path := vapply(path, reg_path, "", cfg = cfg)]
  if (any(!d$build %in% c("hg38", "GRCh38"))) stop("Human context requires hg38 inputs")
  for (p in d$path) if (!file.exists(p) || file.info(p)$size <= 0)
    stop("Requested human context input unavailable: ", p)
  d
}
reg_context_qtl <- function(spec) {
  d <- reg_read(spec$path)
  col <- function(key, alternatives, required = FALSE) {
    override <- paste0(key, "_column")
    named <- if (override %in% names(spec) && !is.na(spec[[override]])) spec[[override]] else NULL
    use <- if (!is.null(named)) named else intersect(alternatives, names(d))[1]
    if (is.na(use) || is.null(use) || !use %in% names(d)) {
      if (required) stop(spec$input_id, ": missing ", key, " column; declare ", override)
      return(rep(NA_character_, nrow(d)))
    }
    d[[use]]
  }
  require_columns(spec, "kind", "QTL manifest")
  if (!spec$kind %in% c("lead", "pip", "credible_set", "coloc"))
    stop("QTL kind must be lead, pip, credible_set or coloc")
  x <- data.table(input_id = spec$input_id, source_row = seq_len(nrow(d)),
    kind = spec$kind, modality = as.character(col("modality", "modality", TRUE)),
    phenotype_id = as.character(col("phenotype", c("phenotype_id", "trait"), TRUE)),
    ancestry = as.character(col("ancestry", c("ancestry", "sig_in"))),
    variant_id = as.character(col("variant", c("variant_id", "snp", "lead_snp"))),
    chrom = as.character(col("chrom", c("chrom", "chr"))),
    pos = suppressWarnings(as.numeric(col("pos", c("pos", "lead_pos")))),
    pip_all = suppressWarnings(as.numeric(col("pip", "pip_all"))),
    pip_cs = suppressWarnings(as.numeric(col("pip_cs", "pip_cs"))),
    cs_index = as.character(col("cs_index", c("cs_index", "CSIndex"))),
    beta = suppressWarnings(as.numeric(col("beta", c("beta", "slope")))),
    se = suppressWarnings(as.numeric(col("se", c("se", "slope_se")))),
    ref = as.character(col("ref", "ref")),
    alt = as.character(col("alt", "alt")),
    allele0 = as.character(col("allele0", "a0")),
    allele1 = as.character(col("allele1", "a1")),
    cs_start = suppressWarnings(as.numeric(col("cs_start", "cs_start"))),
    cs_end = suppressWarnings(as.numeric(col("cs_end", "cs_end"))))
  if (anyNA(x[, .(modality, phenotype_id)])) stop("Missing QTL phenotype/modality keys")
  if ("ancestry" %in% names(spec) && !is.na(spec$ancestry))
    x[is.na(ancestry), ancestry := spec$ancestry]
  if (anyNA(x$ancestry))
    stop("Declare ancestry (e.g. joint or unspecified) in QTL manifest")
  if (spec$kind %in% c("lead", "pip")) {
    if (!"coordinate_system" %in% names(spec) ||
        !identical(spec$coordinate_system, "one_based"))
      stop("Variant coordinates require explicit coordinate_system=one_based")
    if (anyNA(x[, .(chrom, pos, variant_id)]))
      stop("QTL variant coordinates/ID missing; supply explicit column mappings, not ID guessing")
    if (any(!is.finite(x$pos) | x$pos < 1 | x$pos != floor(x$pos)))
      stop("Invalid QTL variant position")
    assert_unique(x, c("modality", "phenotype_id", "ancestry", "variant_id"),
      "QTL rows (declare ancestry for joint or unspecified inputs)")
  }
  for (v in c("pip_all", "pip_cs")) if (any(!is.na(x[[v]]) &
      (!is.finite(x[[v]]) | x[[v]] < 0 | x[[v]] > 1))) stop("Invalid PIP")
  # Keep every original field (including ancestry-specific effects and
  # colocalization probabilities) in the auditable long-form evidence.
  original <- copy(d); setnames(original, paste0("source_", names(original)))
  cbind(x, original)
}
reg_context_map <- function(q, map, events) {
  require_columns(map, c("modality", "phenotype_id", "gene_id", "event_id",
    "mapping_level", "mapping_status"), "Phenotype crosswalk")
  if (any(!map$mapping_level %in% c("event", "gene")) ||
      any(!map$mapping_status %in% c("verified", "ambiguous", "unmapped")))
    stop("Invalid crosswalk level/status")
  map <- unique(copy(map))
  map[, n_mapping := .N, by = .(modality, phenotype_id)]
  map[, mapping_usable := mapping_status == "verified" & n_mapping == 1L]
  valid <- map[mapping_usable == TRUE]
  if (anyNA(valid$gene_id)) stop("Verified crosswalk requires gene_id")
  ev <- valid[mapping_level == "event"]
  if (anyNA(ev$event_id) || any(!ev$event_id %in% events$event_id))
    stop("Verified event mapping refers to absent event")
  if (any(ev$gene_id != events$gene_id[match(ev$event_id, events$event_id)]))
    stop("Crosswalk event/gene disagreement")
  mapping_counts <- unique(map[, .(modality, phenotype_id, n_mapping)])
  coverage <- merge(unique(q[, .(input_id, modality, phenotype_id)]),
    mapping_counts, by = c("modality", "phenotype_id"), all.x = TRUE)
  flags <- map[, .(mapping_status = if (all(mapping_usable)) "verified" else
    if (.N > 1L || any(mapping_status == "ambiguous")) "ambiguous" else "unmapped"),
    by = .(modality, phenotype_id)]
  coverage <- merge(coverage, flags, by = c("modality", "phenotype_id"), all.x = TRUE)
  coverage[is.na(mapping_status), mapping_status := "unmapped"]
  a <- merge(q, valid, by = c("modality", "phenotype_id"), allow.cartesian = TRUE)
  direct <- a[mapping_level == "event"]
  direct[, link_level := "verified_phenotype_event"]
  gene <- a[mapping_level == "gene"]
  gene[, event_id := NULL]
  gene <- merge(gene, events[, .(gene_id, event_id)], by = "gene_id",
    allow.cartesian = TRUE)
  gene[, link_level := "gene_context_not_TSS_phenotype_validation"]
  list(links = rbindlist(list(direct, gene), fill = TRUE), coverage = coverage)
}
reg_context_overlap <- function(tx, specs, target, aliases = NULL) {
  empty <- data.table(event_id = character(), transcript_id = character(),
    input_id = character(), evidence_type = character(), peak_row = integer(),
    annotation_selection_used = logical())
  if (is.null(specs) || !nrow(specs)) return(list(hits = empty, coverage = data.table(
    input_id = character(), evidence_type = character(), n_peaks = integer(),
    n_tss_hits = integer())))
  require_columns(specs, c("input_id", "path", "evidence_type", "selection_used"))
  rows <- list(); coverage <- list()
  for (i in seq_len(nrow(specs))) {
    d <- reg_bed(specs$path[i])
    d[, chrom := reg_alias(chrom, target, aliases)]
    dropped <- sum(is.na(d$chrom)); d <- d[!is.na(chrom)]
    if (!nrow(d)) stop("No peak contigs match the human reference: ", specs$input_id[i])
    h <- reg_hits(tx[, .(chrom, start = tss, end = tss)], d)
    if (nrow(h)) {
      q <- cbind(tx[h$query, .(event_id, transcript_id)],
        d[h$subject, .(peak_chrom = chrom, peak_start = start, peak_end = end)])
      q[, `:=`(input_id = specs$input_id[i], evidence_type = specs$evidence_type[i],
        peak_row = h$subject, annotation_selection_used = specs$selection_used[i])]
      rows[[length(rows) + 1L]] <- q
    }
    coverage[[i]] <- data.table(input_id = specs$input_id[i],
      evidence_type = specs$evidence_type[i], n_peaks = nrow(d),
      unmapped_peak_contigs_records = dropped, n_tss_hits = nrow(h),
      annotation_selection_used = specs$selection_used[i])
  }
  list(hits = if (length(rows)) rbindlist(rows, fill = TRUE) else empty,
    coverage = rbindlist(coverage, fill = TRUE))
}
reg_context_variant_links <- function(links, tx, repeats, target, aliases = NULL) {
  if (!nrow(links)) return(links)
  links <- copy(links)
  links[, `:=`(variant_at_event_TSS = NA, variant_overlaps_TSS_TE = NA)]
  use <- which(links$kind %in% c("lead", "pip") & !is.na(links$pos))
  if (!length(use)) return(links)
  links[use, chrom := reg_alias(chrom, target, aliases)]
  if (anyNA(links$chrom[use])) stop("QTL variant contig does not match target reference")
  r <- copy(repeats); r[, chrom := reg_alias(chrom, target, aliases)]
  r <- r[!is.na(chrom) & is_te %in% TRUE]
  points <- tx[, .(chrom, start = tss, end = tss)]
  h <- reg_hits(points, r)
  er <- if (nrow(h)) unique(data.table(event_id = tx$event_id[h$query],
    repeat_row = h$subject)) else data.table(event_id = character(), repeat_row = integer())
  q <- links[use, .(chrom, start = pos, end = pos)]
  ht <- reg_hits(q, points)
  at <- integer()
  if (nrow(ht)) at <- unique(ht$query[
    links$event_id[use[ht$query]] == tx$event_id[ht$subject]])
  hr <- reg_hits(q, r)
  inside <- integer()
  if (nrow(hr) && nrow(er)) {
    rr <- data.table(query = hr$query, repeat_row = hr$subject,
      event_id = links$event_id[use[hr$query]])
    inside <- unique(merge(rr, er, by = c("event_id", "repeat_row"))$query)
  }
  links[use, `:=`(variant_at_event_TSS = seq_along(use) %in% at,
    variant_overlaps_TSS_TE = seq_along(use) %in% inside)]
  links
}
reg_context_worker <- function(cfg, z, root) {
  if (cfg$species != "human" || cfg$profile != "human_trio")
    stop("human-context requires a human trio profile")
  if (!reg_valid(z, "te-tss")) stop("Complete current te-tss analysis required")
  qspec <- reg_context_manifest(cfg, z, "qtl_tables")
  beds <- reg_context_manifest(cfg, z, "beds")
  required <- reg_default(z$human_context$require_qtl, TRUE)
  if (required && (is.null(qspec) || !nrow(qspec)))
    stop("human-context blocked: configure local qtl_tables and a verified phenotype_map; ",
      "public peaks cannot replace cis-QTL results")
  map_path <- if (!is.null(z$human_context$phenotype_map))
    reg_path(z$human_context$phenotype_map, cfg) else NULL
  if (!is.null(qspec) && is.null(map_path)) stop("QTL integration requires phenotype_map")
  # Peaks come only from the configured beds manifest (no automatic downloads).
  if (is.null(beds)) beds <- data.table(input_id = character(), path = character(),
    build = character(), evidence_type = character(), selection_used = logical())
  if (nrow(beds)) {
    require_columns(beds, c("evidence_type", "selection_used"))
    if (!is.logical(beds$selection_used) || anyNA(beds$selection_used))
      stop("Peak selection_used must be explicit TRUE/FALSE")
  }
  anno_path <- file.path(z$output_dir, "event_TE_annotation.tsv")
  tx_path <- file.path(z$output_dir, "transcript_TE_annotation.tsv")
  repeat_path <- file.path(z$output_dir, "repeat_intervals.tsv.gz")
  pairs_path <- file.path(z$output_dir, "te_tss_event_pairs.tsv")
  outputs <- file.path(z$output_dir, c("human_regulatory_context.tsv",
    "human_qtl_links.tsv.gz", "human_peak_overlaps.tsv.gz",
    "human_mapping_coverage.tsv", "human_resource_coverage.tsv",
    "human_context_provenance.tsv"))
  if (!reg_begin(z, "human-context", outputs)) return(invisible(TRUE))
  events <- reg_read(anno_path); tx <- reg_read(tx_path)
  aliases <- if (!is.null(z$contig_aliases)) reg_read(z$contig_aliases) else NULL
  target <- fread(file = cfg$paths$genome_fai, header = FALSE, select = 1L)[[1]]
  peak <- reg_context_overlap(tx, beds, target, aliases)
  links <- data.table(event_id = character(), kind = character(),
    link_level = character(), pip_all = numeric(), variant_overlaps_TSS_TE = logical())
  coverage <- data.table(input_id = character(), modality = character(),
    phenotype_id = character(), mapping_status = character())
  resources <- copy(peak$coverage)[, availability := "assessed_peak_resource"]
  if (!is.null(qspec) && nrow(qspec)) {
    q <- rbindlist(lapply(seq_len(nrow(qspec)), function(i) reg_context_qtl(qspec[i])),
      fill = TRUE)
    mapped <- reg_context_map(q, reg_read(map_path), events)
    links <- reg_context_variant_links(mapped$links, tx, reg_read(repeat_path),
      target, aliases)
    coverage <- mapped$coverage
    resources <- rbindlist(list(resources, q[, .(evidence_type = paste0("cis_", kind[1]),
      availability = "assessed_QTL_resource_not_all_events_measured",
      n_source_rows = .N, n_phenotypes = uniqueN(paste(modality, phenotype_id))),
      by = input_id]), fill = TRUE)
  } else {
    resources <- rbindlist(list(resources, data.table(input_id = "qtl_bundle",
      evidence_type = "cis_QTL", availability = "not_provided_explicit_chromatin_only")),
      fill = TRUE)
  }
  summary <- copy(events)
  pairs <- reg_read(pairs_path)
  if (nrow(pairs)) {
    representative <- pairs[!duplicated(gene_id)]
    require_columns(representative, c("event_te", "event_nonTE", "status"))
    keep <- intersect(c("gene_id", "event_te", "event_nonTE", "status", "n_units",
      "estimate_te_log2", "estimate_nonTE_log2", "magnitude_difference_log2",
      "signed_difference_log2", "signed_q_value", "magnitude_q_value"),
      names(representative))
    rep <- representative[, ..keep]
    setnames(rep, setdiff(names(rep), "gene_id"),
      paste0("representative_pair_", setdiff(names(rep), "gene_id")))
    summary <- merge(summary, rep, by = "gene_id", all.x = TRUE)
    summary[, representative_pair_event_role := fifelse(
      event_id == representative_pair_event_te, "TE",
      fifelse(event_id == representative_pair_event_nonTE, "nonTE", "other_event"))]
  }
  summary[, `:=`(qtl_assessed = !is.null(qspec) && nrow(qspec) > 0L,
    chromatin_assessed = nrow(beds) > 0L,
    n_event_qtl_links = if (!is.null(qspec)) 0L else NA_integer_,
    n_gene_context_links = if (!is.null(qspec)) 0L else NA_integer_,
    n_peak_overlaps = if (nrow(beds)) 0L else NA_integer_,
    max_event_link_pip = NA_real_, n_qtl_variants_in_TSS_TE = if (!is.null(qspec)) 0L else NA_integer_)]
  if (nrow(links)) {
    a <- links[, .(n_event_qtl_links = sum(link_level == "verified_phenotype_event"),
      n_gene_context_links = sum(link_level != "verified_phenotype_event"),
      max_event_link_pip = if (any(link_level == "verified_phenotype_event" &
          is.finite(pip_all))) max(pip_all[link_level == "verified_phenotype_event"],
            na.rm = TRUE) else NA_real_,
      n_qtl_variants_in_TSS_TE = uniqueN(paste(input_id, variant_id)[
        variant_overlaps_TSS_TE %in% TRUE])), by = event_id]
    summary[a, on = "event_id", `:=`(n_event_qtl_links = i.n_event_qtl_links,
      n_gene_context_links = i.n_gene_context_links,
      max_event_link_pip = i.max_event_link_pip,
      n_qtl_variants_in_TSS_TE = i.n_qtl_variants_in_TSS_TE)]
  }
  if (nrow(peak$hits)) {
    a <- peak$hits[, .(n_peak_overlaps = .N), by = event_id]
    summary[a, on = "event_id", n_peak_overlaps := i.n_peak_overlaps]
  }
  summary[, interpretation := "cis_and_chromatin_context_not_POE_QTL_colocalization_or_causal_TE"]
  summary[, event_QTL_mapping_status := fifelse(!qtl_assessed, "not_assessed",
    fifelse(n_event_qtl_links > 0, "verified_event_mapping",
      "no_verified_event_mapping_not_a_negative_QTL_result"))]
  objects <- list(summary, links, peak$hits, coverage, resources,
    data.table(receipt_policy = "immutable_inputs_output_presence",
      version = z$version, build = z$build,
      overlap_unit = "exact_transcript_TSS_or_same_TE_insertion",
      candidate_filter = "none_no_PIP_or_POE_significance_cutoff",
      qtl_required = required,
      caveat = "shared_GUSTO_or_annotation_selection_not_independent_replication"))
  for (i in seq_along(outputs)) reg_write(objects[[i]], outputs[i])
  reg_finish(z, "human-context", outputs)
}
