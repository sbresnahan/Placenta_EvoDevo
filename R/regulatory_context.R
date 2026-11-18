# cis-QTL, chromatin-peak and ppBMI GxE context helpers for the baseline context stage.
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
  links[, `:=`(variant_at_event_TSS = NA, variant_overlaps_TSS_TE = NA,
    variant_overlaps_promoter_TE = NA, variant_overlaps_proximal_only_TE = NA)]
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
  if (all(c("promoter_start", "promoter_end") %in% names(tx))) {
    wp <- tx[, .(chrom, start = promoter_start, end = promoter_end)]
    wh <- reg_hits(wp, r)
    wr <- if (nrow(wh)) unique(data.table(event_id = tx$event_id[wh$query],
      repeat_row = wh$subject)) else data.table(event_id = character(), repeat_row = integer())
    proximal <- if (nrow(wr)) fsetdiff(wr, er) else wr
    matched <- function(table) {
      if (!nrow(hr) || !nrow(table)) return(integer())
      rr <- data.table(query = hr$query, repeat_row = hr$subject,
        event_id = links$event_id[use[hr$query]])
      unique(merge(rr, table, by = c("event_id", "repeat_row"))$query)
    }
    links[use, `:=`(variant_overlaps_promoter_TE = seq_along(use) %in% matched(wr),
      variant_overlaps_proximal_only_TE = seq_along(use) %in% matched(proximal))]
  }
  links
}
