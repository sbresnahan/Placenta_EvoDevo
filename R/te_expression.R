ex_gr <- function(d, start = "start", end = "end", strand = FALSE) {
  GenomicRanges::GRanges(d$chrom, IRanges::IRanges(d[[start]], d[[end]]),
    strand = if (strand) d$strand else "*")
}
ex_gene_spans <- function(annotation) {
  ex_need(annotation, c("gene_id", "chrom", "start", "end", "transcript_id", "tss"))
  z <- annotation[, .(chroms = data.table::uniqueN(chrom), chrom = chrom[1],
    start = min(start), end = max(end)), by = gene_id]
  if (any(z$chroms != 1L)) stop("A gene maps to multiple contigs")
  z[, chroms := NULL]; z
}
ex_match_targets <- function(annotation, loci) {
  ex_need(annotation, c("transcript_id", "gene_id", "chrom", "tss", "tss_te_label"))
  active <- "promoter_te_label" %in% names(annotation)
  tx <- if (active) annotation[promoter_te_label == "TE"] else annotation[tss_te_label == "TE"]
  if (!nrow(tx)) return(data.table::data.table(level = character(), feature_id = character(),
    gene_id = character(), family_id = character()))
  start <- if (active) tx$promoter_start else tx$tss
  end <- if (active) tx$promoter_end else tx$tss
  hits <- GenomicRanges::findOverlaps(GenomicRanges::GRanges(tx$chrom, IRanges::IRanges(start, end)), ex_gr(loci))
  if (!length(hits)) stop("No TE promoter-window overlaps; check build and contig names")
  q <- S4Vectors::queryHits(hits); h <- S4Vectors::subjectHits(hits)
  insertion <- data.table::data.table(level = "transcript", feature_id = tx$transcript_id[q],
    gene_id = tx$gene_id[q], family_id = loci$family_id[h], locus_id = loci$locus_id[h],
    insertion_chrom = loci$chrom[h], insertion_start = loci$start[h], insertion_end = loci$end[h],
    transcript_tss = tx$tss[q], promoter_start = start[q], promoter_end = end[q],
    distance = pmax(loci$start[h] - tx$tss[q], tx$tss[q] - loci$end[h], 0),
    te_tss_window_bp = if (active) tx$te_tss_window_bp[q] else 0L)
  summarize <- function(d) d[, .(insertion_ids = paste(sort(unique(locus_id)), collapse = ";"),
    n_insertions = data.table::uniqueN(locus_id), exact_TSS = any(distance == 0),
    proximal_only = all(distance > 0), min_distance_bp = min(distance),
    te_tss_window_bp = te_tss_window_bp[1]), by = .(level, feature_id, gene_id, family_id)]
  genes <- data.table::copy(insertion)[, `:=`(level = "gene", feature_id = gene_id)]
  out <- data.table::rbindlist(list(summarize(insertion), summarize(genes)))
  insertion[, `:=`(exact_TSS = distance == 0, proximal_only = distance > 0)]
  attr(out, "insertion_links") <- insertion
  out
}
ex_catalogue <- function(repeats, cfg) {
  ex_need(repeats, c("chrom", "start", "end", "repeat_strand", "repeat_id", "repeat_name",
    "repeat_family", "repeat_class", "is_te", "is_unknown", "source"))
  x <- data.table::copy(repeats[is_te & !is_unknown])
  x[, locus_id := paste(cfg$species, chrom, start, end, repeat_strand, repeat_id, sep = "|")]
  d <- unique(x[, .(locus_id, species = cfg$species, build = cfg$reference$build,
    chrom, start, end, strand = repeat_strand, family_id = repeat_family,
    repeat_name, repeat_class, source)])
  ex_validate_loci(d, cfg); d
}
ex_write_te_gtf <- function(loci, path) {
  # Disjoint pieces prevent duplicate segment ownership. Conflicting families
  # or strands are excluded; nested same-family pieces use the shortest locus.
  gr <- ex_gr(loci)
  pieces <- GenomicRanges::disjoin(gr, ignore.strand = TRUE)
  h <- GenomicRanges::findOverlaps(pieces, gr)
  map <- data.table::data.table(piece = S4Vectors::queryHits(h), row = S4Vectors::subjectHits(h))
  map <- merge(map, loci[, .(row = .I, locus_id, family_id, strand, width = end - start + 1L)], by = "row")
  valid <- map[, .(families = data.table::uniqueN(family_id), strands = data.table::uniqueN(strand)), by = piece]
  map <- map[piece %in% valid[families == 1L & strands == 1L]$piece]
  data.table::setorder(map, piece, width, locus_id)
  owner <- map[!duplicated(piece)]
  if (!nrow(owner)) stop("No unambiguous TE annotation pieces")
  pp <- pieces[owner$piece]
  safe_id <- function(x) {
    for (bad in c('"', "\\", ";", "\n", "\t")) x <- gsub(bad, "_", x, fixed = TRUE)
    x
  }
  gtf <- data.table::data.table(chrom = as.character(GenomicRanges::seqnames(pp)),
    source = "exposure_nb_v1", feature = "exon", start = GenomicRanges::start(pp),
    end = GenomicRanges::end(pp), score = ".", strand = owner$strand, phase = ".",
    attributes = paste0('locus "', safe_id(owner$locus_id), '"; transcript_id "',
      safe_id(owner$locus_id), '"; gene_id "', safe_id(owner$locus_id), '";'))
  if (!identical(safe_id(owner$locus_id), owner$locus_id)) stop("Unsafe locus IDs for GTF")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(gtf, path, sep = "\t", quote = FALSE, col.names = FALSE)
  data.table::data.table(n_loci = nrow(loci), n_segments = nrow(gtf),
    ambiguous_segments_excluded = sum(valid$families > 1 | valid$strands > 1),
    n_loci_represented = data.table::uniqueN(owner$locus_id))
}
ex_predictors <- function(counts, qc, loci, targets, spans, s) {
  samples <- qc$sample_id
  if (!length(samples) || !nrow(targets)) return(data.table::data.table())
  if (isTRUE(s$family_only)) {
    base <- data.table::copy(counts)
    modes <- "inclusive_family_only"
  } else {
    base <- merge(counts, loci[, .(locus_id, family_id, chrom, start, end)], by = "locus_id")
    intergenic <- !seq_len(nrow(loci)) %in% S4Vectors::queryHits(GenomicRanges::findOverlaps(ex_gr(loci), ex_gr(spans)))
    base[, intergenic := locus_id %in% loci$locus_id[intergenic]]
    modes <- c("target_excluded", "inclusive", "intergenic")
    if ("unique_count" %in% names(base)) modes <- c(modes, "unique_target_excluded")
  }
  keys <- unique(targets[, .(gene_id, family_id)])
  out <- lapply(seq_len(nrow(keys)), function(i) {
    gene <- keys$gene_id[i]; fam <- keys$family_id[i]
    span <- spans[gene_id == gene]
    if (nrow(span) != 1L) stop("Missing/ambiguous target gene span")
    z <- base[family_id == fam]
    if (!s$family_only) z[, local := chrom == span$chrom & start <= span$end + s$local_flank &
                                          end >= max(1, span$start - s$local_flank)]
    data.table::rbindlist(lapply(modes, function(mode) {
      zz <- switch(mode, target_excluded = z[local == FALSE], inclusive = z,
        intergenic = z[intergenic == TRUE], unique_target_excluded = z[local == FALSE], inclusive_family_only = z)
      nm <- if (mode == "unique_target_excluded") "unique_count" else "count"
      sums <- zz[, .(family_count = sum(get(nm))), by = sample_id]
      v <- sums$family_count[match(samples, sums$sample_id)]; v[is.na(v)] <- 0
      cpm <- 1e6 * v / qc$mapped_fragments
      logv <- log2(cpm + .5); center <- mean(logv); scale <- stats::sd(logv)
      n_expr <- sum(cpm >= s$min_cpm)
      status <- if (n_expr < max(s$min_expressed, ceiling(s$min_prevalence * length(samples))))
        "insufficient_TE_expression" else if (!is.finite(scale) || scale <= 0) "constant_TE_predictor" else "eligible"
      data.table::data.table(sample_id = samples, gene_id = gene, family_id = fam, predictor = mode,
        family_count = v, cpm = cpm, log_cpm = logv,
        TE_abundance = if (status == "eligible") (logv - center) / scale else NA_real_,
        center = center, scale = scale, n_expressed = n_expr, predictor_status = status)
    }), fill = TRUE)
  })
  data.table::rbindlist(out, fill = TRUE)
}
ex_te_plan <- function(cfg, s) {
  dir.create(ex_path(s, "quantification_plan"), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(s$te_loci)) {
    if (is.null(s$repeat_intervals)) stop("Supply te_loci or repeat_intervals; fish catalogue must be complete")
    loci <- ex_catalogue(ex_read(s$repeat_intervals), cfg)
    ex_write(loci, ex_path(s, "quantification_plan", "te_loci.tsv"))
  } else {loci <- ex_read(s$te_loci); ex_validate_loci(loci, cfg)}
  gtf <- ex_path(s, "quantification_plan", "te_loci.gtf")
  stats <- ex_write_te_gtf(loci, gtf)
  ex_write(stats, ex_path(s, "quantification_plan", "annotation_QC.tsv"))
  lib <- ex_read(s$libraries)
  ex_need(lib, c("sample_id", "bam", "fastq1", "fastq2", "layout", "strandedness",
    "star_index", "telescope_report", "mapped_fragments"), "quantification libraries")
  ex_unique(lib, "sample_id")
  if (any(!grepl("^[A-Za-z0-9_.-]+$", lib$sample_id))) stop("Unsafe sample ID in job scripts")
  if (any(!lib$layout %in% c("PE", "SE")) || any(!lib$strandedness %in% c("unstranded", "forward", "reverse")))
    stop("Declare SE/PE layout and strandedness")
  scripts <- lapply(seq_len(nrow(lib)), function(i) {
    l <- lib[i]; sample <- l$sample_id
    work <- ex_path(s, "private", "quantification", sample)
    bam <- l$bam
    header <- c("#!/bin/bash", "set -euo pipefail", "umask 077",
      paste0("WORK=", shQuote(work)), 'mkdir -p "$WORK"',
      paste0('[[ "$(telescope --version 2>&1)" == *', shQuote(s$quantifier_version), '* ]] || { echo "Telescope version mismatch" >&2; exit 1; }'),
      paste0('[[ "$(git -C "${EX_TELESCOPE_SOURCE:?Set pinned Telescope source checkout}" rev-parse HEAD)" == ',
        shQuote(s$quantifier_commit), ' ]] || { echo "Telescope source commit mismatch" >&2; exit 1; }'),
      "# Inspect the alignment and annotation QC before accepting counts.")
    if (is.na(bam) || !nzchar(bam)) {
      if (!ex_text(l$fastq1) || !ex_text(l$star_index)) stop("Need a BAM or FASTQ plus STAR index")
      if (l$layout == "PE" && !ex_text(l$fastq2)) stop("Paired library lacks FASTQ2")
      bam <- file.path(work, "Aligned.out.bam")
      reads <- paste(shQuote(l$fastq1), if (l$layout == "PE") shQuote(l$fastq2) else "")
      header <- c(header,
        paste0('[[ "$(STAR --version)" == ', shQuote(s$star_version), ' ]] || { echo "STAR version mismatch" >&2; exit 1; }'),
        paste("STAR --runThreadN", '${EX_ALIGN_CPUS:?Set EX_ALIGN_CPUS after pilot}',
          "--genomeDir", shQuote(l$star_index), "--readFilesIn", reads,
          if (grepl("\\.gz$", l$fastq1)) "--readFilesCommand zcat" else "",
          "--outFilterMultimapNmax", s$multimap_limit, "--winAnchorMultimapNmax", s$multimap_limit,
          "--outSAMtype BAM Unsorted --outSAMattributes NH HI AS nM",
          "--outFileNamePrefix", shQuote(paste0(work, "/"))))
    }
    stranded <- switch(l$strandedness, unstranded = "", forward = if (l$layout == "PE") " --stranded_mode FR" else " --stranded_mode F",
      reverse = if (l$layout == "PE") " --stranded_mode RF" else " --stranded_mode R")
    script <- c(header,
      paste("samtools quickcheck", shQuote(bam)),
      paste('Rscript "${EX_MODULE_ROOT:?Set module checkout root}/bin/check_exposure_bam.R"',
        shQuote(bam), s$multimap_limit, shQuote(file.path(work, "BAM_probe.tsv"))),
      paste("samtools collate -o", shQuote(file.path(work, "collated.bam")), shQuote(bam)),
      paste0("telescope assign ", shQuote(file.path(work, "collated.bam")), " ", shQuote(gtf),
        " --attribute locus --reassign_mode average --overlap_threshold ", s$annotation_min_overlap,
        stranded, " --outdir ", shQuote(work), " --exp_tag ", shQuote(sample)),
      "# Import the report using an explicit QC manifest after verifying BAM tags,",
      "# multimapping limits, mapped-fragment denominator, strand and declared annotation ID.")
    if ("unique_telescope_report" %in% names(l) && ex_text(l$unique_telescope_report)) {
      unique_bam <- file.path(work, "NH1.bam")
      script <- c(script,
        "# Strict genomic unique-mapper sensitivity; NOT Telescope annotation-unique_count.",
        paste("samtools view -h", shQuote(file.path(work, "collated.bam")),
          "| awk", shQuote('/^@/ || /\\tNH:i:1(\\t|$)/'),
          "| samtools view -b -o", shQuote(unique_bam), "-"),
        paste0("telescope assign ", shQuote(unique_bam), " ", shQuote(gtf),
          " --attribute locus --reassign_mode average --overlap_threshold ", s$annotation_min_overlap,
          stranded, " --outdir ", shQuote(work), " --exp_tag ", shQuote(paste0(sample, "_NH1"))))
    }
    dest <- ex_path(s, "quantification_plan", paste0(sample, ".sh"))
    writeLines(script, dest)
    data.table::data.table(sample_id = sample, script = dest, executed = FALSE,
      annotation_id = s$annotation_id %||ex% basename(gtf))
  })
  ex_write(data.table::rbindlist(scripts), ex_path(s, "quantification_plan", "jobs.tsv"))
  invisible(stats)
}
ex_bam_probe <- function(path, multimap_limit, yield_size = 10000L) {
  bf <- Rsamtools::BamFile(path, yieldSize = yield_size)
  on.exit(close(bf), add = TRUE); open(bf)
  d <- Rsamtools::scanBam(bf, param = Rsamtools::ScanBamParam(
    what = c("flag", "qname"), tag = c("NH", "HI", "AS", "nM")))[[1]]
  mapped <- bitwAnd(d$flag, 4L) == 0L
  if (!any(mapped)) stop("BAM probe contains no mapped alignments")
  for (tag in c("NH", "HI", "AS", "nM")) {
    v <- d$tag[[tag]]
    if (length(v) != length(mapped) || any(!is.finite(v[mapped]))) stop("BAM probe missing tag: ", tag)
  }
  nh <- d$tag$NH[mapped]
  if (any(nh < 1L | nh > multimap_limit)) stop("BAM NH outside declared multimapping limit")
  n_multi <- sum(nh > 1L)
  n_secondary <- sum(bitwAnd(d$flag[mapped], 256L) != 0L)
  if (n_multi > 0L && n_secondary == 0L) stop("Probe has multimappers but no secondary alignments")
  data.table::data.table(n_alignments_probed = length(mapped), n_mapped = sum(mapped),
    n_multimapping_alignments = n_multi, n_secondary = n_secondary, max_NH = max(nh),
    status = "probe_passed_full_library_review_required")
}
ex_import_telescope <- function(manifest, loci) {
  ex_need(manifest, c("sample_id", "telescope_report", "mapped_fragments"))
  ex_unique(manifest, "sample_id")
  data.table::rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
    x <- data.table::fread(manifest$telescope_report[i], skip = "transcript\t",
      na.strings = c("", "NA"), showProgress = FALSE)
    # Pinned strand-aware 1.0.3.1 report headers; no position-based guessing.
    ex_need(x, c("transcript", "final_count", "unique_count"), "Telescope report")
    x <- x[transcript != "__no_feature"]
    if (any(!x$transcript %in% loci$locus_id)) stop("Telescope report has unknown loci")
    out <- x[, .(sample_id = manifest$sample_id[i], locus_id = transcript,
      count = as.numeric(final_count), annotation_unique_count = as.numeric(unique_count))]
    if ("unique_telescope_report" %in% names(manifest) && ex_text(manifest$unique_telescope_report[i])) {
      u <- data.table::fread(manifest$unique_telescope_report[i], skip = "transcript\t", showProgress = FALSE)
      ex_need(u, c("transcript", "final_count"), "NH1 Telescope report")
      u <- u[transcript != "__no_feature"]
      ex_unique(u, "transcript", "NH1 report")
      if (any(!u$transcript %in% loci$locus_id)) stop("NH1 report has unknown loci")
      out <- merge(out, u[, .(locus_id = transcript, unique_count = as.numeric(final_count))],
        by = "locus_id", all = TRUE)
      out[, sample_id := manifest$sample_id[i]]
      out[is.na(count), count := 0]; out[is.na(unique_count), unique_count := 0]
    }
    out
  }), fill = TRUE)
}
