suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})
gtf_attribute <- function(x, key) {
  pattern <- paste0('(?:^|;)[[:space:]]*', key, ' "([^"]*)"')
  ans <- sub(paste0(".*?", pattern, ".*"), "\\1", x, perl = TRUE)
  ans[!grepl(pattern, x, perl = TRUE)] <- NA_character_
  ans
}
read_gtf <- function(path) {
  lead <- readLines(path, n = 100L, warn = FALSE)
  skip <- which(!grepl("^#", lead) & nzchar(lead))[1] - 1L
  if (is.na(skip)) stop("No GTF records within first 100 lines")
  z <- fread(path, skip = skip, header = FALSE, sep = "\t",
             quote = "", col.names = c("chrom", "source", "type", "start",
               "end", "score", "strand", "frame", "attributes"))
  if (ncol(z) != 9L) stop("GTF must have nine columns")
  for (a in c("gene_id", "gene_name", "transcript_id", "exon_number"))
    z[, (a) := gtf_attribute(attributes, a)]
  z[, exon_number := suppressWarnings(as.integer(exon_number))]
  z
}
audit_gtf <- function(z, chrom_sizes = NULL) {
  ex <- copy(z[type == "exon"])
  if (!nrow(ex) || anyNA(ex[, .(chrom, start, end, strand, gene_id, transcript_id)]))
    stop("Missing exons or required exon attributes")
  if (any(ex$start < 1L | ex$end < ex$start) ||
      any(!ex$strand %in% c("+", "-"))) stop("Invalid exon coordinates/strand")
  key <- c("transcript_id", "chrom", "start", "end")
  if (anyDuplicated(ex[, ..key])) stop("Repeated exon coordinates within a transcript")
  bad <- ex[, .(n_gene = uniqueN(gene_id), n_chrom = uniqueN(chrom),
                 n_strand = uniqueN(strand)), by = transcript_id][
                   n_gene != 1 | n_chrom != 1 | n_strand != 1]
  if (nrow(bad)) stop("Conflicting transcript gene/contig/strand mappings")
  setorder(ex, transcript_id, start, end)
  overlaps <- ex[, any(start[-1] <= head(end, -1)), by = transcript_id]
  if (any(overlaps$V1)) stop("Overlapping exons within a transcript")
  ex[, tx_exon_order := if (strand[1] == "+") seq_len(.N) else rev(seq_len(.N)),
     by = transcript_id]
  tx <- ex[, {
    k <- order(tx_exon_order)
    list(gene_id = gene_id[1], gene_name = gene_name[1], chrom = chrom[1],
         strand = strand[1], start = min(start), end = max(end),
         n_exons = .N, spliced_length = sum(end - start + 1),
         tss = if (strand[1] == "+") start[k[1]] else end[k[1]],
         first_start = start[k[1]], first_end = end[k[1]],
         first_donor = if (.N > 1) {
           if (strand[1] == "+") end[k[1]] else start[k[1]]
         } else NA_integer_,
         first_acceptor = if (.N > 1) {
           if (strand[1] == "+") start[k[2]] else end[k[2]]
         } else NA_integer_)
  }, by = transcript_id]
  declared <- z[type == "transcript", .(transcript_id, declared_start = start,
                                         declared_end = end)]
  assert_unique(declared, "transcript_id", "transcript feature")
  check <- merge(tx, declared, by = "transcript_id", all.x = TRUE)
  if (anyNA(check$declared_start) ||
      any(check$start != check$declared_start | check$end != check$declared_end))
    stop("Transcript feature bounds do not match exons")
  if (!is.null(chrom_sizes)) {
    cs <- fread(chrom_sizes, header = FALSE, col.names = c("chrom", "length"))
    cc <- merge(tx, cs, by = "chrom", all.x = TRUE)
    if (anyNA(cc$length) || any(cc$end > cc$length))
      stop("Annotation contigs/bounds incompatible with supplied genome sizes")
  }
  audit <- data.table(
    metric = c("gtf_rows", "genes_with_exons", "transcripts", "exons",
      "negative_strand_transcripts", "exons_with_nontranscriptional_number",
      "transcripts_with_nontranscriptional_numbers", "unspliced_transcripts"),
    value = c(nrow(z), uniqueN(ex$gene_id), nrow(tx), nrow(ex),
      sum(tx$strand == "-"),
      sum(ex$exon_number != ex$tx_exon_order, na.rm = TRUE),
      uniqueN(ex[!is.na(exon_number) & exon_number != tx_exon_order]$transcript_id),
      sum(tx$n_exons == 1L)))
  list(exons = ex, transcripts = tx, audit = audit)
}
group_promoters <- function(tx, window = 50L) {
  x <- copy(tx[n_exons > 1])
  setorder(x, gene_id, chrom, strand, first_donor, tss, transcript_id)
  x[, anchor := {
    out <- integer(.N); a <- tss[1]
    for (i in seq_len(.N)) {
      if (tss[i] - a > window) a <- tss[i]
      out[i] <- a
    }
    out
  }, by = .(gene_id, chrom, strand, first_donor)]
  x[, event_id := paste(gene_id, chrom, strand, first_donor, anchor, sep = "|")]
  x
}
annotate_repeats <- function(tx_groups, rmsk_path) {
  # UCSC rmsk: bin,swScore,milliDiv,milliDel,milliIns,genoName,genoStart,
  # genoEnd,genoLeft,strand,repName,repClass,repFamily,repStart,repEnd,repLeft,id
  r <- fread(cmd = paste("gzip -dc", shQuote(rmsk_path)), header = FALSE,
             select = c(6, 7, 8, 10, 11, 12, 13),
             col.names = c("chrom", "start0", "end", "repeat_strand",
                           "repeat_name", "repeat_class", "repeat_family"))
  r <- r[chrom %in% tx_groups$chrom]
  gr <- GRanges(r$chrom, IRanges(r$start0 + 1L, r$end))
  tx <- copy(tx_groups)
  pt <- GRanges(tx$chrom, IRanges(tx$tss, tx$tss))
  hit <- findOverlaps(pt, gr, ignore.strand = TRUE)
  rr <- cbind(data.table(i = queryHits(hit)), r[subjectHits(hit)])
  te_classes <- c("DNA", "LINE", "SINE", "LTR", "RC", "Retroposon")
  rr[, is_te := sub("/.*", "", repeat_class) %in% te_classes]
  rr[, is_unknown := repeat_class %in% c("Unknown", "Other")]
  lab <- rr[, .(tss_repeat_overlap = TRUE, tss_te = any(is_te),
    tss_unknown_repeat = any(is_unknown),
    te_class = paste(sort(unique(repeat_class[is_te])), collapse = ";"),
    te_family = paste(sort(unique(repeat_family[is_te])), collapse = ";"),
    te_name = paste(sort(unique(repeat_name[is_te])), collapse = ";"),
    te_orientation = paste(sort(unique(ifelse(repeat_strand[is_te] ==
                        tx$strand[i[1]], "sense", "antisense"))), collapse = ";")),
    by = i]
  tx[, i := .I]
  tx <- merge(tx, lab, by = "i", all.x = TRUE, sort = FALSE)
  tx[is.na(tss_repeat_overlap), tss_repeat_overlap := FALSE]
  tx[is.na(tss_te), tss_te := FALSE]
  tx[is.na(tss_unknown_repeat), tss_unknown_repeat := FALSE]
  lo <- ifelse(tx$strand == "+", pmax(1L, tx$tss - 500L), pmax(1L, tx$tss - 100L))
  hi <- ifelse(tx$strand == "+", tx$tss + 100L, tx$tss + 500L)
  win <- findOverlaps(GRanges(tx$chrom, IRanges(lo, hi)), gr, ignore.strand = TRUE)
  teidx <- which(sub("/.*", "", r$repeat_class) %in% te_classes)
  tx[, window_te := i %in% unique(queryHits(win)[subjectHits(win) %in% teidx])]
  # First-exon overlap is distinct from TSS overlap.
  fe <- findOverlaps(GRanges(tx$chrom, IRanges(tx$first_start, tx$first_end)),
                    gr, ignore.strand = TRUE)
  tx[, first_exon_te := i %in% unique(queryHits(fe)[subjectHits(fe) %in% teidx])]
  setorder(tx, transcript_id)
  tx[, i := NULL]
  tx
}
event_catalog <- function(tx, exons) {
  events <- tx[, .(
    gene_id = gene_id[1], gene_name = gene_name[1], chrom = chrom[1],
    strand = strand[1], first_donor = first_donor[1],
    tss_min = min(tss), tss_max = max(tss), n_transcripts = .N,
    transcript_ids = paste(sort(transcript_id), collapse = ";"),
    tss_te_label = if (any(tss_unknown_repeat)) "unknown_repeat" else
                    if (all(tss_te)) "TE" else if (!any(tss_te)) "nonTE" else "mixed",
    te_class = paste(sort(unique(te_class[!is.na(te_class) & nzchar(te_class)])), collapse = ";"),
    te_family = paste(sort(unique(te_family[!is.na(te_family) & nzchar(te_family)])), collapse = ";"),
    window_te = any(window_te), first_exon_te = any(first_exon_te),
    first_start = min(first_start), first_end = max(first_end)
  ), by = event_id]
  map <- tx[, .(transcript_id, event_id)]
  ee <- merge(exons, map, by = "transcript_id", all.x = TRUE)
  by_gene <- split(seq_len(nrow(ee)), ee$gene_id)
  # Conservative coordinate identifiability: subtract every exon from competing
  # promoter groups in the same gene and every exon from other genes.
  # Reduce per gene (not across genes), preserving gene identities.
  raw <- unique(exons[, .(gene_id, chrom, start, end)])
  reduced <- raw[, {
    q <- reduce(IRanges(start, end))
    list(start = start(q), end = end(q))
  }, by = .(gene_id, chrom)]
  allg <- GRanges(reduced$chrom, IRanges(reduced$start, reduced$end))
  fg <- GRanges(events$chrom, IRanges(events$first_start, events$first_end))
  hits <- findOverlaps(fg, allg, ignore.strand = TRUE)
  cross <- split(subjectHits(hits), queryHits(hits))
  events[, n_groups_gene := .N, by = gene_id]
  result <- vector("list", nrow(events))
  for (i in seq_len(nrow(events))) {
    e <- events[i]
    competing <- ee[by_gene[[e$gene_id]]]
    competing <- competing[is.na(event_id) | event_id != e$event_id]
    other <- cross[[as.character(i)]]
    if (length(other)) other <- other[reduced$gene_id[other] != e$gene_id]
    mask <- IRanges(competing$start, competing$end)
    if (length(other)) mask <- c(mask, IRanges(reduced$start[other], reduced$end[other]))
    unique_part <- setdiff(IRanges(e$first_start, e$first_end), reduce(mask))
    result[[i]] <- data.table(event_id = e$event_id,
      distinguishing_bp = sum(width(unique_part)),
      distinguishing_intervals = paste(paste(start(unique_part), end(unique_part), sep = "-"),
                                        collapse = ";"))
  }
  merge(events, rbindlist(result), by = "event_id", sort = FALSE)
}
select_structural_pairs <- function(events) {
  p <- events[tss_te_label %in% c("TE", "nonTE") & distinguishing_bp > 0]
  splitp <- split(p, p$gene_id)
  ans <- rbindlist(lapply(splitp, function(x) {
    if (nrow(x) < 2L) return(NULL)
    comb <- t(combn(x$event_id, 2L))
    a <- x[match(comb[, 1], event_id)]
    b <- x[match(comb[, 2], event_id)]
    kind <- ifelse(a$tss_te_label != b$tss_te_label, "TE_nonTE",
                   ifelse(a$tss_te_label == "nonTE", "nonTE_nonTE", "TE_TE"))
    cand <- data.table(gene_id = x$gene_id[1], gene_name = x$gene_name[1],
      event_a = a$event_id, event_b = b$event_id, pair_type = kind,
      min_distinguishing_bp = pmin(a$distinguishing_bp, b$distinguishing_bp))
    # A gene can contribute to only one exposure class; no self-matched controls.
    cand <- cand[pair_type != "TE_TE"]
    if (!nrow(cand)) return(NULL)
    cand[, preference := as.integer(pair_type == "TE_nonTE")]
    setorder(cand, -preference, -min_distinguishing_bp, event_a, event_b)
    cand[1, !c("preference")]
  }), fill = TRUE)
  if (!ncol(ans)) ans <- data.table(gene_id=character(),gene_name=character(),
    event_a=character(),event_b=character(),pair_type=character(),
    min_distinguishing_bp=integer())
  ans
}
