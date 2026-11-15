# All per-person inputs and outputs of these functions are private.
parse_gt <- function(gt) {
  gt <- sub(":.*", "", gt)
  if (is.na(gt) || !grepl("^[01][/|][01]$", gt)) return(c(NA_integer_, NA_integer_))
  as.integer(strsplit(gt, "[/|]")[[1]])
}
transmitted_maternal_allele <- function(child, mother, father) {
  c <- parse_gt(child); m <- parse_gt(mother); p <- parse_gt(father)
  if (anyNA(c(c, m, p))) return(list(maternal = NA_integer_, status = "missing_or_nonbiallelic"))
  possible <- expand.grid(m = unique(m), p = unique(p))
  possible <- possible[vapply(seq_len(nrow(possible)), function(i)
    identical(sort(as.integer(unlist(possible[i, ], use.names = FALSE))), sort(c)), logical(1)), , drop = FALSE]
  if (!nrow(possible)) return(list(maternal = NA_integer_, status = "mendelian_conflict"))
  if (length(unique(c)) != 2L) return(list(maternal = NA_integer_, status = "child_not_heterozygous"))
  if (length(unique(possible$m)) != 1L)
    return(list(maternal = NA_integer_, status = "ambiguous_transmission"))
  list(maternal = as.integer(possible$m[1]), status = "typed_trio")
}
read_ornament_counts <- function(path, semantics) {
  if (!identical(semantics, "ref_alt"))
    stop("Explicit ref_alt output semantics required; L/R is not assumed")
  x <- data.table::fread(path, header = FALSE, na.strings = c("NA", "."))
  if (ncol(x) != 8L) stop("Unrecognized Ornament count schema; expected 8 columns")
  data.table::setnames(x, c("transcript_id", "tx_pos", "tx_ref", "tx_alt",
                            "unused", "ref_count", "alt_count", "other_count"))
  assert_unique(x, c("transcript_id", "tx_pos"), "Ornaments count rows")
  if (any(!is.finite(x$ref_count)) || any(!is.finite(x$alt_count)) ||
      any(!is.finite(x$other_count)) ||
      any(x$ref_count < 0 | x$alt_count < 0 | x$other_count < 0) ||
      anyNA(x$tx_pos) || any(x$tx_pos < 1 | x$tx_pos != floor(x$tx_pos)))
    stop("Invalid expected counts/positions")
  x[, output_semantics := "ref_alt"]
  x
}
map_tx_loci <- function(counts, exons) {
  e <- data.table::copy(exons)
  data.table::setorder(e, transcript_id, tx_exon_order)
  e[, tx_end := cumsum(end - start + 1L), by = transcript_id]
  e[, tx_start := tx_end - (end - start)]
  c <- data.table::copy(counts)
  c[, count_row := .I]
  # Restrict to SNPs: indel reverse complementation requires renormalization.
  snp <- c[nchar(tx_ref) == 1L & nchar(tx_alt) == 1L &
             tx_ref %in% c("A","C","G","T") & tx_alt %in% c("A","C","G","T")]
  idx <- e[snp, on = .(transcript_id, tx_start <= tx_pos, tx_end >= tx_pos),
           nomatch = NULL, allow.cartesian = TRUE,
           .(count_row = i.count_row, chrom, exon_start = start, exon_end = end,
             exon_tx_start = x.tx_start, strand)]
  assert_unique(idx, "count_row", "Transcript position to exon mapping")
  out <- merge(c, idx, by = "count_row", all.x = TRUE)
  out[, genomic_pos := ifelse(strand == "+",
    exon_start + tx_pos - exon_tx_start, exon_end - (tx_pos - exon_tx_start))]
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  out[, ref := ifelse(strand == "-", unname(comp[tx_ref]), tx_ref)]
  out[, alt := ifelse(strand == "-", unname(comp[tx_alt]), tx_alt)]
  out[, mapping_status := ifelse(is.na(genomic_pos), "nonSNP_or_unmapped", "mapped_SNP")]
  out
}
orient_counts <- function(mapped_counts, trio_sites) {
  required <- c("chrom", "genomic_pos", "ref", "alt", "child_gt", "mother_gt", "father_gt")
  require_columns(trio_sites, required, "Typed trio sites")
  keys <- required[1:4]
  assert_unique(trio_sites, keys, "Typed trio sites")
  anchor <- data.table::copy(trio_sites)
  for (nm in c("chrom","ref","alt","child_gt","mother_gt","father_gt"))
    data.table::set(anchor,j=nm,value=as.character(anchor[[nm]]))
  anchor[,genomic_pos:=as.integer(genomic_pos)]
  lab <- lapply(seq_len(nrow(anchor)), function(i)
    transmitted_maternal_allele(anchor$child_gt[i], anchor$mother_gt[i], anchor$father_gt[i]))
  anchor[, maternal_index := vapply(lab, `[[`, integer(1), "maternal")]
  anchor[, parent_status := vapply(lab, `[[`, character(1), "status")]
  out <- merge(mapped_counts, anchor, by = keys, all.x = TRUE, sort = FALSE)
  out[is.na(parent_status), parent_status := "no_typed_anchor"]
  out[, maternal_count := ifelse(parent_status == "typed_trio",
    ifelse(maternal_index == 0L, ref_count, alt_count), NA_real_)]
  out[, paternal_count := ifelse(parent_status == "typed_trio",
    ifelse(maternal_index == 0L, alt_count, ref_count), NA_real_)]
  out[, ref_maternal := ifelse(parent_status == "typed_trio", maternal_index == 0L, NA)]
  out[, locus_id := paste(chrom, genomic_pos, ref, alt, sep = ":")]
  out
}
event_locus_counts <- function(parental_counts, tx_map, sample_id, family_id) {
  assert_unique(tx_map, "transcript_id", "Transcript/event map")
  c <- merge(parental_counts, tx_map[, .(transcript_id, event_id, gene_id)],
             by = "transcript_id", all.x = TRUE)
  if (anyNA(c$event_id)) stop("Counts contain transcripts absent from reference event map")
  out <- c[, .(maternal_count = if (all(!is.na(maternal_count))) sum(maternal_count) else NA_real_,
    paternal_count = if (all(!is.na(paternal_count))) sum(paternal_count) else NA_real_,
    ref_maternal = if (uniqueN(ref_maternal) == 1L) ref_maternal[1] else NA,
    n_allocated_transcripts = uniqueN(transcript_id),
    distinguishing_coordinate = if ("distinguishing_coordinate" %in% names(c))
      all(distinguishing_coordinate) else FALSE), by = .(gene_id, event_id, locus_id)]
  out[, `:=`(sample_id = sample_id, family_id = family_id,
              evidence_tier = "model_based_expected_counts")]
  out
}

mark_direct_coordinates <- function(mapped, unique_intervals, tx_map) {
  # A locus in distinguishing sequence is a coordinate support flag, not proof
  # of unique k-mers, mapping accuracy, or a jointly informative raw molecule.
  m <- merge(mapped, tx_map[, .(transcript_id, event_id)], by = "transcript_id", all.x = TRUE)
  require_columns(unique_intervals, c("event_id", "chrom", "start", "end"))
  m[, count_index := .I]
  hits <- unique_intervals[m[!is.na(genomic_pos)],
    on = .(event_id, chrom, start <= genomic_pos, end >= genomic_pos),
    nomatch = NULL, allow.cartesian = TRUE, .(i.count_index)]
  m[, distinguishing_coordinate := count_index %in% hits$i.count_index]
  m[, c("event_id", "count_index") := NULL]
  m
}
