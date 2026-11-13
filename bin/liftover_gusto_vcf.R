#!/usr/bin/env Rscript
# Lift the GUSTO transcriptomic VCF from ESPRESSO transcripts to a different
# transcript assembly (HPLRv2.0) using gffcompare best matches.
#
# The two assemblies are independent: transcript IDs are never trusted, even
# when an ENST ID appears in both. Each ESPRESSO transcript is linked to its
# gffcompare best-matching target transcript (.tmap ref_id), restricted to
# accepted class codes. Each variant is then re-projected from its hg38 locus,
# taken from VCF ID (chrom:pos:ref:alt + vN, forward-strand alleles), onto the
# target transcript's exons. Nothing is carried over from the source CHROM/POS.
#
#   SNV/MNP  must lie inside one target exon.
#   Indel    (VCF-anchored, shared first base) the anchor through the base after
#            REF must lie inside one target exon, so no indel touches a target
#            splice junction. On minus-strand targets the anchor becomes the
#            base following REF on the genome, as in the source VCF.
#   REF      every lifted REF is checked against the target transcript FASTA.
#
# Sample columns (GT:DS:GP...) are copied verbatim: genotypes and Left/Right
# phase are genomic and do not depend on transcript strand. vN is renumbered
# per genomic locus over target transcripts. Run bcftools norm -c e afterward
# (see the runbook) to left-align indels in transcript space.
#
# --project-all also projects every locus onto EVERY target transcript whose
# exon contains it, in addition to the best-match pairs. Without it, target
# isoforms with no ESPRESSO best match stay unpersonalized; the summary reports
# how many transcript/locus pairs that leaves out.
#
# Usage:
#   bin/Rscript_sif bin/liftover_gusto_vcf.R \
#     --vcf SOURCE.vcf (uncompressed) --tmap ESP_vs_TARGET.tmap \
#     --target-gtf TARGET.gtf --target-fasta TARGET_transcripts.fa (indexed) \
#     --out-prefix OUT [--class-codes "=,c,k,m,n,j"] [--project-all]
f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
f <- gsub("~+~", " ", f, fixed = TRUE)
root <- normalizePath(file.path(dirname(f), ".."), mustWork = TRUE)
source(file.path(root, "R/library_paths.R"))
source(file.path(root, "R/common.R"))
source(file.path(root, "R/annotation.R"))
suppressPackageStartupMessages({
  library(Rsamtools); library(Biostrings); library(GenomicRanges)
})

# ---- arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opt <- list(class_codes = "=,c,k,m,n,j", project_all = FALSE)
i <- 1L
while (i <= length(args)) {
  k <- args[i]
  if (k == "--project-all") { opt$project_all <- TRUE; i <- i + 1L; next }
  if (!startsWith(k, "--") || i == length(args)) stop("Bad argument: ", k)
  opt[[gsub("-", "_", sub("^--", "", k))]] <- args[i + 1L]; i <- i + 2L
}
for (k in c("vcf", "tmap", "target_gtf", "target_fasta", "out_prefix"))
  if (is.null(opt[[k]])) stop("Missing --", gsub("_", "-", k))
for (p in c(opt$vcf, opt$tmap, opt$target_gtf, opt$target_fasta,
            paste0(opt$target_fasta, ".fai")))
  if (!file.exists(p)) stop("Missing input: ", p)
codes <- strsplit(opt$class_codes, ",", fixed = TRUE)[[1]]
out_vcf <- paste0(opt$out_prefix, ".unsorted.vcf")
if (file.exists(out_vcf)) stop("Refusing to overwrite: ", out_vcf)
dir.create(dirname(out_vcf), recursive = TRUE, showWarnings = FALSE)
key2 <- function(a, b) paste(a, b, sep = "\r")
rc <- function(x) {
  out <- rep("", length(x)); k <- nzchar(x)
  if (any(k)) out[k] <- as.character(reverseComplement(DNAStringSet(x[k])))
  out
}

# ---- source VCF ----------------------------------------------------------------
h <- readLines(opt$vcf, n = 500000L, warn = FALSE)
hi <- which(startsWith(h, "#CHROM"))[1]
if (is.na(hi)) stop("No #CHROM header line in the first 500000 lines")
meta <- h[seq_len(hi - 1L)]; colline <- h[hi]; rm(h)
v <- fread(opt$vcf, skip = hi - 1L, sep = "\t", header = TRUE,
           colClasses = "character", quote = "")
fixed <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
if (ncol(v) < 10L) stop("VCF has no sample columns")
setnames(v, 1:9, fixed)
samples <- names(v)[-(1:9)]
v[, gstr := do.call(paste, c(.SD, sep = "\t")), .SDcols = samples]
v[, (samples) := NULL]
v[, src_row := .I]
message("Source records: ", nrow(v), "; samples: ", length(samples))

v[, valid_id := grepl("^chr([0-9]+|X|Y|M|MT):[0-9]+:[ACGT]+:[ACGT]+v[0-9]+$", ID)]
v[, locus := sub("v[0-9]+$", "", ID)]
parts <- tstrsplit(v$locus, ":", fixed = TRUE, fill = NA)
while (length(parts) < 4L) parts[[length(parts) + 1L]] <- rep(NA_character_, nrow(v))
v[, `:=`(chrom = fifelse(valid_id, parts[[1]], NA_character_),
         gpos = fifelse(valid_id, suppressWarnings(as.integer(parts[[2]])), NA_integer_),
         gref = fifelse(valid_id, parts[[3]], NA_character_),
         galt = fifelse(valid_id, parts[[4]], NA_character_))]
rm(parts)
v[, status := fifelse(valid_id & !is.na(gpos) & gpos >= 1L, NA_character_,
                      "invalid_or_nonsequence_id")]
# Copies of one genomic locus on different source transcripts must carry identical genotypes.
v[is.na(status), n_geno := uniqueN(gstr), by = locus]
v[is.na(status) & n_geno > 1L, status := "conflicting_genotypes_across_copies"]
v[, `:=`(ref_len = nchar(gref), alt_len = nchar(galt))]
v[is.na(status), vtype := fifelse(ref_len == 1L & alt_len == 1L, "snv",
  fifelse(ref_len == alt_len, "mnp",
    fifelse(substr(gref, 1L, 1L) == substr(galt, 1L, 1L), "indel", "complex")))]
v[is.na(status) & vtype == "complex", status := "complex_allele"]
v[, `:=`(allele_ok = is.na(status), tgt_pos = NA_integer_)]

# ---- target annotation / FASTA ---------------------------------------------------
ta <- audit_gtf(read_gtf(opt$target_gtf))
ex <- copy(ta$exons)[, .(transcript_id, chrom, start, end, strand, tx_exon_order)]
setorder(ex, transcript_id, tx_exon_order)
ex[, offset := c(0L, head(cumsum(end - start + 1L), -1L)), by = transcript_id]
fai <- fread(paste0(opt$target_fasta, ".fai"), header = FALSE, select = 1:2,
             col.names = c("transcript_id", "length"))
if (!setequal(fai$transcript_id, ta$transcripts$transcript_id))
  stop("Target FASTA transcript IDs do not match the target GTF")
m_len <- merge(fai, ta$transcripts[, .(transcript_id, spliced_length)], by = "transcript_id")
if (any(m_len$length != m_len$spliced_length)) stop("Target FASTA lengths do not match the target GTF")
fa <- FaFile(opt$target_fasta); open(fa)
fetch <- function(tx, pos, w)
  toupper(as.character(scanFa(fa, GRanges(tx, IRanges(pos, width = w)))))

# ---- gffcompare best matches ----------------------------------------------------
tm <- fread(opt$tmap, sep = "\t", header = TRUE, colClasses = "character")
require_columns(tm, c("ref_id", "class_code", "qry_id"), "gffcompare tmap")
assert_unique(tm, "qry_id", "gffcompare tmap")
v <- merge(v, tm[, .(CHROM = qry_id, best_tx = ref_id, class_code)],
           by = "CHROM", all.x = TRUE, sort = FALSE)
setorder(v, src_row)
v[is.na(status) & is.na(class_code), status := "source_transcript_absent_from_tmap"]
v[is.na(status) & !class_code %in% codes, status := paste0("class_code_rejected:", class_code)]
v[is.na(status) & !best_tx %in% fai$transcript_id, status := "best_match_absent_from_target"]

# ---- candidate (locus, target transcript) pairs ----------------------------------
# All sequence-resolved loci are candidates for projection, whatever their best-match fate.
g <- unique(v[allele_ok == TRUE], by = "locus")[,
  .(locus, chrom, gpos, gref, galt, ref_len, alt_len, vtype, QUAL, FILTER, INFO, FORMAT, gstr)]
g[, span_end := fifelse(vtype == "indel", gpos + ref_len, gpos + ref_len - 1L)]
ov <- suppressWarnings(findOverlaps(GRanges(g$chrom, IRanges(g$gpos, g$span_end)),
  GRanges(ex$chrom, IRanges(ex$start, ex$end)), type = "within"))
exon_pairs <- unique(data.table(locus = g$locus[queryHits(ov)],
                                tgt_tx = ex$transcript_id[subjectHits(ov)]))
pairs <- unique(v[is.na(status), .(locus, tgt_tx = best_tx)])[, origin := "best_match"]
if (isTRUE(opt$project_all))
  pairs <- rbind(pairs, copy(exon_pairs)[, origin := "exon_projection"])
pairs <- pairs[, .(origin = paste(sort(unique(origin)), collapse = "+")), by = .(locus, tgt_tx)]

# ---- projection ------------------------------------------------------------------
m <- merge(pairs, g, by = "locus")
m <- merge(m, ex[, .(tgt_tx = transcript_id, e_chrom = chrom, start, end, strand, offset)],
           by = "tgt_tx", allow.cartesian = TRUE)
m <- m[e_chrom == chrom & start <= gpos & end >= span_end]
assert_unique(m, c("locus", "tgt_tx"), "Lifted locus/transcript pairs")
m[strand == "+", `:=`(tx_pos = offset + (gpos - start + 1L), tx_ref = gref, tx_alt = galt)]
m[strand == "-" & vtype != "indel",
  `:=`(tx_pos = offset + (end - (gpos + ref_len - 1L) + 1L), tx_ref = rc(gref), tx_alt = rc(galt))]
m[strand == "-" & vtype == "indel", tx_pos := offset + (end - (gpos + ref_len) + 1L)]
k <- m[, which(strand == "-" & vtype == "indel")]
if (length(k)) {
  anchor <- fetch(m$tgt_tx[k], m$tx_pos[k], 1L)
  m[k, `:=`(tx_ref = paste0(anchor, rc(substring(gref, 2L))),
            tx_alt = paste0(anchor, rc(substring(galt, 2L))))]
}
m[, target_seq := fetch(tgt_tx, tx_pos, nchar(tx_ref))]
m_ok <- m[target_seq == tx_ref]
m_ok[, n_dup := .N, by = .(tgt_tx, tx_pos, tx_ref, tx_alt)]
ok <- m_ok[n_dup == 1L]
setorder(ok, locus, tgt_tx)
ok[, ID := paste0(locus, "v", seq_len(.N) - 1L), by = locus]

pairs[, key := key2(locus, tgt_tx)]
pairs[, pair_status := "site_not_in_target_exon"]
pairs[key %in% key2(m$locus, m$tgt_tx), pair_status := "target_ref_mismatch"]
pairs[key %in% key2(m_ok$locus, m_ok$tgt_tx), pair_status := "duplicate_target_record"]
pairs[key %in% key2(ok$locus, ok$tgt_tx), pair_status := "lifted"]
v[is.na(status), status := pairs$pair_status[match(key2(locus, best_tx), pairs$key)]]
v[status == "lifted", tgt_pos := ok$tx_pos[match(key2(locus, best_tx), key2(ok$locus, ok$tgt_tx))]]

# ---- write VCF -----------------------------------------------------------------
keep_meta <- meta[!startsWith(meta, "##contig=") & !startsWith(meta, "##bcftools_")]
contigs <- paste0("##contig=<ID=", fai$transcript_id, ",length=", fai$length, ">")
prov <- paste0("##GUSTO_liftover=source:", opt$vcf, " tmap:", opt$tmap,
               " target_gtf:", opt$target_gtf, " target_fasta:", opt$target_fasta,
               " class_codes:", opt$class_codes, " project_all:", opt$project_all,
               " date:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
ok[, contig_order := match(tgt_tx, fai$transcript_id)]
setorder(ok, contig_order, tx_pos)
writeLines(c(keep_meta[1], contigs, keep_meta[-1], prov, colline), out_vcf)
fwrite(ok[, .(tgt_tx, tx_pos, ID, tx_ref, tx_alt, QUAL, FILTER, INFO, FORMAT, gstr)],
       out_vcf, append = TRUE, sep = "\t", quote = FALSE, col.names = FALSE)

# ---- audits ---------------------------------------------------------------------
write_tsv(v[, .(src_row, src_tx = CHROM, src_pos = POS, ID, locus, vtype, best_tx,
                class_code, status, tgt_pos)], paste0(opt$out_prefix, ".source_records.tsv"))
write_tsv(ok[, .(tgt_tx, tx_pos, ID, tx_ref, tx_alt, locus, vtype, strand, origin)],
          paste0(opt$out_prefix, ".lifted_records.tsv"))
write_tsv(tm[, .N, by = class_code][order(-N)], paste0(opt$out_prefix, ".tmap_class_codes.tsv"))
unpersonalized <- exon_pairs[!key2(locus, tgt_tx) %in% key2(ok$locus, ok$tgt_tx)]
in_loci <- uniqueN(v[allele_ok == TRUE]$locus)
summary <- rbind(
  data.table(metric = "source_records", value = nrow(v)),
  data.table(metric = "source_transcripts", value = uniqueN(v$CHROM)),
  data.table(metric = "source_genomic_loci_liftable_alleles", value = in_loci),
  v[, .(value = .N), by = .(metric = paste0("source_status:", status))][order(metric)],
  data.table(metric = c("lifted_records", "lifted_genomic_loci", "lifted_target_transcripts",
                        "genomic_loci_lost", "lifted_records_from_projection_only",
                        "exonic_locus_target_pairs_not_personalized",
                        "target_transcripts_with_unpersonalized_exonic_loci"),
             value = c(nrow(ok), uniqueN(ok$locus), uniqueN(ok$tgt_tx),
                       in_loci - uniqueN(ok$locus), sum(ok$origin == "exon_projection"),
                       nrow(unpersonalized), uniqueN(unpersonalized$tgt_tx))))
write_tsv(summary, paste0(opt$out_prefix, ".summary.tsv"))
close(fa)
print(summary)
message("Wrote ", out_vcf)
