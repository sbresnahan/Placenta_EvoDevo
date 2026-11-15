#!/usr/bin/env Rscript
# Build a spliced transcript FASTA from a genome FASTA and a GTF, using the
# pipeline's own GTF parser/audit (R/annotation.R) and the same splicing rules
# as the uncurated branch of build_reference(). The output therefore passes
# validate_source_transcripts(): identical transcript IDs and spliced lengths.
#
# Usage (via the container wrapper):
#   bin/Rscript_sif bin/build_transcript_fasta.R GENOME_FA GTF OUT_FA [GENOME_FAI]
#
# GENOME_FAI defaults to GENOME_FA.fai; pass a separate index when the genome
# directory is read-only. Sequences are written in upper case. The output is
# written to OUT_FA.partial and renamed only on success; an existing OUT_FA is
# never overwritten.
f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
f <- gsub("~+~", " ", f, fixed = TRUE)
root <- normalizePath(file.path(dirname(f), ".."), mustWork = TRUE)
source(file.path(root, "R/library_paths.R"))
source(file.path(root, "R/common.R"))
source(file.path(root, "R/annotation.R"))
suppressPackageStartupMessages({
  library(Rsamtools)
  library(Biostrings)
})

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 3L || length(a) > 4L)
  stop("Usage: build_transcript_fasta.R GENOME_FA GTF OUT_FA [GENOME_FAI]")
genome <- a[1]; gtf <- a[2]; out <- a[3]
fai <- if (length(a) == 4L) a[4] else paste0(genome, ".fai")
for (p in c(genome, gtf, fai)) if (!file.exists(p)) stop("Missing input: ", p)
if (file.exists(out)) stop("Refusing to overwrite existing FASTA: ", out)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

cs <- fread(fai, header = FALSE, select = 1:2, col.names = c("chrom", "length"))
sizes <- tempfile(fileext = ".sizes")
fwrite(cs, sizes, sep = "\t", col.names = FALSE)

z <- read_gtf(gtf)
x <- audit_gtf(z, sizes)   # same exon/transcript/contig checks as the reference stage
message("Transcripts: ", nrow(x$transcripts), "; exons: ", nrow(x$exons),
        "; genes: ", uniqueN(x$transcripts$gene_id))

fa <- FaFile(genome, index = fai)
open(fa); on.exit(close(fa), add = TRUE)
ids <- sort(x$transcripts$transcript_id)
tmp <- paste0(out, ".partial"); unlink(tmp)
for (i in seq(1L, length(ids), by = 1000L)) {
  e <- x$exons[transcript_id %in% ids[i:min(i + 999L, length(ids))]]
  setorder(e, transcript_id, tx_exon_order)
  seqs <- scanFa(fa, GenomicRanges::GRanges(e$chrom, IRanges::IRanges(e$start, e$end)))
  minus <- which(e$strand == "-")
  if (length(minus)) seqs[minus] <- reverseComplement(seqs[minus])
  e[, sequence := toupper(as.character(seqs))]
  t <- e[, .(sequence = paste0(sequence, collapse = "")), by = transcript_id]
  s <- DNAStringSet(t$sequence); names(s) <- t$transcript_id
  if (any(width(s) != x$transcripts$spliced_length[match(names(s), x$transcripts$transcript_id)]))
    stop("Transcript length mismatch")
  writeXStringSet(s, tmp, append = TRUE)
}
if (!file.rename(tmp, out)) stop("Could not finalize ", out)
message("Wrote ", length(ids), " transcripts to ", out)
