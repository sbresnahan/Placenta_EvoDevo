#!/usr/bin/env Rscript
# Raw adapter regression with observed example counts and artificial labels.
# Tests actual REF/ALT parsing, coordinate mapping and parental orientation.
a <- commandArgs(TRUE)
root <- if (length(a)) normalizePath(a[1]) else normalizePath(".")
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages({library(data.table); library(GenomicRanges)})
for (module in c("common", "parental", "phase_origin", "gusto", "public", "pipeline"))
  source(file.path(root, "R", paste0(module, ".R")))
checks <- 0L
check <- function(name, expr) {stopifnot(isTRUE(expr)); checks <<- checks + 1L; cat("PASS ", name, "\n", sep = "")}
data("Salamanders", package = "glmmTMB", envir = environment())
home <- tempfile("nb-adapter-"); dir.create(home)
tx <- data.table(transcript_id = c("fixture_A", "fixture_B", "fixture_overlap"),
  gene_id = c("fixture_gene", "fixture_gene", "other_gene"))
ex <- data.table(transcript_id = tx$transcript_id, gene_id = tx$gene_id,
  chrom = "chr1", strand = c("+", "-", "+"), start = c(10, 200, 10),
  end = c(49, 239, 10), tx_exon_order = 1L)
em <- data.table(transcript_id = "fixture_A", gene_id = "fixture_gene",
  event_id = "fixture_event", window = 50L)
ev <- data.table(event_id = "fixture_event", chrom = "chr1",
  distinguishing_intervals = "10-49", window = 50L)
markers <- data.table(chrom = "chr1", genomic_pos = c(10, 11, 239),
  ref = c("A", "A", "T"), alt = c("G", "G", "C"), focal_index = 0L, alternate_index = 1L)
counts <- data.table(transcript_id = c(rep("fixture_A", 3), "fixture_B"),
  tx_pos = c(1L, 2L, 3L, 1L), tx_ref = "A", tx_alt = "G", unused = 0,
  ref_count = Salamanders$count[11:14], alt_count = Salamanders$count[111:114], other_count = 0)
p <- setNames(file.path(home, paste0(c("tx_map", "exons", "event_map", "events", "markers"), ".tsv")),
  c("tx_map", "exons", "event_map", "events", "markers"))
for (name in names(p)) fwrite(switch(name, tx_map = tx, exons = ex, event_map = em,
  events = ev, markers = markers), p[[name]], sep = "\t")
raw <- file.path(home, "counts.txt"); fwrite(counts, raw, sep = "\t", col.names = FALSE)
cfg <- list(profile = "public_cross", ornaments = list(semantics = "ref_alt"),
  paths = as.list(p), parents = list(focal = "focal", alternate = "alternate"))
m <- list(libraries = data.table(sample_id = "fixture_sample", unit_id = "fixture_unit",
  strand_mode = "none"), units = data.table(unit_id = "fixture_unit", maternal_parent = "focal"))
cfg_path <- function(c, k, sample, ...) if (k == "allele_counts") raw else
  if (k == "child_tx_map") file.path(home, "child_tx.tsv") else stop("Unexpected path: ", k)
write_tsv <- function(...) stop("Legacy write reached by baseline adapter")
d <- assign_sample(cfg, m, "fixture_sample", return_transcripts = TRUE)
check("actual raw adapter retains transcript outside event map", "fixture_B" %in% d$transcript_id)
check("actual adapter retains disjoint SNPs", uniqueN(d$locus_id) == 4)
check("public REF count is maternal with focal mother",
  d[transcript_id == "fixture_A" & tx_pos == 2]$maternal_count == counts$ref_count[2])
check("negative strand SNP maps with complemented genomic alleles",
  d[transcript_id == "fixture_B"]$locus_id == "chr1:239:T:C")
check("other-gene overlap is missing, not zero",
  is.na(d[transcript_id == "fixture_A" & tx_pos == 1]$maternal_count))
check("missing public marker remains missing",
  is.na(d[transcript_id == "fixture_A" & tx_pos == 3]$paternal_count))
m$units[, maternal_parent := "alternate"]
r <- assign_sample(cfg, m, "fixture_sample", return_transcripts = TRUE)
check("reciprocal public direction swaps parent counts",
  r[transcript_id == "fixture_A" & tx_pos == 2]$maternal_count == counts$alt_count[2])

cfg$profile <- "human_trio"
cfg$gusto_precomputed <- list(enabled = TRUE, origin_matrix = file.path(home, "origin.tsv"))
human_linkage <- function(c, m) data.table(sample_id = "fixture_sample",
  child_id = "fixture_child", family_id = "fixture_family")
anchors <- data.table(transcript_id = counts$transcript_id, tx_pos = counts$tx_pos,
  tx_ref = counts$tx_ref, tx_alt = counts$tx_alt, gt = "0|1",
  locus_id = c("chr1:10:A:G", "chr1:11:A:G", "chr1:12:A:G", "chr1:239:T:C"),
  chrom = "chr1", genomic_pos = c(10, 11, 12, 239))
fwrite(anchors, file.path(home, "child_tx.tsv"), sep = "\t")
origin <- data.table(chr = "chr1", fixture_child_L = "maternal", fixture_child_R = "paternal")
fwrite(origin, cfg$gusto_precomputed$origin_matrix, sep = "\t")
g <- assign_sample(cfg, m, "fixture_sample", return_transcripts = TRUE)
check("GUSTO phased parent origin uses REF/ALT schema",
  g[transcript_id == "fixture_B"]$maternal_count == counts$ref_count[4])
origin[, `:=`(fixture_child_L = "paternal", fixture_child_R = "maternal")]
fwrite(origin, cfg$gusto_precomputed$origin_matrix, sep = "\t")
g <- assign_sample(cfg, m, "fixture_sample", return_transcripts = TRUE)
check("GUSTO reversed phase swaps maternal count",
  g[transcript_id == "fixture_B"]$maternal_count == counts$alt_count[4])
origin[, `:=`(fixture_child_L = "undefined", fixture_child_R = "undefined")]
fwrite(origin, cfg$gusto_precomputed$origin_matrix, sep = "\t")
g <- assign_sample(cfg, m, "fixture_sample", return_transcripts = TRUE)
check("undefined GUSTO origin retains paired missingness",
  all(is.na(g$maternal_count)) && all(is.na(g$paternal_count)))
cat("Completed ", checks, " raw-adapter software checks. No private data used.\n", sep = "")
unlink(home, recursive = TRUE)
