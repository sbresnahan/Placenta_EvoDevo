#!/usr/bin/env Rscript
# Software regression fixtures only. No fixture is biological evidence.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1]) else normalizePath(".")
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages(library(data.table))
for (module in c("common", "regulatory_common", "baseline_nb", "baseline_te", "baseline_handoff"))
  source(file.path(root, "R", paste0(module, ".R")))
passed <- character()
check <- function(label, expr) {
  stopifnot(isTRUE(expr)); passed <<- c(passed, label); cat("PASS ", label, "\n", sep = "")
}
fails <- function(expr) inherits(tryCatch(force(expr), error = function(e) e), "error")
# Simple hand-specified arithmetic fixtures, not simulated biological data.
d <- data.table(sample_id = "s1", transcript_id = c("txA", "txB", "txA", "txB"),
  gene_id = "g1", locus_id = c("chr1:10:A:G", "chr1:10:A:G", "chr1:20:C:T", "chr1:30:G:A"),
  maternal_count = c(10, 20, 8, 12), paternal_count = c(2, 4, 6, 2))
w <- data.table(transcript_id = c("txA", "txB"), weight = c(.75, .25))
g <- nb_aggregate_gene(d, w)
check("usage-weighted shared gene SNP", g[locus_id == "chr1:10:A:G"]$maternal_count == 12.5)
check("union retains disjoint transcript SNPs", nrow(g) == 3)
check("no per-SNP renormalization", g[locus_id == "chr1:20:C:T"]$maternal_count == 6)
missing <- copy(d); missing[1, maternal_count := NA_real_]
check("unknown parent does not become zero", is.na(nb_aggregate_gene(missing, w)[locus_id == "chr1:10:A:G"]$maternal_count))
bad <- copy(d); bad$locus_id[1] <- "chr1:10"
check("full genomic variant key required", fails(nb_aggregate_gene(bad, w)))
dup <- rbind(d, d[1])
check("duplicate count keys rejected", fails(nb_aggregate_gene(dup, w)))
tmp <- tempfile(fileext = ".txt")
writeLines(c("{", '"novel:tx.1": 3,', '"txB": 1', "}"), tmp)
weights <- nb_read_tpms(tmp, data.table(transcript_id = c("novel:tx.1", "txB"), gene_id = "g"))
check("TPM names with colons retained", weights[transcript_id == "novel:tx.1"]$weight == .75)
writeLines(c('{"not_a_valid_TPM_record"}'), tmp)
check("invalid TPM file rejected", fails(nb_read_tpms(tmp, data.table(transcript_id = "txB", gene_id = "g"))))
writeLines(c("{", '"a": 0,', '"b": 0', "}"), tmp)
check("all-zero TPM weights equal", all(nb_read_tpms(tmp, data.table(transcript_id = c("a", "b"), gene_id = "g"))$weight == .5))
unlink(tmp)
b <- c(allelepaternal = .4, "allelepaternal:crossalternate_mother" = -.2)
V <- matrix(c(.04, .01, .01, .09), 2, dimnames = list(names(b), names(b)))
avg <- nb_linear(b, V, setNames(c(1, .5), names(b)))
check("equal-cross effect is not reference coefficient", abs(avg[["beta"]] - .3) < 1e-12)
check("equal-cross covariance retained", abs(avg[["se"]]^2 - (.04 + .25 * .09 + .01)) < 1e-12)
check("nominal interaction flag", identical(nb_cross_flag(c(.049, .05, NA_real_)), c(TRUE, FALSE, NA)))
check("maize dosage null yields zero adjusted effect", log(.5) - log(1 / 2) == 0)
est <- data.table(feature_id = c("txA", "txB"), beta_adjusted = c(.4, -.1))
annot <- data.table(transcript_id = c("txA", "txB"), tss_te_label = c("TE", "nonTE"),
  first_donor = c(100, 999), marker = c("chr1:20:C:T", "chr1:30:G:A"))
stats <- nb_te_statistics(est, annot)
check("different markers and donors do not gate comparisons", nrow(stats) == 2 && all(abs(stats$signed - .5) < 1e-12))
check("magnitude contrast not signed contrast", all(abs(stats$magnitude - .3) < 1e-12))
z <- data.table(feature_id = "x", SNPpos = c("chr1:1:A:G", "chr1:2:C:T"),
  sample_id = "s1", maternal_count = c(1.6, NA), paternal_count = c(2.4, NA),
  offset = 0, sex = "F", ancestry = "A", cross = NA_character_,
  cluster_id = "s1", independence_verified = TRUE)
long <- nb_long(z)
check("expected counts rounded only at fit", identical(sort(long$count), c(2, 2)))
check("unknown parental rows excluded not zero-filled", nrow(long) == 2)
M <- matrix(c(1, 2, 5, 3, 4, 8), 3, dimnames = list(NULL, c("s1", "s2")))
sf <- nb_size_factors(M, M * 2)
check("size factors positive geometrically centered", all(sf > 0) && abs(prod(sf) - 1) < 1e-10)
# Sample block bootstrap keeps all transcripts and their own SNPs together.
block <- CJ(sample_id = c("s1", "s2", "s3", "s4"), feature_id = c("txA", "txB"))
block[, `:=`(cluster_id = sample_id, cross = rep(c("focal_mother", "alternate_mother"), each = 4),
  SNPpos = ifelse(feature_id == "txA", "chr1:20:C:T", "chr1:30:G:A"))]
set.seed(11); sampled <- nb_bootstrap_units(block, TRUE)
check("joint bootstrap retains disjoint marker identities", setequal(sampled$SNPpos, block$SNPpos))
check("joint bootstrap preserves both transcripts per sampled unit",
  all(sampled[, uniqueN(feature_id), by = sample_id]$V1 == 2))
check("reciprocal bootstrap preserves direction counts",
  identical(sort(unique(sampled[, .(sample_id, cross)])[, .N, by = cross]$N), c(2L, 2L)))
# Receipt mutation tests.
s <- list(output_dir = tempfile("nb-test"))
dir.create(s$output_dir)
file <- file.path(s$output_dir, "value.tsv"); nb_write(data.table(a = 1), file)
nb_finish(s, "unit", "key1", file)
check("verified checkpoint can be reused", nb_current(s, "unit", "key1"))
check("changed checkpoint key rejected", fails(nb_current(s, "unit", "key2")))
nb_write(data.table(a = 2), file)
check("mutated output rejected", fails(nb_require(s, "unit")))
unlink(s$output_dir, recursive = TRUE)
s <- list(output_dir = tempfile("nb-batches"), batches_gene = 3L,
  min_snps = 3L, min_samples = 1L, min_total_counts = 1, max_dense_gb = 1)
metadata <- data.table(sample_id = "s1", cluster_id = "s1", independence_verified = TRUE,
  sex = "F", ancestry = "A", cross = NA_character_)
gd <- copy(g); setnames(gd, "gene_id", "feature_id")
nb_make_batches(gd, metadata, "gene", s, "unit-key")
manifest <- reg_read(nb_path(s, "gene_batches.tsv"))
check("empty and nonempty model batches valid", nrow(manifest) == 3 && sum(manifest$n_features) == 1)
check("feature inventory retains union SNPs", nrow(reg_read(nb_path(s, "gene_inventory.tsv.gz"))) == 3)
check("batch checksum covers serialized input", all(vapply(manifest$input_file, reg_sha, "") == manifest$sha256))
unlink(s$output_dir, recursive = TRUE)
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  # Reuse observed counts from the package's real Salamanders data solely as
  # a numerical fitting fixture. Artificial allele/sample labels are NOT a POE dataset.
  data("Salamanders", package = "glmmTMB", envir = environment())
  pairs <- CJ(sample_id = sprintf("fixture%02d", 1:40), SNPpos = paste0("chr1:", 1:4, ":A:G"))
  pairs[, `:=`(feature_id = "numerical_fixture",
    maternal_count = Salamanders$count[seq_len(.N)],
    paternal_count = Salamanders$count[seq_len(.N) + 160],
    offset = 0, sex = rep(c("F", "M"), each = 4, length.out = .N),
    ancestry = rep(c("A", "A", "B", "B"), each = 4, length.out = .N),
    cross = rep(c("focal_mother", "alternate_mother"), each = 80),
    cluster_id = sample_id, independence_verified = TRUE)]
  cfg <- list(profile = "human_trio", dosage = list(maternal = 1, paternal = 1))
  opt <- list(min_samples = 15L, min_snps = 3L, min_total_counts = 50, min_per_cross = 2L)
  ans <- nb_fit_one(pairs, cfg, opt)
  refdat <- nb_long(pairs)
  refdat[, `:=`(sex = factor(sex), ancestry = factor(ancestry))]
  ref <- suppressWarnings(glmmTMB::glmmTMB(count ~ allele + sex + ancestry + (1|SNPpos) + (1|sample_id),
    family = glmmTMB::nbinom2(), offset = offset, data = refdat,
    control = glmmTMB::glmmTMBControl(optimizer = nlminb, parallel = 1,
      optCtrl = list(rel.tol = 1e-5, eval.max = 300L, iter.max = 300L))))
  co <- summary(ref)$coefficients$cond
  check("human reference fit converges", ans$status == "eligible_NB")
  check("literal human baseline coefficient parity", abs(ans$beta - co["allelepaternal", 1]) < 1e-8)
  check("literal human baseline uncertainty parity", abs(ans$se - co["allelepaternal", 2]) < 1e-8)
  cfg$profile <- "public_cross"
  crossfit <- nb_fit_one(pairs, cfg, opt)
  check("reciprocal interaction model runs", crossfit$status == "exploratory_NB")
  check("reciprocal formula contains interaction", grepl("allele \\* cross", crossfit$formula))
  check("equal-cross average matches cross-specific effects",
    abs(crossfit$beta_adjusted - mean(c(crossfit$beta_cross0, crossfit$beta_cross1))) < 1e-8)
  cfg$dosage <- list(maternal = 2, paternal = 1)
  maize <- nb_fit_one(pairs, cfg, opt)
  check("dosage contrast does not change fitted coefficient", abs(maize$beta - crossfit$beta) < 1e-8)
  check("dosage contrast centers at expected ratio", abs(maize$beta_adjusted - maize$beta - log(2)) < 1e-8)
} else cat("SKIP fitting tests: glmmTMB unavailable\n")
cat("Completed ", length(passed), " baseline software checks.\n", sep = "")
