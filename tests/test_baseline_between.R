#!/usr/bin/env Rscript
# Between-gene software fixture. Observed example counts are transformed into
# arbitrary test outcomes and paired covariates. These are NOT POE estimates.
a <- commandArgs(TRUE)
root <- if (length(a)) normalizePath(a[1]) else normalizePath(".")
suppressPackageStartupMessages(library(data.table))
for (module in c("config", "common", "regulatory_common", "regulatory_te_genes"))
  source(file.path(root, "R", paste0(module, ".R")))
checks <- 0L
check <- function(name, expr) {stopifnot(isTRUE(expr)); checks <<- checks + 1L; cat("PASS ", name, "\n", sep = "")}
data("Salamanders", package = "glmmTMB", envir = environment())
home <- tempfile("nb-between-"); dir.create(home)
n <- 24L
d <- data.table(gene_id = sprintf("software_gene%02d", seq_len(2*n)),
  pair = rep(seq_len(n), each = 2), te = rep(c(TRUE, FALSE), n),
  chrom = "chr1", not_strictly_POE = FALSE, independence_verified = TRUE)
d[, `:=`(start = 100 + .I*10000, end = 1000 + .I*10000,
  estimate = log1p(Salamanders$count[seq_len(2*n)]),
  se = 1 / sqrt(1 + Salamanders$count[101:(100+2*n)]),
  block = paste0("block", ceiling(pair/3)))]
for (j in seq_along(TE_GENES_COVARIATES))
  d[, (TE_GENES_COVARIATES[j]) := log1p(Salamanders$count[pair + j*10])]
d[, magnitude_unbiased := estimate^2 - se^2]
tss <- d[, .(gene_id, chrom, tss = start)]
repeats <- d[te == TRUE, .(chrom, start, end = start + 1, is_te = TRUE, is_unknown = FALSE)]
inp <- list(genes = d[, !c("pair", "te")], tss = tss, repeats = repeats,
  te_u = reg_union(repeats), contigs = data.table(chrom = "chr1", length = 1e7),
  scale = "SOFTWARE_FIXTURE_NOT_POE", source_table = "software_fixture", dropped_incomplete = 0L)
z <- list(output_dir = home, seed = 1147L, bootstrap = 20L, version = "software_fixture",
  min_depth = 10, te_genes = list(permutations = 20, placebo_shifts = 5))
s <- reg_te_genes_settings(z)
out <- file.path(home, TE_GENES_FILES)
reg_te_genes_run(list(cohort_id = "software_fixture", species = "mouse"), z, inp, s, out)
r <- reg_read(out[3]); table <- reg_read(out[1])
check("between-gene nonempty adjusted estimator runs",
  r$status %in% c("exploratory_propensity_subclassification", "exploratory_residual_imbalance_gt_0.1"))
check("between-gene outcome provenance is new fixture only", r$gene_POE_source == "software_fixture")
check("between-gene cross and independence flags retained",
  all(c("not_strictly_POE", "independence_verified") %in% names(table)))
check("between-gene uncertainty and placebo fields exported",
  all(c("magnitude_ci_low", "signed_ci_high", "signed_p_placebo",
    "bootstrap_valid_fraction") %in% names(r)))
check("between-gene balance table complete", nrow(reg_read(out[2])) == length(TE_GENES_COVARIATES))
cat("Completed ", checks, " between-gene software checks. Not biological POE results.\n", sep = "")
unlink(home, recursive = TRUE)
