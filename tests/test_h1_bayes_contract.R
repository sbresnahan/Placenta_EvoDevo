#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
suppressPackageStartupMessages(library(data.table))
for (m in c("common", "regulatory_common", "baseline_nb", "h1_bayes"))
  source(file.path(root, "R", paste0(m, ".R")))
checks <- 0L
check <- function(label, expr) {stopifnot(isTRUE(expr)); checks <<- checks + 1L; cat("PASS ", label, "\n", sep = "")}
fails <- function(expr) inherits(tryCatch({force(expr); NULL}, error = identity), "error")
cfg <- list(profile = "public_cross", dosage = list(maternal = 2, paternal = 1),
  h1_bayes = list(count_provenance = "assignment_covariance_unknown",
    producing_pipeline = "deterministic software fixture"))
s <- list(min_snps = 3L, min_samples = 4L, min_total_counts = 50, min_per_cross = 2L)
ann <- data.table(transcript_id = c("T1", "T2"), gene_id = "G", promoter_te_label = c("TE", "nonTE"))
# Fixed format/count fixture used only to exercise support gates.
z <- CJ(feature_id = ann$transcript_id, sample_id = paste0("S", 1:4), j = 1:3)
z[, `:=`(SNPpos = paste0(feature_id, ":marker", j), maternal_count = 4.4,
  paternal_count = 2.2, cluster_id = sample_id, independence_verified = TRUE,
  cross = ifelse(sample_id %in% c("S1", "S2"), "focal_mother", "alternate_mother"))]
r <- h1_counts(z, ann, cfg, s)
check("count support ignores upstream GLMM status", all(r$coverage$eligibility == "eligible"))
check("parental counts rounded separately", all(r$data$M == 4 & r$data$P == 2 & r$data$N == 6))
check("dosage is conditional logit offset", all(r$data$dosage_offset == log(.5)))
check("reciprocal effects centered equally", setequal(r$data$cross_c, c(-.5, .5)))
check("own markers retained", uniqueN(r$data$SNPpos) == 6)
zz <- copy(z); zz[1, maternal_count := NA_real_]; zz[2, `:=`(maternal_count = 0, paternal_count = 0)]
rr <- h1_counts(zz, ann, cfg, s)
check("missing parental assignments excluded", nrow(rr$data) == nrow(z) - 2L)
check("rounding and zero rows audited", sum(rr$coverage$n_missing_parental_rows) == 1 &&
  sum(rr$coverage$n_zero_total_rows) == 1 &&
  sum(rr$coverage$n_fractional_parental_values) > 0)
zz <- copy(z); zz[sample_id == "S2", cluster_id := "S1"]
check("independent clusters gate transcripts", all(h1_counts(zz, ann, cfg, s)$coverage$eligibility ==
  "insufficient_reciprocal_clusters"))
zz <- copy(z); zz[sample_id == "S1", independence_verified := FALSE]
check("unverified independence is explicit", any(h1_counts(zz, ann, cfg, s)$coverage$eligibility == "independence_unverified"))
check("five prespecified prior profiles", nrow(h1_profiles()) == 5 &&
  identical(h1_profiles()$precision_center, c(20, 20, 20, 5, 80)))
hs <- h1_settings(cfg)
check("sampling defaults", hs$chains == 4 && hs$warmup == 1000 && hs$iterations == 1000)
bad <- cfg; bad$h1_bayes$count_provenance <- "known_fragment_duplication"
check("known duplicate fragments blocked", fails(h1_settings(bad)))
bad$h1_bayes$count_provenance <- NULL
check("undeclared count provenance blocked", fails(h1_settings(bad)))
if (requireNamespace("posterior", quietly = TRUE)) {
  variables <- c("b_Intercept", "b_cross_c", "r_transcript[tx000001,Intercept]",
    "r_transcript[tx000001,cross_c]", "r_transcript[tx000002,Intercept]", "r_transcript[tx000002,cross_c]")
  dd <- array(rep(c(-1, 4, -.5, 0, .5, 0), each = 8), c(2, 4, 6),
    dimnames = list(NULL, paste0("chain", 1:4), variables))
  out <- h1_derived(dd, h1_encode(r$data)$map, ann)
  check("magnitude uses balanced beta", all(out[, , "Delta_mag"] == 1))
  check("signed direction is separate", all(out[, , "Delta_signed"] == -1))
  check("cross cancellation recorded", all(out[, , "cancellation:T1"] > 0))
}
cat("H1 contract checks: ", checks, "\n", sep = "")
