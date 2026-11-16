#!/usr/bin/env Rscript
# Deterministic error-path orchestration; every model invocation deliberately fails.
args <- commandArgs(TRUE); root <- normalizePath(if (length(args)) args[1] else ".")
suppressPackageStartupMessages(library(data.table))
for (m in c("common", "config", "regulatory_common", "baseline_nb", "baseline_te", "h1_bayes"))
  source(file.path(root, "R", paste0(m, ".R")))
checks <- 0L
check <- function(label, expr) {stopifnot(isTRUE(expr)); checks <<- checks + 1L; cat("PASS ", label, "\n", sep = "")}
if (!all(vapply(c("brms", "cmdstanr", "posterior"), requireNamespace, TRUE, quietly = TRUE)))
  stop("Bayesian backend packages required for H1 orchestration integration test")
home <- tempfile("h1-failure-contract-"); dir.create(home)
cfg <- list(profile = "public_cross", cohort_id = "software", species = "mouse",
  dosage = list(maternal = 1, paternal = 1), environment_allowlist = character(),
  paths = list(model_gene = file.path(home, "unused.tsv"), scratch = file.path(home, "scratch")),
  reference = list(build = "software"),
  h1_bayes = list(count_provenance = "assignment_covariance_unknown",
    producing_pipeline = "deterministic software error fixture"))
attr(cfg, "file") <- file.path(home, "cfg.yaml")
s <- list(output_dir = home, min_snps = 3L, min_samples = 4L, min_per_cross = 2L, min_total_counts = 50)
ann <- data.table(transcript_id = c("T1", "T2", "absent"), gene_id = c("G", "G", "H"),
  promoter_te_label = c("TE", "nonTE", "nonTE"), status = "upstream_GLMM_failure")
nb_write(ann, nb_path(s, "transcript_TE_annotation.tsv"))
z <- CJ(feature_id = c("T1", "T2"), sample_id = paste0("S", 1:4), marker = 1:3)
z[, `:=`(SNPpos = paste0(feature_id, ":", marker), maternal_count = 4, paternal_count = 2,
  cluster_id = sample_id, independence_verified = TRUE,
  cross = ifelse(sample_id %in% c("S1", "S2"), "focal_mother", "alternate_mother"))]
input <- file.path(home, "counts.rds"); saveRDS(list(data = z), input)
nb_write(data.table(feature_id = c("T1", "T2"), batch = 1L), nb_path(s, "transcript_feature_batches.tsv"))
nb_write(data.table(batch = 1L, input_file = input), nb_path(s, "transcript_batches.tsv"))
nb_finish(s, "prepare", list(software = TRUE), input)
nb_finish(s, "te-inputs", list(software = TRUE), nb_path(s, "transcript_TE_annotation.tsv"))
h1_prepare(cfg, s, root)
tc <- reg_read(nb_path(s, "h1_transcript_coverage.tsv"))
check("upstream model failures remain count-eligible", all(tc[feature_id %in% c("T1", "T2")]$eligibility == "eligible"))
check("unobserved annotation remains in coverage", tc[feature_id == "absent"]$eligibility == "no_informative_counts")
check("one gene has mixed-class support", nrow(reg_read(nb_path(s, "h1_bayes_jobs.tsv"))) == 1)
h1_fit_once <- function(...) stop("deliberate software error-path invocation")
h1_gene(cfg, s, root, 1L)
di <- reg_read(nb_path(s, "h1_bayes", "gene_000001_diagnostics.tsv"))
check("all five prior profile failures retained", sum(di$status == "fit_failure") == 5L)
check("LOCO never relaxes reciprocal support", sum(di$status == "unavailable_LOCO_support") == 4L)
check("omission labels omit private cluster IDs", all(!grepl("^S", di$omitted_cluster)))
h1_combine(cfg, s, root)
cov <- reg_read(nb_path(s, "h1_coverage.tsv"))
check("failed eligible gene survives combine", cov[gene_id == "G"]$fit_status == "fit_failure")
check("ineligible gene survives combine", cov[gene_id == "H"]$fit_status == "ineligible")
check("no posterior states its actual availability", reg_read(nb_path(s, "h1_cohort_descriptive.tsv"))$status ==
  "no_posterior_available")
unlink(home, recursive = TRUE)
cat("H1 workflow checks: ", checks, "\n", sep = "")
