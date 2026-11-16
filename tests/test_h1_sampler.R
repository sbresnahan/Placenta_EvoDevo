#!/usr/bin/env Rscript
# Compiled backend integration check using observed counts with software labels.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("Usage: test_h1_sampler.R SOURCE_ROOT PRIVATE_TEST_DIR OBSERVED_COUNT_FIXTURE.rds")
root <- normalizePath(args[1]); work <- args[2]; input <- readRDS(args[3])
suppressPackageStartupMessages(library(data.table))
for (m in c("common", "regulatory_common", "baseline_nb", "h1_bayes"))
  source(file.path(root, "R", paste0(m, ".R")))
dir.create(work, recursive = TRUE, showWarnings = FALSE)
cfg <- list(profile = "public_cross", dosage = list(maternal = 1, paternal = 1),
  h1_bayes = list(count_provenance = "assignment_covariance_unknown",
    producing_pipeline = "glmmTMB Salamanders observed example counts with software labels",
    warmup = 250L, iterations = 250L, parallel_chains = 1L))
hs <- h1_settings(cfg)
if (identical(Sys.getenv("H1_SMOKE_FULL"), "true")) {
  hs$warmup <- 1000L; hs$iterations <- 1000L
}
s <- list(min_snps = 3L, min_samples = 4L, min_total_counts = 50, min_per_cross = 2L)
counts <- h1_counts(input$data, input$annotation, cfg, s)
stopifnot(all(counts$coverage$eligibility == "eligible"))
encoded <- h1_encode(counts$data)
prior_table <- brms::get_prior(h1_formula(), data = encoded$data, family = brms::beta_binomial())
stopifnot(any(prior_table$class == "b" & prior_table$coef == "Intercept" & prior_table$dpar == ""))
nb_write(as.data.table(prior_table), file.path(work, "available_priors.tsv"))
code <- brms::make_stancode(h1_formula(), data = encoded$data,
  family = brms::beta_binomial(), prior = h1_priors(h1_profiles()[1]))
for (i in seq_len(nrow(h1_profiles()))) stopifnot(nchar(brms::make_stancode(h1_formula(),
  data = encoded$data, family = brms::beta_binomial(), prior = h1_priors(h1_profiles()[i]))) > 0L)
writeLines(code, file.path(work, "h1_model.stan"))
started <- proc.time()[["elapsed"]]
result <- h1_fit_once(counts$data, input$annotation, hs, h1_profiles()[1], 4271L,
  file.path(work, "sampler_fit.rds"))
elapsed <- proc.time()[["elapsed"]] - started
stopifnot(all(c("Delta_mag", "Delta_signed", "TE_magnitude_mean",
  "nonTE_magnitude_mean") %in% result$summary$variable))
stopifnot(all(is.finite(result$summary$mean)),
  nrow(result$residuals) == 2L, nrow(result$diagnostics) == 1L)
saved <- readRDS(file.path(work, "sampler_fit.rds"))
stopifnot(posterior::nchains(saved$derived) == 4L, posterior::niterations(saved$derived) == hs$iterations,
  nrow(saved$diagnostics$energy) == 4L)
nb_write(result$summary, file.path(work, "sampler_posterior_software_fixture.tsv"))
nb_write(result$diagnostics, file.path(work, "sampler_diagnostics.tsv"))
nb_write(result$residuals, file.path(work, "sampler_residuals_software_fixture.tsv"))
receipt <- data.table(check = c("compiled_joint_model", "four_sampled_chains",
  "posterior_contrasts", "derived_diagnostics", "per_chain_energy", "observed_residuals"),
  status = "pass", role = "software_backend_smoke_test", elapsed_seconds = elapsed,
  brms = as.character(packageVersion("brms")), cmdstanr = as.character(packageVersion("cmdstanr")),
  posterior = as.character(packageVersion("posterior")), cmdstan = as.character(cmdstanr::cmdstan_version()),
  warmup = hs$warmup, retained_per_chain = hs$iterations)
nb_write(receipt, file.path(work, "sampler_validation_receipt.tsv"))
print(receipt)
print(result$diagnostics)
