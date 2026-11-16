#!/usr/bin/env Rscript
args_all <- commandArgs()
script <- sub("^--file=", "", args_all[grepl("^--file=", args_all)][1])
root <- normalizePath(dirname(script), mustWork = TRUE)
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages(library(data.table))
for (m in c("config", "exposure_inputs", "te_expression", "exposure_nb",
            "exposure_mediation", "exposure_followup", "exposure_workflow", "exposure_plots"))
  source(file.path(root, "R", paste0(m, ".R")))
data.table::setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE)
usage <- paste("Usage: poe_exposure.R CONFIG",
  "preflight|te-plan|te-import|prepare|pilot [N]|fit BATCH|mediation BATCH|combine|bootstrap INDEX|",
  "enrichment|report-ready|export|status")
if (length(args) < 2L) stop(usage)
cfg <- read_config_v2(args[1]); s <- ex_settings(cfg); action <- args[2]
index <- function() {
  if (length(args) != 3L || is.na(suppressWarnings(as.integer(args[3]))) || as.integer(args[3]) < 1L)
    stop("Supply a positive integer index")
  as.integer(args[3])
}
switch(action,
  preflight = print(ex_preflight(cfg, s)$audit),
  `te-plan` = ex_te_plan(cfg, s),
  `te-import` = {
    loci <- ex_read(s$te_loci); ex_validate_loci(loci, cfg)
    counts <- ex_import_telescope(ex_read(s$libraries), loci)
    ex_write(counts, s$te_counts)
    cat("Imported counts. Supply and validate the completed quantification QC before prepare.\n")
  },
  prepare = print(ex_prepare(cfg, s)),
  pilot = print(ex_pilot(cfg, s, if (length(args) == 3L) index() else 9L)),
  fit = ex_fit_batch(cfg, s, index(), mediation = FALSE),
  mediation = ex_fit_batch(cfg, s, index(), mediation = TRUE),
  combine = ex_combine(cfg, s),
  bootstrap = ex_bootstrap_job(cfg, s, index()),
  enrichment = ex_run_enrichment(cfg, s),
  `report-ready` = {
    d <- ex_read(ex_path(s, "model_results.tsv"))
    ex_write(d[, .(planned = .N, estimable = sum(status == "eligible"),
      discoveries = sum(discovery)), by = .(module, level, predictor)], ex_path(s, "report_summary.tsv"))
    print(ex_status(s))
  },
  export = print(ex_export(s)),
  status = print(ex_status(s)),
  stop(usage))
