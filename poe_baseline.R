#!/usr/bin/env Rscript
# Versioned all-SNP model workflow, called through bin/Rscript_sif on deployment.
f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), fixed = TRUE)
root <- normalizePath(dirname(f), mustWork = TRUE)
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages({library(data.table); library(GenomicRanges)})
for (module in c("config", "common", "workflow", "parental", "phase_origin", "gusto",
    "models", "acquisition", "public", "reference", "pipeline", "stages", "additional",
    "annotation", "regulatory_common", "regulatory_resources", "regulatory_fish",
    "regulatory_te_genes", "regulatory_context", "baseline_nb", "baseline_te", "baseline_handoff"))
  source(file.path(root, "R", paste0(module, ".R")))
setDTthreads(1)
args <- commandArgs(trailingOnly = TRUE)
usage <- paste("Usage: poe_baseline.R CONFIG preflight|prepare|pilot|fit LEVEL BATCH|combine|",
  "te-inputs|te-pilot|te-fit JOB_INDEX|te-combine|between-genes|context|report-ready|status\n",
  "       poe_baseline.R cross-summary OUTPUT.tsv CONFIG [CONFIG ...]")
if (length(args) < 2L) stop(usage)
if (args[1] == "cross-summary") {
  if (length(args) < 3L) stop(usage)
  nb_cross_summary(args[-c(1, 2)], args[2], root)
  quit(status = 0)
}
cfg <- read_config_v2(args[1]); s <- nb_settings(cfg); action <- args[2]
switch(action,
  preflight = {p <- nb_preflight(cfg, s, root); cat("Validated ", nrow(p$metadata),
    " samples and ", nrow(p$files), " input files. No model has run.\n", sep = "")},
  prepare = nb_prepare(cfg, s, root),
  pilot = nb_pilot(cfg, s, root, if (length(args) >= 3L) as.integer(args[3]) else 10L),
  fit = {if (length(args) != 4L || !args[3] %in% c("gene", "transcript")) stop(usage)
    nb_fit_batch(cfg, s, root, args[3], as.integer(args[4]))},
  combine = nb_combine(cfg, s, root),
  `te-inputs` = nb_te_prepare(cfg, s, root),
  `te-pilot` = nb_pilot(cfg, s, root, if (length(args) >= 3L) as.integer(args[3]) else 5L, TRUE),
  `te-fit` = {if (length(args) != 3L) stop(usage); nb_te_gene(cfg, s, root, as.integer(args[3]))},
  `te-combine` = nb_te_combine(cfg, s, root),
  `between-genes` = nb_between_genes(cfg, s, root),
  context = nb_context(cfg, s, root),
  `report-ready` = nb_report_ready(cfg, s, root),
  status = print(nb_status(cfg, s, root)),
  stop(usage))
