#!/usr/bin/env Rscript
# Local full-config and container-wrapper mocks; no cluster access.
args <- commandArgs(TRUE); root <- normalizePath(if (length(args)) args[1] else ".")
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages({library(data.table); library(yaml)})
source(file.path(root, "R/config.R"))
checks <- 0L
check <- function(label, expr) {stopifnot(isTRUE(expr)); checks <<- checks + 1L; cat("PASS ", label, "\n", sep = "")}
home <- tempfile("promoter-deployment-"); dir.create(home)
project <- file.path(home, "project")
cohorts <- c("GUSTO", "mouse_Babr_CAST", "mouse_FVB_CAST_E12", "maize_B73_Ki11",
  "maize_B73_W22", "fish_Pprolifica_AB")
for (id in cohorts) {
  cfg <- yaml::read_yaml(file.path(root, "config", paste0(id, ".example.yaml")))
  cfg$baseline_nb <- list(output_dir = file.path(project, "data", id, "results", "baseline_nb_v1_run1"))
  path <- file.path(project, "data", id, "private", paste0(id, ".baseline_nb.yaml"))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); yaml::write_yaml(cfg, path)
}
op <- file.path(home, "operator.yaml")
yaml::write_yaml(list(h1_bayes = list(count_provenance = "assignment_covariance_unknown",
  producing_pipeline = "software declaration")), op)
dest <- file.path(home, "configs")
out <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(root, "bin/prepare_promoter_configs.R")),
  shQuote(project), shQuote(op), shQuote(dest)), stdout = TRUE, stderr = TRUE)
if (!is.null(attr(out, "status"))) stop(paste(out, collapse = "\n"))
manifest <- fread(file.path(dest, "CONFIGS.tsv"))
check("all six cohorts and three radii configured", nrow(manifest) == 18 && uniqueN(manifest$cohort_id) == 6 &&
  setequal(manifest$radius, c(3000, 0, 1000)))
parsed <- lapply(manifest$destination_config, read_config_v2)
check("all generated configurations satisfy schema", length(parsed) == 18L)
check("independent output directories", uniqueN(vapply(parsed, function(x) x$baseline_nb$output_dir, "")) == 18L)
check("upstream POE pointers remain separate", all(vapply(seq_along(parsed), function(i)
  parsed[[i]]$baseline_nb$output_dir != read_config_v2(manifest$poe_input_config[i])$baseline_nb$output_dir, TRUE)))
lib <- file.path(home, "Rlib"); dir.create(lib)
sif <- file.path(home, "R.sif"); file.create(sif)
capture <- file.path(home, "container_args")
mock <- file.path(home, "apptainer")
writeLines(c("#!/bin/bash", paste0("printf '%s\\n' \"$@\" > ", shQuote(capture))), mock)
Sys.chmod(mock, "0755")
site <- file.path(home, "site.sh")
writeLines(c(paste0("POE_R_SIF=", shQuote(sif)), paste0("POE_R_LIBRARY=", shQuote(lib)),
  paste0("POE_SINGULARITY=", shQuote(mock))), site)
env <- c(paste0("SCRIPTS_DIR=", shQuote(root)), paste0("POE_SITE_FILE=", shQuote(site)),
  "CMDSTAN=/absolute/compiled_cmdstan")
rc <- system2("bash", c(shQuote(file.path(root, "bin/Rscript_sif")), "-e", shQuote("1+1")),
  env = env)
check("container wrapper forwards CmdStan through clean environment",
  rc == 0L && "CMDSTAN=/absolute/compiled_cmdstan" %in% readLines(capture))
texts <- list.files(root, recursive = TRUE, full.names = TRUE)
texts <- texts[grepl("[.](md|R|py|yaml|sh)$", texts)]
pattern <- paste(c(paste0("migra", "tion"), paste0("super", "sed"),
  paste0("lega", "cy"), paste0("rather", " than")), collapse = "|")
hits <- texts[vapply(texts, function(p) any(grepl(pattern, readLines(p, warn = FALSE), ignore.case = TRUE)), TRUE)]
check("maintained prose passes direct-writing guard", !length(hits))
unlink(home, recursive = TRUE)
cat("Promoter deployment checks: ", checks, "\n", sep = "")
