#!/usr/bin/env Rscript
# Standalone entry point; also routed from poe.R and run_pipeline.R.
f <- gsub("~+~", " ", sub("^--file=", "",
  grep("^--file=", commandArgs(), value = TRUE)[1]), fixed = TRUE)
root <- normalizePath(dirname(f), mustWork = TRUE)
source(file.path(root, "R/library_paths.R"))
source(file.path(root, "R/regulatory_cli.R"))
reg_main(commandArgs(trailingOnly = TRUE), root)
