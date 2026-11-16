#!/usr/bin/env Rscript
# Metadata-only setup; no aggregate scans, fits or external job execution.
f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), fixed = TRUE)
root <- normalizePath(file.path(dirname(f), ".."), mustWork = TRUE)
source(file.path(root, "R/library_paths.R"))
args <- commandArgs(TRUE)
if (!length(args) || length(args) > 3L)
  stop("Usage: setup_external_ppbmi.R PROJECT [SOURCE_PATHS.yaml [EXISTING_OPERATOR.yaml]]")
suppressPackageStartupMessages(library(yaml))
source(file.path(root, "R/h1_annotations.R"))
source(file.path(root, "R/external_ppbmi_setup.R"))
setup_external_ppbmi(root, args[1],
  if (length(args) >= 2L) args[2] else file.path(root, "config/external_ppbmi.yaml"),
  if (length(args) >= 3L) args[3] else NULL)
