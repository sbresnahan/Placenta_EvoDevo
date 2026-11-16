#!/usr/bin/env Rscript
f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), fixed = TRUE)
root <- normalizePath(file.path(dirname(f), ".."), mustWork = TRUE)
source(file.path(root, "R/library_paths.R"))
source(file.path(root, "R/h1_backend.R"))
source(file.path(root, "R/rstan_repair.R"))
a <- commandArgs(TRUE)
usage <- "Usage: bin/Rscript_sif bin/repair_rstan_stack.R check | fetch SOURCE_DIR | stage WORK [--sources SOURCE_DIR] | apply WORK --maintenance-confirmed | rollback WORK --maintenance-confirmed"
if (!length(a) || !a[1L] %in% c("check", "fetch", "stage", "apply", "rollback")) stop(usage)
action <- a[1L]
if ((action == "check" && length(a) != 1L) ||
    (action == "fetch" && length(a) != 2L) ||
    (action == "stage" && !(length(a) == 2L || (length(a) == 4L && a[3L] == "--sources"))) ||
    (action %in% c("apply", "rollback") &&
      (length(a) != 3L || a[3L] != "--maintenance-confirmed"))) stop(usage)
result <- poe_rstan_repair(root, action, if (length(a) > 1L) a[2L] else NULL,
  confirmed = length(a) == 3L && a[3L] == "--maintenance-confirmed",
  sources = if (action == "stage" && length(a) == 4L) a[4L] else NULL)
if (action %in% c("check", "fetch")) print(result)
