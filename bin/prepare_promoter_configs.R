#!/usr/bin/env Rscript
# Write independent full cohort configurations; execution is operator-controlled.
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("Usage: prepare_promoter_configs.R PROJECT OPERATOR_SETTINGS.yaml DEST_CONFIG_DIR")
suppressPackageStartupMessages(library(yaml))
project <- normalizePath(args[1], mustWork = TRUE)
settings <- yaml::read_yaml(args[2])
destination <- args[3]
if (!grepl("^/", destination)) stop("Use an absolute destination directory")
if (dir.exists(destination) && length(list.files(destination, all.files = TRUE, no.. = TRUE)))
  stop("Destination config directory must be empty")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
cohorts <- c("GUSTO", "mouse_Babr_CAST", "mouse_FVB_CAST_E12", "maize_B73_Ki11",
  "maize_B73_W22", "fish_Pprolifica_AB")
rows <- list()
for (cohort in cohorts) {
  source <- file.path(project, "data", cohort, "private", paste0(cohort, ".baseline_nb.yaml"))
  cfg <- yaml::read_yaml(source)
  if (!identical(cfg$cohort_id, cohort)) stop("Source cohort ID disagrees: ", source)
  for (k in c("paths", "templates")) {
    values <- unlist(cfg[[k]])
    if (any(!grepl("^/", values))) stop("Input config requires absolute ", k, ": ", cohort)
  }
  for (radius in c(3000L, 0L, 1000L)) {
    d <- cfg
    run <- paste0("promoter_r", radius, "_run1")
    d$baseline_nb$output_dir <- file.path(project, "data", cohort, "results", run)
    d$regulatory$te_tss_window_bp <- radius
    if (cohort != "GUSTO") {
      if (is.null(settings$h1_bayes)) stop("Supply reviewed h1_bayes settings")
      d$h1_bayes <- settings$h1_bayes
    } else {
      for (key in names(settings$gusto_human_context))
        d$regulatory$human_context[[key]] <- settings$gusto_human_context[[key]]
    }
    if (!is.null(settings$exposure_by_cohort[[cohort]]))
      d$exposure_nb <- settings$exposure_by_cohort[[cohort]]
    if (!is.null(d$exposure_nb)) {
      d$exposure_nb$baseline_dir <- d$baseline_nb$output_dir
      d$exposure_nb$output_dir <- file.path(project, "data", cohort, "results",
        paste0("exposure_r", radius, "_run1"))
    }
    file <- file.path(destination, paste0(cohort, ".r", radius, ".yaml"))
    yaml::write_yaml(d, file)
    rows[[length(rows) + 1L]] <- data.frame(cohort_id = cohort, radius = radius,
      destination_config = file, poe_input_config = source)
  }
}
write.table(do.call(rbind, rows), file.path(destination, "CONFIGS.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE)
