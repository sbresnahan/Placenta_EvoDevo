#!/usr/bin/env Rscript
# Write independent full cohort configurations; execution is operator-controlled.
f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), fixed = TRUE)
root <- normalizePath(file.path(dirname(f), ".."), mustWork = TRUE)
source(file.path(root, "R/library_paths.R"))
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("Usage: prepare_promoter_configs.R PROJECT OPERATOR_SETTINGS.yaml DEST_CONFIG_DIR")
suppressPackageStartupMessages(library(yaml))
source(file.path(root, "R/h1_annotations.R"))
project <- normalizePath(args[1], mustWork = TRUE)
mapping <- function(x, label) {
  if (!is.list(x) || !length(x) || is.null(names(x)) ||
      anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)))
    stop(label, " must be a nonempty YAML mapping")
  invisible(x)
}
read_mapping <- function(path, label, operator = FALSE) {
  if (!file.exists(path) || dir.exists(path) || file.access(path, 4L) != 0L) {
    guidance <- if (operator) paste0(
      "\nCreate the private operator file from config/promoter_operator_settings.example.yaml, ",
      "then review its required fields; see runbook section 1.") else ""
    stop(label, " is missing or unreadable: ", path, guidance)
  }
  x <- tryCatch(yaml::read_yaml(path), error = function(e)
    stop("Cannot parse ", label, " at ", path, ": ", conditionMessage(e)))
  mapping(x, label)
  x
}
settings <- read_mapping(args[2], "Operator settings YAML", operator = TRUE)
allowed <- c("h1_bayes", "gusto_human_context", "exposure_by_cohort")
unknown <- setdiff(names(settings), allowed)
if (length(unknown)) stop("Unknown operator settings field: ", paste(unknown, collapse = ", "))
settings$h1_bayes <- h1_normalize_annotations(settings$h1_bayes, "operator.h1_bayes")
reviewed <- function(x, field) {
  if (is.list(x)) {
    labels <- names(x)
    if (is.null(labels)) labels <- as.character(seq_along(x))
    for (i in seq_along(x)) reviewed(x[[i]], paste(field, labels[i], sep = "."))
  } else if (is.character(x) && any(grepl(
      "REVIEW_REQUIRED|DECLARE_|TODO|UNSPECIFIED|/actual/path|/absolute/path", x,
      ignore.case = TRUE), na.rm = TRUE))
    stop(field, " contains an unresolved template placeholder; edit and review it first")
}
reviewed(settings, "operator")
h <- settings$h1_bayes
h1_validate_controls(h, "operator.h1_bayes")
if (!is.null(settings$gusto_human_context)) {
  mapping(settings$gusto_human_context, "operator.gusto_human_context")
  path_keys <- intersect(names(settings$gusto_human_context),
    c("phenotype_map", "gxe_tables", "exposure_manifest", "qtl_tables", "beds"))
  for (key in path_keys) {
    value <- settings$gusto_human_context[[key]]
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !grepl("^/", value) || !nzchar(trimws(value)))
      stop("operator.gusto_human_context.", key, " must be one absolute manifest/map path")
  }
}
cohorts <- c("GUSTO", "mouse_Babr_CAST", "mouse_FVB_CAST_E12",
  "maize_B73_Ki11", "maize_B73_W22", "fish_Pprolifica_AB")
if (!is.null(settings$exposure_by_cohort)) {
  mapping(settings$exposure_by_cohort, "operator.exposure_by_cohort")
  if (length(setdiff(names(settings$exposure_by_cohort), cohorts)))
    stop("Unknown cohort in operator.exposure_by_cohort")
  for (cohort in names(settings$exposure_by_cohort))
    mapping(settings$exposure_by_cohort[[cohort]], paste0("operator.exposure_by_cohort.", cohort))
}
destination <- args[3]
if (!grepl("^/", destination)) stop("Use an absolute destination directory")
if (file.exists(destination) && !dir.exists(destination)) stop("Destination is not a directory")
if (dir.exists(destination) && length(list.files(destination, all.files = TRUE, no.. = TRUE)))
  stop("Destination config directory must be empty")
# Read all six sources before creating any destination or writing configs.
sources <- setNames(file.path(project, "data", cohorts, "private",
  paste0(cohorts, ".baseline_nb.yaml")), cohorts)
configs <- lapply(cohorts, function(cohort) {
  cfg <- read_mapping(sources[[cohort]], paste0("Upstream configuration for ", cohort))
  if (!identical(cfg$cohort_id, cohort)) stop("Source cohort ID disagrees: ", sources[[cohort]])
  for (k in c("paths", "templates")) {
    values <- unlist(cfg[[k]])
    if (anyNA(values) || any(!grepl("^/", values)))
      stop("Input config requires absolute ", k, ": ", cohort)
  }
  if (cohort != "GUSTO") {
    combined <- h1_merge_overrides(cfg$h1_bayes, h, paste0(cohort, ".h1_bayes"))
    reviewed(combined, paste0(cohort, ".h1_bayes"))
    cfg$h1_bayes <- if (length(combined)) combined else NULL
  }
  cfg
})
names(configs) <- cohorts
if (!dir.exists(destination) && !dir.create(destination, recursive = TRUE, showWarnings = FALSE))
  stop("Cannot create destination config directory: ", destination)
rows <- list()
for (cohort in cohorts) {
  cfg <- configs[[cohort]]
  for (radius in c(3000L, 0L, 1000L)) {
    d <- cfg
    run <- paste0("promoter_r", radius, "_run1")
    d$baseline_nb$output_dir <- file.path(project, "data", cohort, "results", run)
    d$regulatory$te_tss_window_bp <- radius
    if (cohort == "GUSTO") {
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
      destination_config = file, poe_input_config = sources[[cohort]])
  }
}
write.table(do.call(rbind, rows), file.path(destination, "CONFIGS.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE)
cat("Wrote ", length(rows), " configurations. Review CONFIGS.tsv and downstream input checks before analysis.\n", sep = "")
cat("Absent H1 count/pipeline annotations remain unrecorded; existing counts are unchanged.\n")
