#!/usr/bin/env Rscript
a <- commandArgs(TRUE)
root <- normalizePath(if (length(a)) a[1] else ".")
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(data.table))
for (m in c("common", "regulatory_common", "baseline_nb", "h1_bayes"))
  source(file.path(root, "R", paste0(m, ".R")))
main <- function() {
  checks <- character()
  check <- function(label, value) {
    if (!isTRUE(value)) stop("FAIL ", label)
    checks <<- c(checks, label); cat("PASS ", label, "\n", sep = "")
  }
  fails <- function(x) inherits(tryCatch({force(x); NULL}, error = identity), "error")
  raw <- function(p) readBin(p, "raw", n = file.info(p)$size)
  home <- tempfile("optional H1 annotations "); dir.create(home)
  on.exit(unlink(home, recursive = TRUE), add = TRUE)
  hs <- h1_settings(list())
  numeric_keys <- c("chains", "warmup", "iterations", "adapt_delta", "max_treedepth", "seed", "parallel_chains")
  expected <- list(chains = 4L, warmup = 1000L, iterations = 1000L,
    adapt_delta = .99, max_treedepth = 12L, seed = 1147L, parallel_chains = 4L,
    count_provenance = NA_character_, producing_pipeline = NA_character_)
  check("absent H1 blocks initialize with default sampling settings and missing annotations",
    identical(hs, expected))
  variants <- list(
    list(count_provenance = NULL, producing_pipeline = NULL),
    list(count_provenance = "", producing_pipeline = " "),
    list(count_provenance = NA_character_, producing_pipeline = NA_character_),
    list(count_provenance = "REVIEW_REQUIRED_COUNT_PROVENANCE",
      producing_pipeline = "REVIEW_REQUIRED_COUNT_PRODUCER"),
    list(producing_pipeline = "DECLARE_REVIEWED_COUNT_PRODUCER"))
  for (i in seq_along(variants)) check(paste("optional missing annotation form", i),
    identical(h1_settings(list(h1_bayes = variants[[i]])), expected))
  valid <- list(chains = 6L, parallel_chains = 2L, warmup = 1200L, iterations = 1500L,
    adapt_delta = .97, max_treedepth = 13L, seed = 2718L,
    count_provenance = "upstream_duplicate_corrected", producing_pipeline = "already recorded count producer")
  # Build the expected settings list for valid annotations to verify task identity.
  prior_settings <- expected
  for (key in names(valid)) prior_settings[[key]] <- valid[[key]]
  check("valid existing settings retain exact saved-task identity",
    identical(h1_settings(list(h1_bayes = valid)), prior_settings))
  placeholder <- variants[[4]]
  merged <- h1_merge_overrides(valid, c(placeholder, list(seed = 1147L)))
  check("annotation placeholders never erase inherited valid annotations",
    identical(merged$count_provenance, valid$count_provenance) &&
      identical(merged$producing_pipeline, valid$producing_pipeline) &&
      merged$chains == 6L && merged$parallel_chains == 2L && merged$seed == 1147L)
  for (key in c("count_provenance", "producing_pipeline")) {
    bad <- list(); bad[[key]] <- 123
    check(paste("invalid non-text annotation rejected", key), fails(h1_settings(list(h1_bayes = bad))))
  }
  for (value in c("guessed", "known_fragment_duplication")) check(paste("explicit invalid assignment rejected", value),
    fails(h1_settings(list(h1_bayes = list(count_provenance = value)))))
  check("unrecognized producer placeholder is not accepted as real metadata",
    fails(h1_settings(list(h1_bayes = list(producing_pipeline = "TODO choose producer")))))
  for (item in list(list(chains = 3L), list(parallel_chains = 5L), list(warmup = 0L),
      list(iterations = "1000"), list(adapt_delta = .9), list(adapt_delta = 1))) {
    check(paste("sampling constraint remains enforced", names(item)), fails(h1_settings(list(h1_bayes = item))))
  }
  project <- file.path(home, "project"); dir.create(project)
  cohorts <- c("GUSTO", "mouse_Babr_CAST", "mouse_FVB_CAST_E12",
    "maize_B73_Ki11", "maize_B73_W22", "fish_Pprolifica_AB")
  sources <- setNames(file.path(project, "data", cohorts, "private",
    paste0(cohorts, ".baseline_nb.yaml")), cohorts)
  context <- list(phenotype_map = file.path(project, "existing/phenotype_map.tsv"),
    qtl_tables = file.path(project, "existing/qtl_inputs.tsv"),
    beds = file.path(project, "existing/beds.tsv"), require_qtl = TRUE)
  for (id in cohorts) {
    cfg <- yaml::read_yaml(file.path(root, "config", paste0(id, ".example.yaml")))
    cfg$h1_bayes <- NULL
    cfg$regulatory$human_context <- context
    cfg$exposure_nb <- list(fixture_marker = id)
    dir.create(dirname(sources[[id]]), recursive = TRUE); yaml::write_yaml(cfg, sources[[id]])
  }
  private <- file.path(project, "data/private")
  dir.create(private, recursive = TRUE)
  count_file <- file.path(private, "unchanged_count_batch.fixture")
  # Sentinel is not scientific count data; setup/generation must not open it.
  writeLines("IMMUTABLE_SOFTWARE_SENTINEL_NOT_A_COUNT_MATRIX", count_file)
  manifest_dir <- file.path(project, "data/GUSTO/private/qtl_context_v1")
  dir.create(manifest_dir)
  gxe <- file.path(manifest_dir, "gxe_tables.ppBMI_v1.tsv")
  exposure <- file.path(manifest_dir, "ppBMI_exposure.ppBMI_v1.tsv")
  receipt <- file.path(manifest_dir, "ppBMI_setup.ppBMI_v1.yaml")
  for (p in c(gxe, exposure, receipt)) writeLines("EXISTING_SETUP_OUTPUT_SENTINEL", p)
  operator <- file.path(private, "promoter_operator_settings.ppBMI_v1.yaml")
  op <- yaml::read_yaml(file.path(root, "config/promoter_operator_settings.example.yaml"))
  op$h1_bayes$count_provenance <- placeholder$count_provenance
  op$h1_bayes$producing_pipeline <- placeholder$producing_pipeline
  op$gusto_human_context <- c(context["phenotype_map"],
    list(gxe_tables = gxe, exposure_manifest = exposure))
  op$exposure_by_cohort <- list(GUSTO = list(fixture_marker = "operator exposure"))
  yaml::write_yaml(op, operator)
  protected <- c(unname(sources), operator, count_file, gxe, exposure, receipt)
  before <- lapply(protected, raw)
  invoke <- function(op_path, destination) {
    out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      shQuote(c(file.path(root, "bin/prepare_promoter_configs.R"), project, op_path, destination)),
      stdout = TRUE, stderr = TRUE))
    list(code = if (is.null(attr(out, "status"))) 0L else attr(out, "status"),
      text = paste(out, collapse = "\n"))
  }
  dest <- file.path(private, "promoter_configs_ppBMI_v1")
  result <- invoke(operator, dest)
  if (result$code != 0L) stop(result$text)
  inventory <- read.delim(file.path(dest, "CONFIGS.tsv"))
  configs <- lapply(inventory$destination_config, yaml::read_yaml)
  reciprocal <- configs[inventory$cohort_id != "GUSTO"]
  initialized <- lapply(reciprocal, h1_settings)
  check("placeholder operator generates all eighteen configs directly",
    nrow(inventory) == 18L && setequal(inventory$cohort_id, cohorts) &&
      setequal(inventory$radius, c(0L, 1000L, 3000L)))
  check("all fifteen reciprocal configs initialize H1 with missing metadata and no compiler",
    length(initialized) == 15L &&
      all(vapply(initialized, function(x) is.na(x$count_provenance) &&
        is.na(x$producing_pipeline) && identical(x[numeric_keys], expected[numeric_keys]), TRUE)))
  check("generated reciprocal YAML contains no invented assignment or pipeline labels",
    all(vapply(reciprocal, function(x) is.null(x$h1_bayes$count_provenance) &&
      is.null(x$h1_bayes$producing_pipeline), TRUE)))
  check("generated GUSTO context inherits original QTL BED and completed ppBMI paths",
    all(vapply(configs[inventory$cohort_id == "GUSTO"], function(x) {
      z <- x$regulatory$human_context
      identical(z$qtl_tables, context$qtl_tables) && identical(z$beds, context$beds) &&
        identical(z$gxe_tables, gxe) && identical(z$exposure_manifest, exposure)
    }, TRUE)))
  check("existing per-cohort and optional exposure settings remain inherited",
    all(vapply(seq_along(configs), function(i) identical(configs[[i]]$exposure_nb$fixture_marker,
      if (inventory$cohort_id[i] == "GUSTO") "operator exposure" else inventory$cohort_id[i]), TRUE)))
  check("no source count upstream YAML operator manifest or setup receipt is modified",
    identical(before, lapply(protected, raw)))
  check("generator does not run setup or require valid sentinel aggregate contents",
    result$code == 0L && grepl("annotations remain unrecorded", result$text, fixed = TRUE))
  # Missing H1 overrides must keep each cohort's existing controls and labels.
  cohort <- cohorts[2]
  cfg <- yaml::read_yaml(sources[[cohort]]); cfg$h1_bayes <- valid
  yaml::write_yaml(cfg, sources[[cohort]])
  no_override <- op; no_override$h1_bayes <- NULL
  op2 <- file.path(home, "no_h1_override.yaml"); yaml::write_yaml(no_override, op2)
  dest2 <- file.path(home, "inherited configs")
  result <- invoke(op2, dest2)
  inherited <- yaml::read_yaml(file.path(dest2, paste0(cohort, ".r3000.yaml")))
  check("omitted operator H1 preserves cohort-specific controls and annotations",
    result$code == 0L && identical(inherited$h1_bayes, valid))
  control_only <- op
  control_only$h1_bayes <- c(placeholder, list(seed = 1147L))
  op3 <- file.path(home, "placeholder_overlay.yaml"); yaml::write_yaml(control_only, op3)
  dest3 <- file.path(home, "merged configs")
  result <- invoke(op3, dest3)
  combined <- yaml::read_yaml(file.path(dest3, paste0(cohort, ".r3000.yaml")))$h1_bayes
  check("generator merges controls without erasing inherited annotations",
    result$code == 0L && identical(combined$count_provenance, valid$count_provenance) &&
      identical(combined$producing_pipeline, valid$producing_pipeline) &&
      combined$chains == valid$chains && combined$seed == 1147L)
  bad <- op; bad$h1_bayes$unrelated_field <- "REVIEW_REQUIRED_OTHER_SETTING"
  bad_file <- file.path(home, "unrelated_placeholder.yaml"); yaml::write_yaml(bad, bad_file)
  bad_dest <- file.path(home, "blocked configs")
  result <- invoke(bad_file, bad_dest)
  check("unrelated required placeholders still block before creating configs",
    result$code != 0L && grepl("operator.h1_bayes.unrelated_field", result$text, fixed = TRUE) &&
      !dir.exists(bad_dest))
  result <- invoke(operator, dest)
  check("nonempty config destinations remain protected",
    result$code != 0L && grepl("Destination config directory must be empty", result$text, fixed = TRUE))
  if (length(a) >= 2L) {
    dir.create(a[2], recursive = TRUE, showWarnings = FALSE)
    write.table(data.frame(check = checks, status = "pass", evidence = "local software fixture"),
      file.path(a[2], "h1_optional_annotations_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  cat("Optional H1 annotation checks: ", length(checks),
    "; unchanged count sentinels, no fitting or compiler access\n", sep = "")
}
main()
