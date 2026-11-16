#!/usr/bin/env Rscript
a <- commandArgs(TRUE)
root <- normalizePath(if (length(a)) a[1] else ".")
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages(library(yaml))
main <- function() {
  checks <- character()
  check <- function(label, value) {
    if (!isTRUE(value)) stop("FAIL ", label)
    checks <<- c(checks, label); cat("PASS ", label, "\n", sep = "")
  }
  home <- tempfile("operator setup "); dir.create(home)
  on.exit(unlink(home, recursive = TRUE), add = TRUE)
  project <- file.path(home, "project with spaces"); dir.create(project)
  cohorts <- c("GUSTO", "mouse_Babr_CAST", "mouse_FVB_CAST_E12",
    "maize_B73_Ki11", "maize_B73_W22", "fish_Pprolifica_AB")
  raw_file <- function(path) readBin(path, "raw", n = file.info(path)$size)
  upstream <- setNames(file.path(project, "data", cohorts, "private",
    paste0(cohorts, ".baseline_nb.yaml")), cohorts)
  context <- list(phenotype_map = file.path(project, "existing/phenotype_map.tsv"),
    qtl_tables = file.path(project, "existing/qtl_inputs.tsv"),
    beds = file.path(project, "existing/beds.tsv"), require_qtl = TRUE)
  for (id in cohorts) {
    cfg <- yaml::read_yaml(file.path(root, "config", paste0(id, ".example.yaml")))
    cfg$baseline_nb <- list(output_dir = file.path(project, "data", id, "results", "baseline_nb_v1_run1"))
    cfg$regulatory$human_context <- context
    cfg$exposure_nb <- list(baseline_dir = cfg$baseline_nb$output_dir,
      output_dir = file.path(project, "old_exposure", id), fixture_marker = id)
    dir.create(dirname(upstream[[id]]), recursive = TRUE)
    yaml::write_yaml(cfg, upstream[[id]])
  }
  original <- lapply(upstream, raw_file)
  template_path <- file.path(root, "config/promoter_operator_settings.example.yaml")
  template <- yaml::read_yaml(template_path)
  expected_map <- paste0("/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/",
    "data/GUSTO/private/qtl_context_v1/phenotype_map.tsv")
  check("template uses the phenotype-map location recorded in the supplied tree",
    identical(template$gusto_human_context$phenotype_map, expected_map))
  # Execute the exact marked creation block from the delivered runbook.
  doc <- readLines(file.path(root, "docs/Seadragon_analysis_runbook.md"), warn = FALSE)
  marker <- which(doc == "<!-- operator-create -->")
  stopifnot(length(marker) == 1L)
  begin <- which(seq_along(doc) > marker & doc == "```bash")[1L]
  end <- which(seq_along(doc) > begin & doc == "```")[1L]
  shell <- file.path(home, "documented_create.sh")
  writeLines(c('set -e', 'cd "$1"', 'PROJECT="$2"',
    doc[seq.int(begin + 1L, end - 1L)]), shell)
  create <- function() {
    out <- suppressWarnings(system2("bash", shQuote(c(shell, root, project)),
      stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(out, "status"))) stop(paste(out, collapse = "\n"))
    out
  }
  create()
  operator <- file.path(project, "data/private/promoter_operator_settings.yaml")
  check("documented commands create the missing private file from the template",
    file.exists(operator) && identical(raw_file(operator), raw_file(template_path)))
  check("new operator file is private",
    bitwAnd(as.integer(file.info(operator)$mode), 511L) == 384L)
  writeLines(c(readLines(operator), "# USER_EDIT_PRESERVED"), operator)
  edited <- raw_file(operator); create_output <- create()
  check("rerunning documented creation preserves an existing operator file",
    identical(raw_file(operator), edited) &&
      any(grepl("Keeping existing operator file", create_output, fixed = TRUE)))
  invoke <- function(op, dest) {
    out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      shQuote(c(file.path(root, "bin/prepare_promoter_configs.R"), project, op, dest)),
      stdout = TRUE, stderr = TRUE))
    list(code = if (is.null(attr(out, "status"))) 0L else attr(out, "status"),
      text = paste(out, collapse = "\n"))
  }
  rejected <- function(label, op, expected, dest = file.path(home, paste0("bad-output-", length(checks)))) {
    result <- invoke(op, dest)
    ok <- result$code != 0L && grepl(expected, result$text, fixed = TRUE) && !dir.exists(dest)
    if (!ok) message(result$text)
    check(label, ok)
  }
  rejected("missing operator file gives creation guidance before writing",
    file.path(home, "missing.yaml"), "config/promoter_operator_settings.example.yaml")
  rejected("a directory cannot be used as the operator file",
    home, "Operator settings YAML is missing or unreadable")
  rejected("unreviewed template stops before creating the config destination",
    operator, "operator.gusto_human_context.gxe_tables contains an unresolved template placeholder")
  complete <- template
  complete$h1_bayes$count_provenance <- "assignment_covariance_unknown"
  complete$h1_bayes$producing_pipeline <- "explicit software fixture declaration, not cohort evidence"
  complete$gusto_human_context$phenotype_map <- context$phenotype_map
  complete$gusto_human_context$gxe_tables <- file.path(project, "later/gxe_tables.tsv")
  complete$gusto_human_context$exposure_manifest <- file.path(project, "later/ppBMI_exposure.tsv")
  pending <- template
  pending$h1_bayes <- complete$h1_bayes
  pending_path <- file.path(home, "pending.yaml"); yaml::write_yaml(pending, pending_path)
  rejected("GxE placeholder is reported by its actual operator field",
    pending_path, "operator.gusto_human_context.gxe_tables contains an unresolved template placeholder")
  cases <- list(
    list(label = "invalid H1 mapping", x = list(h1_bayes = "not a mapping", gusto_human_context = context),
      expected = "operator.h1_bayes must be a YAML mapping"),
    list(label = "invalid count provenance", field = "count_provenance", value = "guessed",
      expected = "operator.h1_bayes.count_provenance must be one of"),
    list(label = "non-text count producer", field = "producing_pipeline", value = 123,
      expected = "operator.h1_bayes.producing_pipeline must be one optional text annotation"),
    list(label = "count-producer placeholder recognized by H1 runtime",
      field = "producing_pipeline", value = "DECLARE actual pipeline",
      expected = "operator.h1_bayes.producing_pipeline must identify"),
    list(label = "too few chains", field = "chains", value = 3L,
      expected = "requires chains >= 4"),
    list(label = "zero warmup", field = "warmup", value = 0L,
      expected = "Invalid operator.h1_bayes.warmup"),
    list(label = "text iterations", field = "iterations", value = "1000",
      expected = "Invalid operator.h1_bayes.iterations"),
    list(label = "excess parallel chains", field = "parallel_chains", value = 5L,
      expected = "parallel_chains <= chains"),
    list(label = "invalid adapt delta", field = "adapt_delta", value = .90,
      expected = "operator.h1_bayes.adapt_delta"),
    list(label = "relative GxE manifest", section = "gusto_human_context", field = "gxe_tables",
      value = "relative.tsv", expected = "operator.gusto_human_context.gxe_tables must be one absolute"),
    list(label = "unknown top-level key", x = c(complete, list(misspelled_context = TRUE)),
      expected = "Unknown operator settings field"),
    list(label = "unknown exposure cohort", x = c(complete, list(exposure_by_cohort =
      list(UNKNOWN = list(enabled = TRUE)))), expected = "Unknown cohort in operator.exposure_by_cohort"),
    list(label = "nested exposure placeholder", x = c(complete, list(exposure_by_cohort =
      list(GUSTO = list(counts = "REVIEW_REQUIRED_COUNTS")))),
      expected = "operator.exposure_by_cohort.GUSTO.counts contains an unresolved template placeholder")
  )
  for (i in seq_along(cases)) {
    item <- cases[[i]]
    x <- if (!is.null(item$x)) item$x else complete
    if (!is.null(item$field)) {
      section <- if (is.null(item$section)) "h1_bayes" else item$section
      x[[section]][[item$field]] <- item$value
    }
    path <- file.path(home, paste0("invalid_", i, ".yaml"))
    yaml::write_yaml(x, path)
    rejected(paste("invalid operator rejected:", item$label), path, item$expected)
  }
  malformed <- file.path(home, "malformed.yaml"); writeLines("h1_bayes: [unterminated", malformed)
  rejected("malformed YAML has a parse-specific diagnostic", malformed, "Cannot parse Operator settings YAML")
  empty <- file.path(home, "empty.yaml"); writeLines("", empty)
  rejected("empty YAML cannot masquerade as reviewed settings", empty,
    "Operator settings YAML must be a nonempty YAML mapping")
  done <- file.path(home, "reviewed.yaml"); yaml::write_yaml(complete, done)
  last <- upstream[[tail(cohorts, 1L)]]
  stopifnot(file.rename(last, paste0(last, ".held")))
  result <- invoke(done, file.path(home, "missing-upstream-output"))
  stopifnot(file.rename(paste0(last, ".held"), last))
  check("all six upstream files are checked before any config is written",
    result$code != 0L && grepl("Upstream configuration for fish_Pprolifica_AB", result$text, fixed = TRUE) &&
      !dir.exists(file.path(home, "missing-upstream-output")))
  dest <- file.path(home, "complete configs")
  result <- invoke(done, dest)
  if (result$code != 0L) stop(result$text)
  inventory <- read.delim(file.path(dest, "CONFIGS.tsv"))
  check("reviewed software fixture generates six cohorts times three radii",
    nrow(inventory) == 18L && setequal(inventory$cohort_id, cohorts) &&
      setequal(inventory$radius, c(0L, 1000L, 3000L)))
  outputs <- lapply(inventory$destination_config, yaml::read_yaml)
  gusto <- outputs[inventory$cohort_id == "GUSTO"]
  check("GUSTO QTL and BED configuration is preserved while adding overrides",
    all(vapply(gusto, function(x) identical(x$regulatory$human_context$qtl_tables, context$qtl_tables) &&
      identical(x$regulatory$human_context$beds, context$beds) &&
      identical(x$regulatory$human_context$gxe_tables, complete$gusto_human_context$gxe_tables), TRUE)))
  check("later manifest contents are not falsely required or certified during generation",
    !file.exists(complete$gusto_human_context$gxe_tables) &&
      !file.exists(complete$gusto_human_context$exposure_manifest))
  check("existing per-cohort H6/H7 settings inherit with fresh radius-specific directories",
    all(vapply(seq_along(outputs), function(i) {
      x <- outputs[[i]]
      identical(x$exposure_nb$fixture_marker, inventory$cohort_id[i]) &&
        identical(x$exposure_nb$baseline_dir, x$baseline_nb$output_dir) &&
        endsWith(x$exposure_nb$output_dir, paste0("exposure_r", inventory$radius[i], "_run1"))
    }, TRUE)))
  minimal <- file.path(home, "inherit.yaml")
  yaml::write_yaml(list(h1_bayes = complete$h1_bayes), minimal)
  inherit_dest <- file.path(home, "inherited configs")
  result <- invoke(minimal, inherit_dest)
  check("omitted GUSTO override inherits upstream context without deleting fields",
    result$code == 0L && identical(yaml::read_yaml(file.path(inherit_dest,
      "GUSTO.r3000.yaml"))$regulatory$human_context, context))
  before_output <- lapply(inventory$destination_config, raw_file)
  result <- invoke(done, dest)
  check("nonempty output destination is rejected without overwriting configs",
    result$code != 0L && grepl("Destination config directory must be empty", result$text, fixed = TRUE) &&
      identical(before_output, lapply(inventory$destination_config, raw_file)))
  check("source configurations remain byte-identical",
    identical(original, lapply(upstream, raw_file)))
  if (length(a) >= 2L) {
    dir.create(a[2], recursive = TRUE, showWarnings = FALSE)
    write.table(data.frame(check = checks, status = "pass", evidence = "local software fixture"),
      file.path(a[2], "operator_setup_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  cat("Operator setup checks: ", length(checks), "; no real cohort inputs or model fitting\n", sep = "")
}
main()
