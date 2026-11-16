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
  home <- tempfile("ppBMI setup software "); dir.create(home)
  on.exit(unlink(home, recursive = TRUE), add = TRUE)
  raw <- function(p) readBin(p, "raw", n = file.info(p)$size)
  tsv <- function(x, p) write.table(x, p, sep = "\t", row.names = FALSE, quote = FALSE)
  fixture <- function(label, existing = TRUE) {
    base <- file.path(home, label); dir.create(base)
    project <- file.path(base, "project with spaces"); dir.create(project)
    gx <- file.path(base, "external gxe")
    dir.create(file.path(gx, "aggregate"), recursive = TRUE)
    dir.create(file.path(gx, "inputs"))
    producer <- file.path(base, "producer.tsv")
    p <- data.frame(exposure_id = c("GA", "ppBMI"), column = c("GA", "ppBMI"),
      transform = "none", min = c(20, 14), max = c(45, 70), enabled = c(0, 1))
    tsv(p, producer)
    runtime <- file.path(gx, "inputs/exposure_manifest.tsv")
    e <- data.frame(exposure_id = "ppBMI", column = "ppBMI", binary = "False",
      n = 417, n_missing = 0, mean = "24.939616306954434", sd = "5.520665327307657")
    tsv(e, runtime)
    # Observed summary counts supplied by the operator, used only as setup records.
    modalities <- c("RNA_editing", "alt_TSS", "alt_polyA", "expression",
      "intron_retention", "isoform_expression", "isoforms", "splicing", "stability")
    summary <- data.frame(ancestry = rep(c("EAS", "EUR"), each = 9),
      modality = rep(modalities, 2), exposure = "ppBMI",
      n_tested = c(3120, 29057, 25552, 12568, 3812, 36441, 20673, 56149, 4493,
        17472, 35773, 31816, 13128, 10341, 39514, 21884, 75424, 7365),
      n_sig = c(3, 46, 61, 14, 0, 61, 4, 109, 2, 0, 57, 33, 3, 0, 1, 0, 1673, 0))
    # Explicit shard-label software record; no shard effect values are invented.
    summary <- rbind(summary, transform(summary[1, ], exposure = "ppBMI.chr1"))
    summary_path <- file.path(gx, "aggregate/gxe_summary.tsv"); tsv(summary, summary_path)
    headers <- list(
      gxe_tier1_results.tsv.gz = c("ancestry", "modality", "exposure", "phenotype_id",
        "variant_id", "b_gi", "b_gi_se", "b_g", "b_e", "pval_nominal", "pval_beta", "qval"),
      gxe_tier1_meta.tsv.gz = c("modality", "exposure", "phenotype_id", "beta_meta",
        "se_meta", "meta_variant_id", "scope", "lead_concordant", "p_feature_acat",
        "q_feature_acat", "n_ancestries"),
      gxe_tier2_results.tsv.gz = c("modality", "exposure", "phenotype_id",
        "beta_meta", "se_meta", "variant_id", "p_meta", "q_meta"))
    write_gz <- function(name, columns) {
      con <- gzfile(file.path(gx, "aggregate", name), "wt")
      writeLines(c(paste(columns, collapse = "\t"), "UNPARSED_SOFTWARE_SENTINEL"), con); close(con)
    }
    for (name in names(headers)) write_gz(name, headers[[name]])
    settings <- yaml::read_yaml(file.path(root, "config/external_ppbmi.yaml"))
    settings$gxe_root <- gx; settings$producing_config <- producer
    source <- file.path(base, "sources.yaml"); yaml::write_yaml(settings, source)
    context <- list(phenotype_map = file.path(base, "phenotype_map.tsv"),
      qtl_tables = file.path(base, "qtl_inputs.tsv"), beds = file.path(base, "beds.tsv"),
      require_qtl = TRUE)
    for (resource in unlist(context[c("phenotype_map", "qtl_tables", "beds")]))
      writeLines("SOFTWARE_FIXTURE_NOT_A_VALIDATED_RESOURCE", resource)
    cfg <- yaml::read_yaml(file.path(root, "config/GUSTO.example.yaml"))
    cfg$regulatory$human_context <- context
    upstream <- file.path(project, "data/GUSTO/private/GUSTO.baseline_nb.yaml")
    dir.create(dirname(upstream), recursive = TRUE); yaml::write_yaml(cfg, upstream)
    op <- yaml::read_yaml(file.path(root, "config/promoter_operator_settings.example.yaml"))
    op$h1_bayes$count_provenance <- "assignment_covariance_unknown"
    op$h1_bayes$producing_pipeline <- "software fixture declaration; not cohort provenance"
    op$gusto_human_context$qtl_tables <- file.path(base, "operator_qtl.tsv")
    op$gusto_human_context$beds <- file.path(base, "operator_beds.tsv")
    op$exposure_by_cohort <- list(GUSTO = list(enabled = TRUE, fixture_marker = "retain me"))
    operator <- file.path(project, "data/private/promoter_operator_settings.yaml")
    if (existing) {
      dir.create(dirname(operator), recursive = TRUE); yaml::write_yaml(op, operator)
    }
    list(base = base, project = project, gx = gx, producer = producer, p = p,
      runtime = runtime, e = e, summary_path = summary_path, summary = summary,
      headers = headers, write_gz = write_gz, source = source, settings = settings,
      upstream = upstream, context = context, operator = operator, op = op,
      output = file.path(project, "data/private/promoter_operator_settings.ppBMI_v1.yaml"))
  }
  invoke <- function(x, extra = character()) {
    out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      shQuote(c(file.path(root, "bin/setup_external_ppbmi.R"), x$project, x$source, extra)),
      stdout = TRUE, stderr = TRUE))
    list(code = if (is.null(attr(out, "status"))) 0L else attr(out, "status"),
      text = paste(out, collapse = "\n"))
  }
  x <- fixture("complete")
  originals <- lapply(c(x$producer, x$runtime, x$summary_path, x$upstream, x$operator), raw)
  result <- invoke(x)
  if (result$code != 0L) stop(result$text)
  op <- yaml::read_yaml(x$output)
  manifest <- read.delim(op$gusto_human_context$gxe_tables, check.names = FALSE)
  exposure <- read.delim(op$gusto_human_context$exposure_manifest, check.names = FALSE)
  receipt_path <- file.path(dirname(op$gusto_human_context$gxe_tables), "ppBMI_setup.ppBMI_v1.yaml")
  receipt <- yaml::read_yaml(receipt_path)
  outputs <- c(op$gusto_human_context$gxe_tables, op$gusto_human_context$exposure_manifest,
    x$output, receipt_path)
  check("actual supplied ppBMI metadata derives an enabled continuous local row",
    exposure$exposure_type == "continuous" && exposure$transformation == "none" &&
      isTRUE(exposure$enabled))
  check("pooled population-SD scaling and inclusive source range are explicit",
    grepl("pooled across retained ancestries", exposure$scale_definition, fixed = TRUE) &&
      grepl("ddof=0", exposure$scale_definition, fixed = TRUE) &&
      grepl("14 <= ppBMI <= 70", exposure$range_definition, fixed = TRUE))
  check("producer and runtime provenance paths are real fixture inputs",
    identical(exposure$producing_config, x$producer) &&
      identical(exposure$producing_provenance, x$runtime))
  check("all three authoritative aggregate manifests use explicit coordinate contracts",
    identical(manifest$tier, c("tier1", "tier1_meta", "tier2")) &&
      all(manifest$build == "GRCh38") && all(manifest$coordinate_system == "one_based") &&
      all(manifest$source_format == "CHROM:POS:REF:ALT"))
  check("operator H1 and optional exposure settings are preserved",
    identical(op$h1_bayes, x$op$h1_bayes) &&
      identical(op$exposure_by_cohort, x$op$exposure_by_cohort))
  check("operator QTL and BED overrides survive while live upstream map is reused",
    identical(op$gusto_human_context$qtl_tables, x$op$gusto_human_context$qtl_tables) &&
      identical(op$gusto_human_context$beds, x$op$gusto_human_context$beds) &&
      identical(op$gusto_human_context$phenotype_map, x$context$phenotype_map))
  check("setup reads headers without parsing aggregate sentinel rows",
    result$code == 0L && grepl("headers checked only", result$text, fixed = TRUE) &&
      identical(receipt$aggregate_validation, "headers only; stage row checks still required"))
  check("setup receipt separates exact merged groups from chromosome labels",
    receipt$metadata$merged_groups == 18L && receipt$chromosome_labelled_summary_entries == 1L &&
      identical(receipt$selection, "exact exposure == ppBMI; never strip chromosome suffixes"))
  check("receipt preserves full numeric exposure moments and descriptive source sizes",
    receipt$metadata$n == 417 && receipt$metadata$mean == as.numeric(x$e$mean) &&
      receipt$metadata$sd == as.numeric(x$e$sd) && all(unlist(receipt$source_bytes) > 0))
  check("complete declarations print generator command without claiming full readiness",
    grepl("Next (requires an empty destination)", result$text, fixed = TRUE) &&
      grepl("generator preflight is still required", result$text, fixed = TRUE))
  check("private setup outputs have mode 0600",
    all(bitwAnd(as.integer(file.info(outputs)$mode), 511L) == 384L))
  before <- lapply(outputs, raw); repeated <- invoke(x)
  check("identical rerun reuses all four outputs",
    repeated$code == 0L && identical(before, lapply(outputs, raw)))
  check("all source metadata and original operator files stay byte-identical",
    identical(originals, lapply(c(x$producer, x$runtime, x$summary_path, x$upstream, x$operator), raw)))
  writeLines(c(readLines(x$output), "# OPERATOR_EDIT"), x$output)
  edited <- lapply(outputs, raw)
  conflict <- invoke(x)
  check("edited output is never overwritten on rerun",
    conflict$code != 0L && grepl("Refusing to overwrite differing output", conflict$text, fixed = TRUE) &&
      identical(edited, lapply(outputs, raw)))
  pending <- fixture("template", existing = FALSE)
  result <- invoke(pending)
  copy <- yaml::read_yaml(pending$output)
  check("missing original operator uses template without an optional H1 annotation blocker",
    result$code == 0L && grepl("Optional H1 annotations unrecorded", result$text, fixed = TRUE) &&
      !grepl("Config generation remains blocked", result$text, fixed = TRUE) &&
      grepl("Next (requires", result$text, fixed = TRUE))
  check("template fallback never invents H1 evidence or QTL BED overrides",
    is.null(copy$h1_bayes$count_provenance) && is.null(copy$h1_bayes$producing_pipeline) &&
      is.null(copy$gusto_human_context$qtl_tables) && is.null(copy$gusto_human_context$beds) &&
      !file.exists(pending$operator))
  explicit <- fixture("explicit")
  chosen <- file.path(explicit$base, "chosen.yaml")
  chosen_op <- explicit$op; chosen_op$h1_bayes$seed <- 2468L; yaml::write_yaml(chosen_op, chosen)
  result <- invoke(explicit, chosen)
  check("explicit existing operator path is honored",
    result$code == 0L && yaml::read_yaml(explicit$output)$h1_bayes$seed == 2468L)
  cohorts <- c("GUSTO", "mouse_Babr_CAST", "mouse_FVB_CAST_E12",
    "maize_B73_Ki11", "maize_B73_W22", "fish_Pprolifica_AB")
  for (id in cohorts) {
    path <- file.path(explicit$project, "data", id, "private", paste0(id, ".baseline_nb.yaml"))
    if (id == "GUSTO") next
    cfg <- yaml::read_yaml(file.path(root, "config", paste0(id, ".example.yaml")))
    dir.create(dirname(path), recursive = TRUE); yaml::write_yaml(cfg, path)
  }
  destination <- file.path(explicit$base, "generated configs")
  generated <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(file.path(root, "bin/prepare_promoter_configs.R"), explicit$project,
      explicit$output, destination)), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(generated, "status"))) stop(paste(generated, collapse = "\n"))
  inventory <- read.delim(file.path(destination, "CONFIGS.tsv"))
  configs <- lapply(inventory$destination_config, yaml::read_yaml)
  check("automatic setup operator feeds real generator for all eighteen configs",
    nrow(inventory) == 18L && setequal(inventory$cohort_id, cohorts) &&
      setequal(inventory$radius, c(0L, 1000L, 3000L)))
  check("generated GUSTO configs retain original QTL BED overrides and new ppBMI manifests",
    all(vapply(configs[inventory$cohort_id == "GUSTO"], function(cfg)
      identical(cfg$regulatory$human_context$qtl_tables, explicit$op$gusto_human_context$qtl_tables) &&
        identical(cfg$regulatory$human_context$beds, explicit$op$gusto_human_context$beds) &&
        endsWith(cfg$regulatory$human_context$gxe_tables, "gxe_tables.ppBMI_v1.tsv"), TRUE)))
  check("generated reciprocal configs retain exact H1 declarations",
    all(vapply(configs[inventory$cohort_id != "GUSTO"], function(cfg)
      identical(cfg$h1_bayes, chosen_op$h1_bayes), TRUE)))
  reject <- function(label, mutate, expected) {
    z <- fixture(paste0("invalid-", length(checks)))
    mutate(z)
    result <- invoke(z)
    if (result$code == 0L || !grepl(expected, result$text, fixed = TRUE)) message(result$text)
    outdir <- file.path(z$project, "data/GUSTO/private/qtl_context_v1")
    check(label, result$code != 0L && grepl(expected, result$text, fixed = TRUE) &&
      !file.exists(z$output) && !dir.exists(outdir))
  }
  reject("missing producer stops before output writes",
    function(z) unlink(z$producer), "Missing or unreadable input")
  reject("disabled external ppBMI cannot be mistaken for local inclusion",
    function(z) {z$p$enabled[2] <- 0; tsv(z$p, z$producer)}, "External producer ppBMI is disabled")
  reject("duplicate ppBMI runtime rows are rejected",
    function(z) tsv(rbind(z$e, z$e), z$runtime), "exactly one ppBMI row")
  reject("producer runtime column mismatch is rejected",
    function(z) {z$e$column <- "wrong"; tsv(z$e, z$runtime)}, "column disagreement")
  reject("binary runtime cannot silently replace reviewed continuous source",
    function(z) {z$e$binary <- "True"; tsv(z$e, z$runtime)}, "Inconsistent runtime")
  reject("changed source transformation stops for review",
    function(z) {z$p$transform[2] <- "log"; tsv(z$p, z$producer)}, "Unreviewed ppBMI transformation")
  reject("changed source range stops for review",
    function(z) {z$p$min[2] <- 15; tsv(z$p, z$producer)}, "differs from supplied record: min")
  reject("changed runtime sample size stops for review",
    function(z) {z$e$n <- 416; tsv(z$e, z$runtime)}, "differs from supplied record: n")
  reject("nonfinite exposure moments are rejected",
    function(z) {z$e$sd <- "NaN"; tsv(z$e, z$runtime)}, "Invalid numeric metadata")
  reject("chromosome shard summaries cannot substitute for exact merged records",
    function(z) tsv(z$summary[z$summary$exposure != "ppBMI", ], z$summary_path),
    "No exact merged ppBMI summary rows")
  reject("duplicate merged ancestry modality groups are rejected",
    function(z) tsv(rbind(z$summary, z$summary[1, ]), z$summary_path),
    "duplicate merged ppBMI summary group")
  reject("missing aggregate file stops before output writes",
    function(z) unlink(file.path(z$gx, "aggregate/gxe_tier2_results.tsv.gz")), "Missing or unreadable input")
  reject("missing interaction column cannot pass header compatibility",
    function(z) z$write_gz("gxe_tier1_results.tsv.gz",
      setdiff(z$headers$gxe_tier1_results.tsv.gz, "b_gi")), "tier1 missing columns: b_gi")
  reject("duplicate aggregate header columns are rejected",
    function(z) z$write_gz("gxe_tier1_results.tsv.gz",
      c(z$headers$gxe_tier1_results.tsv.gz, "b_gi")), "Empty or duplicate TSV columns")
  reject("missing phenotype map is not silently recreated",
    function(z) unlink(z$context$phenotype_map), "Missing or unreadable input")
  reject("invalid operator mapping is rejected before output writes",
    function(z) writeLines("h1_bayes: [", z$operator), "Parser error")
  # A conflict in the final receipt must prevent creation of all earlier outputs.
  late <- fixture("late-conflict")
  outdir <- file.path(late$project, "data/GUSTO/private/qtl_context_v1"); dir.create(outdir)
  held <- file.path(outdir, "ppBMI_setup.ppBMI_v1.yaml"); writeLines("USER_CONTENT", held)
  result <- invoke(late)
  check("all output conflicts are checked before any earlier output is written",
    result$code != 0L && identical(readLines(held), "USER_CONTENT") &&
      identical(list.files(outdir), basename(held)) && !file.exists(late$output))
  if (length(a) >= 2L) {
    dir.create(a[2], recursive = TRUE, showWarnings = FALSE)
    write.table(data.frame(check = checks, status = "pass", evidence = "local software fixtures"),
      file.path(a[2], "external_ppbmi_setup_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  cat("External ppBMI setup checks: ", length(checks), "; no real cohort files or model fitting\n", sep = "")
}
main()
