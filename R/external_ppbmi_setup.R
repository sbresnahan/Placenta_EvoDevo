# Read small producer records and aggregate headers; never scan aggregate rows.
setup_external_ppbmi <- function(root, project, source_paths, existing_operator = NULL) {
  project <- normalizePath(project, mustWork = TRUE)
  readable <- function(path) {
    if (length(path) != 1L || is.na(path) || !file.exists(path) ||
        dir.exists(path) || file.access(path, 4L) != 0L)
      stop("Missing or unreadable input: ", paste(path, collapse = ", "))
    normalizePath(path, mustWork = TRUE)
  }
  mapping <- function(x, label) {
    if (!is.list(x) || !length(x) || is.null(names(x)) || anyNA(names(x)) ||
        any(!nzchar(names(x))) || anyDuplicated(names(x)))
      stop(label, " must be a nonempty mapping")
    x
  }
  read_yaml <- function(path) mapping(yaml::read_yaml(readable(path)), path)
  paths <- read_yaml(source_paths)
  for (key in c("gxe_root", "producing_config", "producer_revision"))
    if (!is.character(paths[[key]]) || length(paths[[key]]) != 1L ||
        is.na(paths[[key]]) || !nzchar(paths[[key]])) stop("Missing source setting: ", key)
  if (!grepl("^/", paths$gxe_root) || !grepl("^/", paths$producing_config))
    stop("External source paths must be absolute")
  expected <- mapping(paths$expected, "expected metadata")
  expected_keys <- c("column", "exposure_type", "transformation", "min", "max",
    "n", "n_missing", "mean", "sd", "merged_groups")
  if (!setequal(names(expected), expected_keys) ||
      any(vapply(expected, function(x) length(x) != 1L || is.na(x), TRUE)))
    stop("Source settings require all expected metadata fields")
  header <- function(path) {
    path <- readable(path)
    con <- if (endsWith(path, ".gz")) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con))
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) != 1L || !nzchar(line) || endsWith(line, "\t"))
      stop("Empty or invalid TSV header: ", path)
    columns <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (any(!nzchar(columns)) || anyDuplicated(columns))
      stop("Empty or duplicate TSV columns: ", path)
    columns
  }
  require_cols <- function(columns, keys, label) {
    missing <- setdiff(keys, columns)
    if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
  }
  small <- function(path, keys) {
    # Bound input sizes to guard against accidentally supplying a result table.
    path <- readable(path)
    if (file.info(path)$size > 5e6) stop("Expected a small metadata TSV: ", path)
    require_cols(header(path), keys, path)
    read.delim(path, colClasses = "character", check.names = FALSE,
      stringsAsFactors = FALSE, quote = "", comment.char = "", na.strings = c("", "NA"))
  }
  one <- function(x, key, label) {
    z <- x[!is.na(x[[key]]) & x[[key]] == "ppBMI", , drop = FALSE]
    if (nrow(z) != 1L) stop(label, " must have exactly one ppBMI row")
    z
  }
  number <- function(x, label, integer = FALSE, optional = FALSE) {
    if (optional && length(x) == 1L && is.na(x)) return(NA_real_)
    v <- suppressWarnings(as.numeric(x))
    if (length(v) != 1L || !is.finite(v) || (integer && (v < 0 || v != floor(v))))
      stop("Invalid numeric metadata: ", label)
    v
  }
  boolean <- function(x, label) {
    v <- tolower(as.character(x))
    if (length(v) != 1L || is.na(v) || !v %in% c("true", "false", "1", "0"))
      stop("Invalid boolean metadata: ", label)
    v %in% c("true", "1")
  }
  producer <- readable(paths$producing_config)
  runtime <- readable(file.path(paths$gxe_root, "inputs/exposure_manifest.tsv"))
  summary_path <- readable(file.path(paths$gxe_root, "aggregate/gxe_summary.tsv"))
  p <- one(small(producer, c("exposure_id", "column", "transform", "min", "max", "enabled")),
    "exposure_id", "Producer configuration")
  e <- one(small(runtime, c("exposure_id", "column", "binary", "n", "n_missing", "mean", "sd")),
    "exposure_id", "Runtime exposure summary")
  if (!boolean(p$enabled, "producer.enabled")) stop("External producer ppBMI is disabled")
  if (is.na(p$column) || is.na(e$column) || p$column != e$column)
    stop("Producer/runtime ppBMI column disagreement")
  binary <- boolean(e$binary, "runtime.binary")
  type <- if (binary) "binary" else "continuous"
  transform <- p$transform
  # This producer's untransformed ppBMI contract is the reviewed setup target.
  if (is.na(transform) || transform != "none")
    stop("Unreviewed ppBMI transformation; expected none")
  lo <- number(p$min, "producer.min")
  hi <- number(p$max, "producer.max")
  if (lo > hi) stop("Invalid producer range")
  n <- number(e$n, "runtime.n", TRUE)
  missing <- number(e$n_missing, "runtime.n_missing", TRUE)
  mean <- number(e$mean, "runtime.mean")
  sd <- number(e$sd, "runtime.sd")
  if (n < 2 || sd <= 0 || mean < lo || mean > hi ||
      (binary && (mean < 0 || mean > 1)))
    stop("Inconsistent runtime exposure moments or range")
  summary <- small(summary_path, c("ancestry", "modality", "exposure", "n_tested", "n_sig"))
  merged <- summary[!is.na(summary$exposure) & summary$exposure == "ppBMI", , drop = FALSE]
  if (!nrow(merged)) stop("No exact merged ppBMI summary rows; chromosome shards cannot substitute")
  if (anyNA(merged[c("ancestry", "modality")]) ||
      any(!nzchar(merged$ancestry)) || any(!nzchar(merged$modality)) ||
      anyDuplicated(merged[c("ancestry", "modality")]))
    stop("Missing or duplicate merged ppBMI summary group")
  for (i in seq_len(nrow(merged))) {
    nt <- number(merged$n_tested[i], "summary.n_tested", TRUE)
    ns <- number(merged$n_sig[i], "summary.n_sig", TRUE)
    if (ns > nt) stop("Summary n_sig exceeds n_tested")
  }
  shard_count <- sum(grepl("^ppBMI[.]chr", summary$exposure), na.rm = TRUE)
  actual <- list(column = p$column, exposure_type = type, transformation = transform,
    min = lo, max = hi, n = n, n_missing = missing, mean = mean, sd = sd,
    merged_groups = nrow(merged))
  for (key in expected_keys) {
    same <- if (is.numeric(actual[[key]])) isTRUE(all.equal(actual[[key]],
      as.numeric(expected[[key]]), tolerance = 1e-12)) else identical(actual[[key]], expected[[key]])
    if (!same) stop("Live ppBMI metadata differs from supplied record: ", key,
      ". Review source records and config/external_ppbmi.yaml before setup.")
  }
  basenames <- c(tier1 = "gxe_tier1_results.tsv.gz", tier1_meta = "gxe_tier1_meta.tsv.gz",
    tier2 = "gxe_tier2_results.tsv.gz")
  required <- list(
    tier1 = c("ancestry", "modality", "exposure", "phenotype_id", "variant_id",
      "b_gi", "b_gi_se", "b_g", "b_e", "pval_nominal", "pval_beta", "qval"),
    tier1_meta = c("modality", "exposure", "phenotype_id", "beta_meta", "se_meta",
      "meta_variant_id", "scope", "lead_concordant", "p_feature_acat", "q_feature_acat", "n_ancestries"),
    tier2 = c("modality", "exposure", "phenotype_id", "beta_meta", "se_meta",
      "variant_id", "p_meta", "q_meta"))
  aggregates <- setNames(vapply(basenames, function(b)
    readable(file.path(paths$gxe_root, "aggregate", b)), ""), names(basenames))
  headers <- lapply(aggregates, header)
  for (tier in names(required)) require_cols(headers[[tier]], required[[tier]], tier)
  upstream <- readable(file.path(project, "data/GUSTO/private/GUSTO.baseline_nb.yaml"))
  cfg <- read_yaml(upstream)
  if (!identical(cfg$cohort_id, "GUSTO")) stop("Upstream configuration must identify GUSTO")
  context <- mapping(cfg$regulatory$human_context, "GUSTO human_context")
  inherited <- file.path(project, "data/private/promoter_operator_settings.yaml")
  if (is.null(existing_operator)) existing_operator <- if (file.exists(inherited))
    inherited else file.path(root, "config/promoter_operator_settings.example.yaml")
  operator_source <- readable(existing_operator)
  settings <- read_yaml(operator_source)
  if (length(setdiff(names(settings), c("h1_bayes", "gusto_human_context", "exposure_by_cohort"))))
    stop("Unknown operator settings field")
  if (!is.null(settings$gusto_human_context))
    mapping(settings$gusto_human_context, "operator.gusto_human_context")
  # Use the live upstream map, retaining QTL/BED overrides and all other settings.
  map <- readable(context$phenotype_map)
  if (!grepl("^/", context$phenotype_map)) stop("GUSTO phenotype_map must be absolute")
  outdir <- file.path(project, "data/GUSTO/private/qtl_context_v1")
  destinations <- c(
    gxe = file.path(outdir, "gxe_tables.ppBMI_v1.tsv"),
    exposure = file.path(outdir, "ppBMI_exposure.ppBMI_v1.tsv"),
    operator = file.path(project, "data/private/promoter_operator_settings.ppBMI_v1.yaml"),
    receipt = file.path(outdir, "ppBMI_setup.ppBMI_v1.yaml"))
  settings$gusto_human_context$phenotype_map <- map
  settings$gusto_human_context$gxe_tables <- destinations[["gxe"]]
  settings$gusto_human_context$exposure_manifest <- destinations[["exposure"]]
  format_num <- function(x) format(x, digits = 17L, trim = TRUE, scientific = FALSE)
  scale <- if (binary) "0/1; no z-score; scanner centers the binary exposure" else paste0(
    "z=(x-mean)/sd pooled across retained ancestries; population SD (ddof=0); mean=",
    format_num(mean), "; sd=", format_num(sd), "; nonmissing n=", n)
  exposure <- data.frame(exposure = "ppBMI", exposure_type = type,
    transformation = transform, scale_definition = scale,
    range_definition = paste0("source inclusion: ", format_num(lo), " <= ppBMI <= ",
      format_num(hi), "; before pooled scaling"),
    producing_config = producer, producing_provenance = runtime, enabled = TRUE)
  manifest <- data.frame(input_id = paste0("ppBMI_", names(aggregates)),
    path = unname(aggregates), build = "GRCh38", tier = names(aggregates),
    source_format = "CHROM:POS:REF:ALT", coordinate_system = "one_based")
  pending <- character()
  placeholders <- function(x, field) {
    if (is.list(x)) {
      labels <- names(x); if (is.null(labels)) labels <- seq_along(x)
      for (i in seq_along(x)) placeholders(x[[i]], paste(field, labels[i], sep = "."))
    } else if (is.character(x) && any(grepl(
        "REVIEW_REQUIRED|DECLARE_|TODO|UNSPECIFIED|/actual/path|/absolute/path", x,
        ignore.case = TRUE), na.rm = TRUE)) pending <<- c(pending, field)
  }
  inspected <- settings
  inspected$h1_bayes <- h1_normalize_annotations(settings$h1_bayes, "operator.h1_bayes")
  h1_validate_controls(inspected$h1_bayes, "operator.h1_bayes")
  placeholders(inspected, "operator")
  unrecorded <- setdiff(c("count_provenance", "producing_pipeline"), names(inspected$h1_bayes))
  pending <- unique(pending)
  small_paths <- c(source_paths = readable(source_paths), producer = producer,
    runtime = runtime, summary = summary_path, upstream = upstream, operator = operator_source)
  receipt <- list(setup = "external_ppBMI_metadata_only_v1",
    local_enabled_meaning = "include existing external ppBMI results; not prior local analysis",
    producer_revision = paths$producer_revision,
    scaling_basis = "reviewed 50_build_gxe_inputs.py; continuous pooled z-score with ddof=0",
    metadata = actual, chromosome_labelled_summary_entries = shard_count,
    selection = "exact exposure == ppBMI; never strip chromosome suffixes",
    source_paths = as.list(small_paths), source_bytes = as.list(unname(file.info(small_paths)$size)),
    aggregate_headers = headers, aggregate_validation = "headers only; stage row checks still required",
    outputs = as.list(destinations[c("gxe", "exposure", "operator")]),
    pending_operator_fields = pending, optional_unrecorded_h1_annotations = unrecorded)
  names(receipt$source_bytes) <- names(small_paths)
  staging <- tempfile("ppBMI setup "); dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  staged <- setNames(file.path(staging, basename(destinations)), names(destinations))
  for (key in c("gxe", "exposure")) write.table(get(if (key == "gxe") "manifest" else "exposure"),
    staged[[key]], sep = "\t", quote = FALSE, row.names = FALSE)
  yaml::write_yaml(settings, staged[["operator"]], precision = 17L)
  yaml::write_yaml(receipt, staged[["receipt"]], precision = 17L)
  bytes <- function(p) readBin(p, "raw", n = file.info(p)$size)
  # Check every output before making directories or copying any output.
  for (key in names(destinations)) {
    dest <- destinations[[key]]
    link <- Sys.readlink(dest)
    if ((!is.na(link) && nzchar(link)) || dir.exists(dest))
      stop("Refusing non-regular output: ", dest)
    if (file.exists(dest) && !identical(bytes(dest), bytes(staged[[key]])))
      stop("Refusing to overwrite differing output: ", dest,
        ". Keep existing files and review the versioned setup inputs.")
  }
  mask <- Sys.umask("0077"); on.exit(Sys.umask(mask), add = TRUE)
  for (key in names(destinations)) {
    dest <- destinations[[key]]
    if (file.exists(dest)) next
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(staged[[key]], dest, overwrite = FALSE, copy.mode = FALSE))
      stop("Cannot create setup output: ", dest)
    Sys.chmod(dest, "0600")
    if (!identical(bytes(dest), bytes(staged[[key]]))) stop("Output verification failed: ", dest)
  }
  cat("External ppBMI setup complete: ", nrow(merged), " exact merged summary groups; ",
    shard_count, " chromosome-labelled entries excluded by exact-label selection.\n", sep = "")
  cat("Aggregate headers checked only; context-stage row and mapping checks remain required.\n")
  cat("Operator copy: ", destinations[["operator"]], "\n", sep = "")
  if (length(unrecorded)) cat("Optional H1 annotations unrecorded (not a config blocker): ",
    paste(unrecorded, collapse = ", "), ". No new count or pipeline audit is required.\n", sep = "")
  if (length(pending)) {
    cat("Config generation remains blocked by operator declarations:\n",
      paste0("  ", pending, collapse = "\n"), "\n", sep = "")
    cat("Review these fields in the versioned operator copy; H1 provenance is not inferred from ppBMI.\n")
  } else {
    cat("No unresolved declaration placeholders found; generator preflight is still required.\n")
    cat("Next (requires an empty destination): bash bin/Rscript_sif bin/prepare_promoter_configs.R ",
      shQuote(project), " ", shQuote(destinations[["operator"]]), " ",
      shQuote(file.path(project, "data/private/promoter_configs_ppBMI_v1")), "\n", sep = "")
  }
  invisible(list(outputs = destinations, metadata = actual, pending = pending))
}
