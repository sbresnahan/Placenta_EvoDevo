# Independent aggregate snapshots only. No private model or baseline dependency.
EX_ARCHIVE_CORE <- c("model_results.tsv", "model_coverage.tsv", "job_inventory.tsv", "provenance.tsv")
EX_ARCHIVE_OPTIONAL <- c("pilot_resources.tsv", "gene_evidence.tsv", "functional_enrichment.tsv",
  "H5_support_enrichment.tsv", "frozen_annotation_manifest.tsv", "bootstrap_selection.tsv",
  "report_summary.tsv", "TE_scaling.tsv", "exposure_coverage.tsv", "target_family_links.tsv",
  "target_insertion_links.tsv", "H5_ppBMI_support_enrichment.tsv")
ex_archive_safe <- function(x) {
  is.character(x) && length(x) > 0L && all(!is.na(x) & nzchar(x) &
    !grepl("(^/|\\\\|[\r\n\t]|(^|/)\\.{1,2}(/|$)|//|/$|:)", x))
}
ex_archive_file <- function(path) {
  p <- path
  repeat {
    target <- Sys.readlink(p)
    if (!is.na(target) && nzchar(target)) stop("Symlinked snapshot path rejected: ", basename(path))
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }
  if (!file.exists(path) || dir.exists(path) || !file_test("-f", path) ||
      !is.finite(file.info(path)$size) || file.info(path)$size <= 0)
    stop("Missing/empty/nonregular snapshot file: ", basename(path))
  invisible(path)
}
ex_archive_table <- function(path, required = character()) {
  ex_archive_file(path)
  d <- data.table::fread(path, showProgress = FALSE, check.names = FALSE)
  if (!length(names(d)) || anyDuplicated(names(d)) || any(!nzchar(names(d))))
    stop("Invalid table header: ", basename(path))
  if (!all(required %in% names(d))) stop("Missing required columns in ", basename(path),
    ": ", paste(setdiff(required, names(d)), collapse = ", "))
  d
}
ex_archive_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  id <- cfg$cohort_id
  if (!is.character(id) || length(id) != 1L || is.na(id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", id)) stop("Unsafe or missing cohort_id")
  p <- cfg$exposure_nb$output_dir
  if (!is.character(p) || length(p) != 1L || is.na(p) || !nzchar(p) ||
      grepl("[\r\n\t]|\\$|~", p)) stop("Declare a literal exposure_nb.output_dir")
  if (!grepl("^/", p)) p <- file.path(dirname(normalizePath(path, mustWork = TRUE)), p)
  # Do not resolve away symlinks before screening their components.
  list(id = id, snapshot = file.path(p, "aggregate_export"))
}
ex_archive_snapshot <- function(cfg) {
  p <- cfg$snapshot
  m <- ex_archive_table(file.path(p, "MANIFEST.tsv"), c("file", "bytes"))
  # Retain allowlisted descriptive metadata columns.
  m <- m[, .(file = as.character(file), bytes = suppressWarnings(as.numeric(bytes)))]
  if (!ex_archive_safe(m$file) || anyDuplicated(m$file) ||
      any(!is.finite(m$bytes) | m$bytes <= 0 | m$bytes != floor(m$bytes)))
    stop("Invalid or duplicate snapshot manifest paths/sizes")
  allowed <- m$file %in% c(EX_ARCHIVE_CORE, EX_ARCHIVE_OPTIONAL) |
    grepl("^bootstrap/[A-Za-z0-9_.-]+[.]tsv$", m$file)
  if (any(!allowed)) stop("Unapproved snapshot file: ", paste(m$file[!allowed], collapse = ", "))
  if (length(setdiff(EX_ARCHIVE_CORE, m$file))) stop("Snapshot missing core aggregate files")
  schemas <- list(model_results.tsv = c("job_id", "module", "cohort_id", "level",
      "feature_id", "family_id", "status"),
    model_coverage.tsv = "module", job_inventory.tsv = c("module", "jobs"),
    provenance.tsv = c("version", "cohort_id", "species", "build"))
  for (i in seq_len(nrow(m))) {
    f <- m$file[i]; src <- file.path(p, f)
    d <- ex_archive_table(src, schemas[[f]])
    if (!identical(as.numeric(file.info(src)$size), m$bytes[i])) stop("Snapshot byte size mismatch: ", f)
    # family_id is the native TE family, not a human participant-family field.
    banned <- grepl(paste0("sample_id|subject|cluster_id|unit_id|child_id|mother_id|father_id|",
      "infant_id|run_id|library_name|sample_accession|input_file|source_path|",
      "warning|blocker|^(e|m)$|family_count|maternal_count|paternal_count|te_abundance|",
      "^pf(ba|oa|na|da|unda|bs|hxs|os)_mat$"), names(d), ignore.case = TRUE)
    if (any(banned)) stop("Private columns rejected in ", f, ": ", paste(names(d)[banned], collapse = ", "))
    if ("cohort_id" %in% names(d) && nrow(d) &&
        any(is.na(d$cohort_id) | d$cohort_id != cfg$id)) stop("Cohort mismatch in ", f)
    if (f == "provenance.tsv" && (!nrow(d) || any(d$version != "exposure_nb_v1")))
      stop("Unsupported or empty exposure provenance")
    if (f == "model_results.tsv" && (anyNA(d$job_id) || anyDuplicated(d$job_id)))
      stop("Missing/duplicate result job IDs")
  }
  m
}
ex_archive <- function(configs, out, root) {
  if (!grepl("[.]tar[.]gz$", out) || file.exists(out)) stop("Invalid or existing archive output")
  out <- file.path(normalizePath(dirname(out), mustWork = TRUE), basename(out))
  name <- sub("[.]tar[.]gz$", "", basename(out))
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", name)) stop("Unsafe archive basename")
  cfgs <- lapply(configs, ex_archive_config)
  if (!length(cfgs) || anyDuplicated(vapply(cfgs, `[[`, "", "id"))) stop("Missing/duplicate cohort IDs")
  inventories <- lapply(cfgs, ex_archive_snapshot)
  stage <- tempfile("exposure-archive-"); dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  top <- file.path(stage, name); dir.create(top)
  tmp <- paste0(out, ".partial")
  if (file.exists(tmp)) stop("Remove the incomplete archive staging directory first")
  on.exit(unlink(tmp), add = TRUE)
  availability <- list()
  for (i in seq_along(cfgs)) {
    cfg <- cfgs[[i]]; m <- inventories[[i]]
    base <- file.path(top, "cohorts", cfg$id, "exposure_nb_v1")
    dir.create(base, recursive = TRUE)
    for (f in m$file) {
      src <- file.path(cfg$snapshot, f); dest <- file.path(base, f)
      ex_archive_file(src)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(src, dest) || file.info(dest)$size != m[file == f]$bytes)
        stop("Snapshot copy failed or changed size: ", f)
      # Screen the actual staged table as well, avoiding unvalidated copy content.
    }
    data.table::fwrite(m, file.path(base, "MANIFEST.tsv"), sep = "\t")
    # Validate the staged copy again before packaging.
    ex_archive_snapshot(list(id = cfg$id,
      snapshot = base))
    optional <- c(EX_ARCHIVE_OPTIONAL, "bootstrap/*.tsv")
    found <- c(optional[-length(optional)] %in% m$file, any(startsWith(m$file, "bootstrap/")))
    availability[[i]] <- data.table::data.table(cohort_id = cfg$id, output = optional,
      state = ifelse(found, "present_in_snapshot", "absent_from_snapshot"),
      interpretation = "presence_is_not_success_absence_is_not_a_null_result")
  }
  dir.create(file.path(top, "docs"))
  for (f in c("docs/EXPOSURE_NB_WORKFLOW.md", "docs/H6_H7_MODELS.md", "docs/exposure_output_schema.md")) {
    src <- file.path(root, f)
    if (!file.exists(src) || !file.copy(src, file.path(top, "docs", basename(f))))
      stop("Missing/copy-failed archive documentation: ", f)
  }
  data.table::fwrite(data.table::rbindlist(availability), file.path(top, "AVAILABILITY.tsv"), sep = "\t")
  code <- sort(unique(c("bin/export_results.R", "poe_exposure.R",
    list.files(file.path(root, "R"), pattern = "[.]R$", full.names = FALSE))))
  code <- ifelse(grepl("/", code, fixed = TRUE), code, file.path("R", code))
  data.table::fwrite(data.table::data.table(path = code, bytes = file.info(file.path(root, code))$size),
    file.path(top, "SCRIPTS.tsv"), sep = "\t")
  writeLines(c("# Independent H6/H7 aggregate archive", "",
    paste("Created UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)), "",
    "Contains existing aggregate snapshots only; no model fitting or baseline report dependency.",
    "family_id is a native TE family identifier. Failed/unavailable rows are retained.",
    "Optional-output presence is not certification of success; absence is not a null finding.",
    "Inventories contain filenames and sizes only and do not establish content integrity or freshness.",
    "Review source text and cohort release rules; privacy screening is not anonymization approval.",
    "No private counts, model objects, sample crosswalks, PFAS values or raw bootstrap draws.",
    "See docs/EXPOSURE_NB_WORKFLOW.md for inputs, execution and aggregate export instructions."),
    file.path(top, "README.md"))
  files <- sort(list.files(top, recursive = TRUE, full.names = FALSE))
  data.table::fwrite(data.table::data.table(archive_path = files,
    bytes = file.info(file.path(top, files))$size), file.path(top, "MANIFEST.tsv"), sep = "\t")
  old <- setwd(stage); on.exit(setwd(old), add = TRUE)
  status <- utils::tar(tmp, files = name, compression = "gzip", tar = "internal")
  if ((length(status) && status != 0L) || !file.exists(tmp)) stop("Archive creation failed")
  entries <- utils::untar(tmp, list = TRUE)
  if (!all(file.path(name, c(files, "MANIFEST.tsv")) %in% entries)) stop("Archive listing incomplete")
  verify <- file.path(stage, "verification"); dir.create(verify)
  utils::untar(tmp, exdir = verify)
  vtop <- file.path(verify, name)
  inv <- data.table::fread(file.path(vtop, "MANIFEST.tsv"))
  if (any(!file.exists(file.path(vtop, inv$archive_path))) ||
      any(file.info(file.path(vtop, inv$archive_path))$size != inv$bytes)) stop("Archive extraction size mismatch")
  for (cfg in cfgs) ex_archive_snapshot(list(id = cfg$id,
    snapshot = file.path(vtop, "cohorts", cfg$id, "exposure_nb_v1")))
  if (file.exists(out) || !file.rename(tmp, out)) stop("Cannot finalize archive without overwrite")
  cat("Created ", basename(out), " (", file.info(out)$size, " bytes; ",
    length(cfgs), " independent cohort snapshots)\n", sep = "")
  invisible(out)
}
