# Additive H6/H7 input contracts. Never writes into baseline_nb_v1.
EX_VERSION <- "exposure_nb_v1"
EX_PFAS <- paste0(c("PFBA", "PFOA", "PFNA", "PFDA", "PFUnDA", "PFBS", "PFHxS", "PFOS"), "_mat")
ex_need <- function(x, cols, label = "table") {
  miss <- setdiff(cols, names(x))
  if (length(miss)) stop(label, " missing columns: ", paste(miss, collapse = ", "))
  invisible(x)
}
ex_unique <- function(x, cols, label = "table") {
  ex_need(x, cols, label)
  if (anyNA(x[, ..cols]) || anyDuplicated(x[, ..cols])) stop(label, " has missing/duplicate keys")
  invisible(x)
}
ex_text <- function(x) !is.na(x) & nzchar(trimws(as.character(x))) &
  !tolower(trimws(as.character(x))) %in% c("unknown", "tbd", "todo", "unspecified", "not_available")
ex_numeric <- function(x, nonnegative = FALSE, allow_na = FALSE) {
  is.numeric(x) && all((allow_na & is.na(x)) | (is.finite(x) & (!nonnegative | x >= 0)))
}
ex_read <- function(path) {
  if (is.null(path) || !file.exists(path) || file.info(path)$size == 0) stop("Missing/empty input: ", path)
  data.table::fread(path, na.strings = c("", "NA"), showProgress = FALSE)
}
ex_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = if (grepl("\\.gz$", path)) ".tsv.gz" else ".tsv")
  data.table::fwrite(data.table::as.data.table(x), tmp, sep = "\t")
  if (!file.rename(tmp, path)) stop("Atomic table write failed: ", path)
  invisible(path)
}
ex_save <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path))
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Atomic RDS write failed: ", path)
  invisible(path)
}
ex_path <- function(s, ...) file.path(s$output_dir, ...)
ex_settings <- function(cfg) {
  x <- cfg$exposure_nb
  if (is.null(x)) stop("Add an exposure_nb section to a full schema-2 cohort YAML")
  required <- c("output_dir", "baseline_dir", "sample_crosswalk", "te_loci", "te_counts", "quant_qc")
  ex_need(x, required, "exposure_nb")
  base <- dirname(attr(cfg, "file") %||ex% file.path(getwd(), "config.yaml"))
  for (nm in intersect(c(required, "pfas", "exposure_dictionary", "libraries", "repeat_intervals",
                        "annotation_manifest", "gene_sets", "orthology", "h5_annotations"), names(x))) {
    p <- x[[nm]]
    if (!is.null(p) && nzchar(p)) x[[nm]] <- if (grepl("^/", p)) p else file.path(base, p)
  }
  defaults <- list(batch_jobs = 25L, seed = 1147L, bootstrap = 2000L, bootstrap_top = 20L,
    bootstrap_min_valid = .9, permutations = 2000L, min_cpm = 1, min_prevalence = .2,
    min_expressed = 3L, local_flank = 1000L, family_only = FALSE, common_pfas_samples = FALSE,
    modification = FALSE, max_condition = 1e8, quantifier_version = "1.0.3.1",
    quantifier_commit = "26a6c38df6d143387799edc2f3e9f976505f093d",
    star_version = "2.7.11b", multimap_limit = 100L, annotation_min_overlap = .2,
    count_column = "count", allow_unverified_independence = FALSE)
  for (nm in names(defaults)) if (is.null(x[[nm]])) x[[nm]] <- defaults[[nm]]
  for (nm in c("batch_jobs", "bootstrap", "bootstrap_top", "permutations", "min_expressed", "local_flank"))
    if (length(x[[nm]]) != 1 || !is.finite(x[[nm]]) || x[[nm]] < 1 || x[[nm]] != floor(x[[nm]]))
      stop("Invalid exposure_nb.", nm)
  if (x$bootstrap < 200L || x$bootstrap_min_valid <= 0 || x$bootstrap_min_valid > 1 ||
      x$min_prevalence <= 0 || x$min_prevalence > 1 || x$min_cpm < 0 ||
      x$max_condition <= 1 || x$multimap_limit < 2) stop("Invalid exposure settings")
  b <- normalizePath(x$baseline_dir, mustWork = FALSE)
  o <- normalizePath(x$output_dir, mustWork = FALSE)
  if (identical(b, o) || startsWith(paste0(o, "/"), paste0(b, "/")))
    stop("Exposure output must be separate from baseline_dir")
  x$profile <- cfg$profile; x$cohort_id <- cfg$cohort_id; x$species <- cfg$species
  x$build <- cfg$reference$build; x$version <- EX_VERSION
  x$dosage <- cfg$dosage
  x
}
`%||ex%` <- function(x, y) if (is.null(x)) y else x
ex_require_human <- function(cfg) {
  if (!identical(cfg$cohort_id, "GUSTO") || !identical(cfg$profile, "human_trio"))
    stop("PFAS/mediation is restricted to GUSTO")
}
ex_validate_crosswalk <- function(d, cfg) {
  d <- data.table::copy(d)
  ex_need(d, c("sample_id", "unit_id", "cluster_id", "sex", "ancestry", "independence_verified"), "crosswalk")
  ex_unique(d, "sample_id", "crosswalk")
  if (any(!ex_text(d$unit_id) | !ex_text(d$cluster_id))) stop("Missing biological identifiers")
  ex_unique(d, "unit_id", "one library per biological unit")
  if (!is.logical(d$independence_verified) || anyNA(d$independence_verified)) stop("Explicit logical independence_verified required")
  if (cfg$profile == "human_trio") {
    ex_need(d, "SubjectID", "human crosswalk")
    ex_unique(d, "SubjectID", "human crosswalk")
    if (any(!ex_text(d$sex) | !ex_text(d$ancestry))) stop("Missing human sex/ancestry")
  } else {
    ex_need(d, "cross")
    if (anyNA(d$cross) || any(!d$cross %in% c("focal_mother", "alternate_mother"))) stop("Invalid reciprocal cross")
  }
  d
}
ex_validate_dictionary <- function(d) {
  cols <- c("exposure", "biospecimen", "collection_time", "original_units", "pre_z_transform",
    "z_reference", "z_center", "z_scale", "detection_limit_handling", "prior_imputation")
  ex_need(d, cols, "exposure dictionary"); ex_unique(d, "exposure", "exposure dictionary")
  if (!setequal(d$exposure, EX_PFAS)) stop("Dictionary must declare exactly the maternal eight PFAS")
  for (nm in setdiff(cols, c("z_center", "z_scale")))
    if (any(!ex_text(d[[nm]]))) stop("Undeclared exposure metadata: ", nm)
  if (!ex_numeric(d$z_center) || !ex_numeric(d$z_scale) || any(d$z_scale <= 0))
    stop("Finite pre-standardization z_center and positive z_scale required")
  if (any(!grepl("maternal", d$biospecimen, ignore.case = TRUE))) stop("Maternal biospecimen declaration required")
  invisible(d)
}
ex_join_pfas <- function(crosswalk, pf, dictionary, cfg, common = FALSE) {
  ex_require_human(cfg); ex_validate_dictionary(dictionary)
  ex_need(pf, c("SubjectID", EX_PFAS), "PFAS")
  ex_unique(pf, "SubjectID", "PFAS")
  for (nm in EX_PFAS) if (!ex_numeric(pf[[nm]], allow_na = TRUE)) stop("PFAS must be numeric, finite or NA: ", nm)
  out <- merge(crosswalk, pf[, c("SubjectID", EX_PFAS), with = FALSE], by = "SubjectID", all.x = TRUE, sort = FALSE)
  if (common) out[!complete.cases(out[, ..EX_PFAS]), (EX_PFAS) := NA_real_]
  out
}
ex_validate_loci <- function(loci, cfg) {
  ex_need(loci, c("locus_id", "species", "build", "chrom", "start", "end", "strand",
                  "family_id", "repeat_name", "repeat_class", "source"), "TE loci")
  ex_unique(loci, "locus_id", "TE loci")
  if (any(loci$species != cfg$species) || any(loci$build != cfg$reference$build)) stop("TE species/build mismatch")
  if (!ex_numeric(loci$start) || !ex_numeric(loci$end) || any(loci$start < 1 | loci$end < loci$start) ||
      any(loci$start != floor(loci$start) | loci$end != floor(loci$end))) stop("TE coordinates must be 1-based inclusive integers")
  if (any(!loci$strand %in% c("+", "-")) || anyNA(loci$strand) ||
      any(!ex_text(loci$family_id) | !ex_text(loci$chrom) | !ex_text(loci$source))) stop("Invalid TE strand/family/contig/source")
  if (any(vapply(c(";", "\t", "\n"), function(v) any(grepl(v, loci$family_id, fixed = TRUE)), logical(1))))
    stop("One native family per locus is required")
  invisible(loci)
}
ex_validate_quant <- function(counts, qc, loci, crosswalk, cfg, s) {
  ex_need(qc, c("sample_id", "mapped_fragments", "status", "species", "build", "layout",
    "strandedness", "quantifier", "quantifier_version", "multimap_limit", "annotation_id",
    "count_scope"), "quantification QC")
  ex_unique(qc, "sample_id", "quantification QC")
  if (any(qc$species != cfg$species | qc$build != cfg$reference$build)) stop("Quantification build/species mismatch")
  if (any(qc$status != "complete") || !ex_numeric(qc$mapped_fragments) || any(qc$mapped_fragments <= 0))
    stop("Only complete, positively normalized libraries can supply TE expression")
  if (any(!qc$layout %in% c("SE", "PE")) || any(!qc$strandedness %in% c("unstranded", "forward", "reverse")))
    stop("Unknown library layout/strandedness")
  if (any(!ex_text(qc$annotation_id))) stop("Explicit TE annotation ID/release is required")
  if (data.table::uniqueN(qc$annotation_id) != 1L) stop("Samples declare different TE annotations")
  if (!is.null(s$annotation_id) && !identical(unique(qc$annotation_id), s$annotation_id))
    stop("QC annotation ID differs from configured TE annotation ID")
  if (any(qc$count_scope != "all_catalogue_loci_sparse_zero") || any(qc$multimap_limit < 2) ||
      any(!ex_text(qc$quantifier) | !ex_text(qc$quantifier_version))) stop("Incomplete quantification scope/provenance")
  key <- if (isTRUE(s$family_only)) "family_id" else "locus_id"
  ex_need(counts, c("sample_id", key, "count"), "TE counts")
  ex_unique(counts, c("sample_id", key), "TE counts")
  if (!ex_numeric(counts$count, nonnegative = TRUE)) stop("Invalid TE counts")
  if ("unique_count" %in% names(counts) && !ex_numeric(counts$unique_count, nonnegative = TRUE))
    stop("Invalid unique counts")
  if ("unique_count" %in% names(counts)) {
    ex_need(qc, "unique_count_definition", "unique-mapper QC")
    if (anyNA(qc$unique_count_definition) || any(qc$unique_count_definition != "genomic_NH1"))
      stop("unique_count must come from strict genomic NH=1 alignments, not annotation assignment")
  }
  if (!all(counts$sample_id %in% qc$sample_id)) stop("Counts without completed library QC")
  if (!all(counts[[key]] %in% loci[[key]])) stop("Counts contain unannotated TE IDs")
  if (!all(qc$sample_id %in% crosswalk$sample_id)) stop("Quantification sample absent from explicit crosswalk")
  invisible(TRUE)
}
ex_run_lock <- function(s, cfg, inputs) {
  identity <- list(version = EX_VERSION, identity_schema = "readable_jobs_v2",
    settings = s,
    dosage = cfg$dosage, files = lapply(inputs, function(p) {
      i <- file.info(p); list(path = p, size = i$size, mtime = as.numeric(i$mtime))
    }))
  p <- ex_path(s, "private", "run_identity.rds")
  if (file.exists(p) && !identical(readRDS(p)$identity, identity))
    stop("Inputs/settings changed. Use a new exposure output_dir.")
  if (!file.exists(p)) ex_save(list(identity = identity), p)
  invisible(identity)
}
