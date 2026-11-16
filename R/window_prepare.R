# Reuse of the completed upstream POE layer and immutable count pointers.
nb_poe_contract <- function(cfg, s) {
  keep_paths <- c("libraries", "units", "tx_map", "event_map", "exons", "events",
    "markers", "trios")
  list(cohort_id = cfg$cohort_id, profile = cfg$profile, species = cfg$species,
    reference = cfg$reference, dosage = cfg$dosage, parents = cfg$parents,
    paths = cfg$paths[intersect(keep_paths, names(cfg$paths))],
    preprocessing = cfg[c("fragments", "ornaments", "phase_origin", "gusto_precomputed", "models", "markers")],
    settings = s[c("min_snps", "min_total_counts", "min_samples", "min_per_cross",
      "metadata", "sex_column", "ancestry_column")])
}
nb_prepare_window <- function(cfg, s, source_cfg, root) {
  upstream <- nb_settings(source_cfg)
  nb_require(upstream, "prepare"); nb_require(upstream, "combine")
  if (!identical(nb_poe_contract(cfg, s), nb_poe_contract(source_cfg, upstream)))
    stop("POE reuse requires matching cohort, reference, dosage, preprocessing and model settings")
  if (dir.exists(s$output_dir) && length(list.files(s$output_dir, all.files = TRUE, no.. = TRUE)))
    stop("prepare-window requires a fresh, empty destination")
  if (identical(normalizePath(s$output_dir, mustWork = FALSE),
      normalizePath(upstream$output_dir))) stop("Source and destination must differ")
  allow <- c(NB_POE_LAYER_FILES, "private/metadata.rds",
    unlist(lapply(c("gene", "transcript"), function(level)
      paste0(level, c("_batches.tsv", "_feature_batches.tsv")))))
  allow <- unique(allow)
  src <- nb_path(upstream, allow)
  if (any(!file.exists(src) | file.info(src)$size <= 0)) stop("Incomplete allowlisted POE inputs")
  for (level in c("gene", "transcript")) {
    manifest <- reg_read(nb_path(upstream, paste0(level, "_batches.tsv")))
    require_columns(manifest, c("input_file", "batch"))
    if (any(!grepl("^/", manifest$input_file))) stop("Count-batch pointers must be absolute")
    if (any(!file.exists(manifest$input_file) | file.info(manifest$input_file)$size <= 0))
      stop("Immutable private count batches unavailable")
  }
  for (i in seq_along(allow)) {
    dest <- nb_path(s, allow[i]); dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(src[i], dest)) stop("POE metadata copy failed")
  }
  nb_write(data.table(cohort_id = cfg$cohort_id,
    source_config = attr(source_cfg, "file"), source_output_dir = upstream$output_dir,
    destination_config = attr(cfg, "file"), te_tss_window_bp = reg_settings(cfg)$te_tss_window_bp,
    reused_layer = "completed_gene_transcript_POE_and_immutable_count_batches",
    created_utc = as.character(Sys.time())), nb_path(s, "window_preparation.tsv"))
  key <- list(operation = "prepare-window", source_prepare = nb_require(upstream, "prepare")$key,
    source_combine = nb_require(upstream, "combine")$key,
    contract = nb_poe_contract(cfg, s))
  nb_finish(s, "prepare", key, c(nb_path(s, allow), nb_path(s, "window_preparation.tsv")))
  nb_finish(s, "combine", key, nb_path(s, c("gene_POE.tsv", "transcript_POE.tsv")))
  invisible(TRUE)
}
