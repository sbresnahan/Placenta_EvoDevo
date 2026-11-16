# Batch-wise H1 coverage; count rules remain in h1_counts.
h1_prep_number <- function(name, default = NULL, integer = FALSE) {
  x <- Sys.getenv(name, unset = "")
  if (!nzchar(x)) return(default)
  n <- suppressWarnings(as.numeric(x))
  if (length(n) != 1L || !is.finite(n) || n <= 0 ||
      (integer && (n != floor(n) || n > .Machine$integer.max)))
    stop(name, " must be a positive ", if (integer) "integer" else "number")
  n
}
h1_prep_resources <- function(manifest, annotation, map) {
  requested <- h1_prep_number("H1_PREP_WORKERS", 1L, integer = TRUE)
  budget <- h1_prep_number("H1_PREP_MEMORY_GB")
  workers <- min(requested, max(1L, nrow(manifest)))
  for (name in c("LSB_DJOB_NUMPROC", "LSB_MAX_NUM_PROCESSORS")) {
    slots <- h1_prep_number(name, integer = TRUE)
    if (!is.null(slots)) workers <- min(workers, slots)
  }
  sizes_known <- "object_bytes" %in% names(manifest) &&
    all(is.finite(manifest$object_bytes) & manifest$object_bytes > 0)
  base <- 1024^3 + 4 * as.numeric(object.size(annotation)) +
    2 * as.numeric(object.size(map)) + as.numeric(object.size(manifest))
  per_worker <- NA_real_
  if (nrow(manifest) && sizes_known) {
    per_worker <- 4 * max(manifest$object_bytes) + 256 * 1024^2
    if (!is.null(budget)) {
      capacity <- floor((budget * 1024^3 - base) / per_worker)
      if (capacity < 1) stop("H1_PREP_MEMORY_GB is below the estimated one-worker requirement (",
        sprintf("%.2f", (base + per_worker) / 1024^3), " GB)")
      workers <- min(workers, capacity)
    } else if (workers > 1L) {
      message("H1 preparation: no H1_PREP_MEMORY_GB supplied; using batch-wise serial loading")
      workers <- 1L
    }
  } else if (nrow(manifest)) {
    message("H1 preparation: batch object_bytes unavailable; using batch-wise serial loading; ",
      "memory use cannot be estimated from the manifest")
    workers <- 1L
  }
  if (.Platform$OS.type != "unix" && workers > 1L) {
    message("H1 preparation: process parallelism requires Unix; using one worker")
    workers <- 1L
  }
  message(sprintf("H1 preparation: requested workers=%d; effective workers=%d; batches=%d; memory budget=%s GB",
    requested, workers, nrow(manifest), if (is.null(budget)) "unset" else format(budget)))
  if (is.finite(per_worker)) message(sprintf(
    "H1 preparation: estimated coordinator=%.2f GB; each worker=%.2f GB; budget is not an RSS limit",
    base / 1024^3, per_worker / 1024^3))
  list(requested = requested, workers = as.integer(workers),
    base_bytes = base, worker_bytes = per_worker)
}
h1_prep_read_batch <- function(path) readRDS(path)
h1_prepare_coverage <- function(a, cfg, s) {
  started <- proc.time()[["elapsed"]]
  require_columns(a, c("transcript_id", "gene_id", "promoter_te_label"))
  # split() in the original preparation excludes missing gene IDs.
  gene_order <- names(split(a, a$gene_id))
  ann <- copy(a[!is.na(gene_id)])
  assert_unique(ann, "transcript_id")
  map <- reg_read(nb_path(s, "transcript_feature_batches.tsv"))
  manifest <- reg_read(nb_path(s, "transcript_batches.tsv"))
  require_columns(map, c("feature_id", "batch"))
  require_columns(manifest, c("batch", "input_file"))
  if (anyNA(map$feature_id) || anyNA(map$batch) || anyNA(manifest$batch))
    stop("Missing H1 batch index identifiers")
  assert_unique(map, "feature_id")
  assert_unique(manifest, "batch")
  if (any(!map$batch %in% manifest$batch)) stop("H1 feature index refers to an absent batch")
  required <- map[feature_id %in% ann$transcript_id]
  tasks <- manifest[batch %in% required$batch]
  if (anyNA(tasks$input_file) || anyDuplicated(tasks$input_file) ||
      any(!file.exists(tasks$input_file)) || any(dir.exists(tasks$input_file)) ||
      any(file.access(tasks$input_file, 4L) != 0L))
    stop("H1 batches require unique readable input files")
  if ("object_bytes" %in% names(tasks))
    tasks[, object_bytes := suppressWarnings(as.numeric(object_bytes))]
  resources <- h1_prep_resources(tasks, ann, map)
  # Only compact metadata for genes spanning batches must leave a worker.
  required[, gene_id := ann$gene_id[match(feature_id, ann$transcript_id)]]
  spanning <- required[, .(n_batches = uniqueN(batch)), by = gene_id][n_batches > 1L]$gene_id
  worker <- function(i) tryCatch({
    data.table::setDTthreads(1L)
    t0 <- proc.time()[["elapsed"]]
    row <- tasks[i]
    message(sprintf("H1 batch %d/%d (%s): start", i, nrow(tasks), row$batch))
    inp <- h1_prep_read_batch(row$input_file)
    if (!is.list(inp) || !is.data.frame(inp$data)) stop("Invalid count-batch data")
    z <- as.data.table(inp$data)
    require_columns(z, "feature_id")
    indexed <- map[batch == row$batch]$feature_id
    if (anyNA(z$feature_id) || any(!z$feature_id %in% indexed) ||
        any(!indexed %in% z$feature_id)) stop("Count-batch features disagree with the H1 batch index")
    wanted <- required[batch == row$batch]
    z <- z[feature_id %in% wanted$feature_id]
    z[, .h1_gene := wanted$gene_id[match(feature_id, wanted$feature_id)]]
    setindexv(z, ".h1_gene")
    ba <- ann[transcript_id %in% wanted$feature_id]
    cov <- rbindlist(lapply(split(ba, ba$gene_id), function(ga) {
      gd <- z[.(ga$gene_id[1L]), on = ".h1_gene"]
      gd[, .h1_gene := NULL]
      h1_counts(gd, ga, cfg, s)$coverage
    }), fill = TRUE)
    meta <- z[.h1_gene %in% spanning &
      is.finite(maternal_count) & is.finite(paternal_count) &
      round(maternal_count) + round(paternal_count) > 0,
      .(gene_id = .h1_gene, sample_id, cluster_id, cross, independence_verified)]
    meta <- unique(meta)
    elapsed <- proc.time()[["elapsed"]] - t0
    message(sprintf("H1 batch %d/%d (%s): complete; transcripts=%d; elapsed=%.2fs",
      i, nrow(tasks), row$batch, nrow(cov), elapsed))
    list(ok = TRUE, coverage = cov, meta = meta)
  }, error = function(e) list(ok = FALSE, error = sprintf(
    "H1 batch %s failed: %s", tasks$batch[i], conditionMessage(e))))
  result <- if (resources$workers > 1L) parallel::mclapply(seq_len(nrow(tasks)), worker,
    mc.cores = resources$workers, mc.preschedule = FALSE, mc.set.seed = FALSE) else
    lapply(seq_len(nrow(tasks)), worker)
  ok <- vapply(result, function(x) is.list(x) && isTRUE(x$ok), TRUE)
  if (any(!ok)) {
    errors <- vapply(result[!ok], function(x) if (is.list(x) && !is.null(x$error))
      x$error else "H1 worker terminated without a valid result", "")
    stop(paste(errors, collapse = "\n"))
  }
  meta <- rbindlist(lapply(result, `[[`, "meta"))
  if (nrow(meta)) {
    meta <- unique(meta)
    if (any(meta[, .N, by = .(gene_id, sample_id)]$N > 1L) ||
        any(meta[, uniqueN(cross), by = .(gene_id, cluster_id)]$V1 > 1L))
      stop("Inconsistent biological clusters across H1 batches")
  }
  missing <- ann[!transcript_id %in% required$feature_id]
  # Empty counts have no per-gene metadata to reconcile; annotate them together.
  empty <- if (nrow(missing)) list(h1_counts(data.table(), missing, cfg, s)$coverage) else list()
  coverage <- rbindlist(c(lapply(result, `[[`, "coverage"), empty), fill = TRUE)
  if (nrow(coverage)) {
    coverage[, .h1_order := match(gene_id, gene_order)]
    setorderv(coverage, c(".h1_order", "feature_id"))
    coverage[, .h1_order := NULL]
  }
  message(sprintf("H1 coverage complete: batches read=%d; transcripts=%d; elapsed=%.2fs",
    nrow(tasks), nrow(coverage), proc.time()[["elapsed"]] - started))
  coverage
}
