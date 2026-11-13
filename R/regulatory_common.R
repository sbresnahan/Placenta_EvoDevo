# Opt-in regulatory analyses. No primary/additional_v1 file is modified.
reg_default <- function(x, default) if (is.null(x)) default else x
reg_text <- function(x, label) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) ||
      grepl("[\r\n\t]", x)) stop("Invalid ", label)
  x
}
reg_path <- function(x, cfg) {
  resolve_value(x, dirname(attr(cfg, "file")), unlist(cfg$environment_allowlist))
}
reg_settings <- function(cfg) {
  z <- reg_default(cfg$regulatory, list())
  base <- dirname(cfg$paths$model_gene)
  z$output_dir <- reg_path(reg_default(z$output_dir, file.path(base, "additional_v2")), cfg)
  z$cache_dir <- reg_path(reg_default(z$cache_dir,
    file.path(base, "reference", "regulatory_resources")), cfg)
  z$work_dir <- reg_path(reg_default(z$work_dir,
    file.path(cfg$paths$scratch, "regulatory_v2")), cfg)
  z$bootstrap <- as.integer(reg_default(z$bootstrap, 2000L))
  z$seed <- as.integer(reg_default(z$seed, 1147L))
  z$min_depth <- as.numeric(reg_default(z$min_depth, 10))
  z$max_pairs_per_gene <- as.integer(reg_default(z$max_pairs_per_gene, 10L))
  z$version <- "regulatory_v2.1"
  z$build <- reg_text(cfg$reference$build, "reference.build")
  for (k in c("bootstrap", "seed", "min_depth", "max_pairs_per_gene"))
    if (length(z[[k]]) != 1L || !is.finite(z[[k]]) || z[[k]] < 1)
      stop("Invalid regulatory.", k)
  if (z$bootstrap < 200L) stop("regulatory.bootstrap must be >=200")
  if (is.null(z$human_context)) z$human_context <- list()
  if (is.null(z$fish)) z$fish <- list()
  z$fish$image <- reg_path(reg_default(z$fish$image,
    "/risapps/apptainer/repo/tetools/20240523/tetools-202405023.sif"), cfg)
  z$fish$runtime <- reg_default(z$fish$runtime, "container")
  if (z$fish$runtime != "container" || !is.null(z$fish$tools) ||
      !is.null(z$fish$famdb_dir))
    stop("Only the container route is supported; remove native fish tool settings")
  z$fish$executable <- reg_text(reg_default(z$fish$executable, "apptainer"), "fish.executable")
  if (!is.null(z$fish$cpus)) {
    positive(z$fish$cpus, "regulatory.fish.cpus")
    if (z$fish$cpus != floor(z$fish$cpus))
      stop("regulatory.fish.cpus must be a whole allocated slot count")
    z$fish$cpus <- as.integer(z$fish$cpus)
  }
  z$fish$seed <- as.integer(reg_default(z$fish$seed, 1147L))
  if ((!is.null(z$fish$cpus) && (length(z$fish$cpus) != 1L ||
      !is.finite(z$fish$cpus) || z$fish$cpus < 4L)) ||
      !is.finite(z$fish$seed) || z$fish$seed < 1L)
    stop("Fish repeats require explicit >=4 allocated CPUs and a positive seed")
  if (!is.null(z$repeat_file)) z$repeat_file <- reg_path(z$repeat_file, cfg)
  if (!is.null(z$contig_aliases)) z$contig_aliases <- reg_path(z$contig_aliases, cfg)
  z$config_file <- attr(cfg, "file")
  z
}
reg_read <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) stop("Missing/empty input: ", path)
  if (grepl("\\.parquet$", path, ignore.case = TRUE)) {
    if (!requireNamespace("arrow", quietly = TRUE))
      stop("Parquet input requires R package arrow, or export the table as TSV")
    return(as.data.table(arrow::read_parquet(path)))
  }
  fread(file = path, na.strings = c("", "NA"), showProgress = FALSE)
}
reg_sha <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) stop("Missing fingerprint input: ", path)
  digest::digest(file = path, algo = "sha256")
}
reg_atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".reg-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Atomic receipt rename failed")
}
reg_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".reg-", tmpdir = dirname(path),
    fileext = if (endsWith(path, ".gz")) ".gz" else ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  fwrite(x, tmp, sep = "\t", na = "NA")
  if (!file.rename(tmp, path)) stop("Atomic table rename failed")
}
reg_receipt_path <- function(z, stage) file.path(z$output_dir, paste0(stage, ".rds"))
# Regulatory stages are immutable once completed. Receipts track completion and
# required output presence, not mutable input/config/source fingerprints.
# Older receipts with SHA256-named outputs remain readable.
reg_valid <- function(z, stage) {
  parents <- switch(stage, `fish-discovery` = "regulatory-resources",
    `fish-masking` = "fish-discovery", `te-tss` = "regulatory-resources",
    `human-context` = "te-tss", character())
  if (length(parents) && !all(vapply(parents, function(s) reg_valid(z, s), logical(1))))
    return(FALSE)
  if (stage == "te-tss" && file.exists(reg_receipt_path(z, "fish-masking")) &&
      !reg_valid(z, "fish-masking")) return(FALSE)
  p <- reg_receipt_path(z, stage)
  if (!file.exists(p)) return(FALSE)
  r <- tryCatch(readRDS(p), error = function(e) NULL)
  if (is.null(r) || !identical(r$state, "complete")) return(FALSE)
  outputs <- names(r$outputs)
  if (!length(outputs)) return(FALSE)
  info <- file.info(outputs)
  all(!is.na(info$size) & info$size > 0)
}
reg_begin <- function(z, stage, outputs, recover = FALSE) {
  r <- reg_receipt_path(z, stage)
  if (reg_valid(z, stage)) return(FALSE)
  if (!recover && (file.exists(r) || any(file.exists(outputs))))
    stop(stage, ": incomplete receipt or existing output; use a new output_dir",
      if (stage == "fish-discovery") " or explicit --resume RECOVERY_DIRECTORY" else "")
  reg_atomic_rds(list(state = "running", started = as.character(Sys.time())), r)
  TRUE
}
reg_finish <- function(z, stage, outputs) {
  p <- reg_receipt_path(z, stage)
  if (!file.exists(p)) stop("Missing running receipt: ", stage)
  old <- readRDS(p)
  if (!identical(old$state, "running")) stop("Stage is not running: ", stage)
  info <- file.info(outputs)
  if (anyNA(info$size) || any(info$size <= 0))
    stop("Missing/empty required stage output: ", stage)
  # Named output paths are the sole status contract. Do not rehash large
  # datasets or invalidate unrelated completed stages after a code edit.
  reg_atomic_rds(list(state = "complete", outputs = stats::setNames(
    rep(TRUE, length(outputs)), outputs), finished = as.character(Sys.time()),
    version = z$version), p)
}
reg_status <- function(z) {
  rbindlist(lapply(c("regulatory-resources", "fish-discovery", "fish-masking", "te-tss", "human-context"), function(s) {
    p <- reg_receipt_path(z, s)
    r <- if (file.exists(p)) tryCatch(readRDS(p), error = function(e) NULL) else NULL
    data.table(stage = s, state = if (is.null(r)) "not_run" else
      if (reg_valid(z, s)) "complete_outputs_present" else "incomplete_or_outputs_missing",
      receipt = p)
  }))
}
reg_alias <- function(chrom, target, alias = NULL) {
  chrom <- as.character(chrom); target <- as.character(target)
  canon <- function(v) sub("^chr", "", tolower(v))
  if (anyDuplicated(canon(target))) stop("Ambiguous target chromosome aliases")
  out <- target[match(chrom, target)]
  if (!is.null(alias)) {
    require_columns(alias, c("source", "target"), "Contig aliases")
    assert_unique(alias, "source", "Contig aliases")
    a <- alias$target[match(chrom, alias$source)]
    out[is.na(out)] <- a[is.na(out)]
  }
  miss <- which(is.na(out))
  out[miss] <- target[match(canon(chrom[miss]), canon(target))]
  if (any(!is.na(out) & !out %in% target)) stop("Alias maps outside target reference")
  mapping <- unique(data.table(source = chrom, target = out))[!is.na(target)]
  if (anyDuplicated(mapping$target)) stop("Multiple source contigs collapse onto one target")
  out
}
reg_intervals <- function(x, label = "intervals") {
  require_columns(x, c("chrom", "start", "end"), label)
  if (!is.numeric(x$start) || !is.numeric(x$end) || anyNA(x[, .(chrom, start, end)]) ||
      any(!is.finite(x$start) | !is.finite(x$end) | x$start < 1 |
      x$end < x$start | x$start != floor(x$start) | x$end != floor(x$end)))
    stop(label, ": require nonmissing 1-based closed integer coordinates")
  x
}
reg_hits <- function(query, subject) {
  reg_intervals(query); reg_intervals(subject)
  empty <- data.table(query = integer(), subject = integer())
  if (!nrow(query) || !nrow(subject)) return(empty)
  q <- copy(query[, .(chrom, start, end)])[, query := .I]
  s <- copy(subject[, .(chrom, start, end)])[, subject := .I]
  setkey(s, chrom, start, end)
  h <- foverlaps(q, s, nomatch = 0L)
  unique(h[, .(query, subject)])
}
reg_run <- function(command, args, stdout = "", stderr = "") {
  rc <- system2(command, vapply(args, shQuote, character(1)),
    stdout = stdout, stderr = stderr)
  if (!identical(as.integer(rc), 0L)) stop("Command failed (", rc, "): ", command)
  invisible(TRUE)
}
