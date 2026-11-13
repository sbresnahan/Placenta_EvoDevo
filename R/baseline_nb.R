# All-SNP baseline NB models. No sentinel-marker result is an input.
NB_VERSION <- "baseline_nb_v1"
NB_REFERENCE <- "c285f131071becec7347c830d263e78a88ac45de"
nb_default <- function(x, y) if (is.null(x)) y else x
nb_settings <- function(cfg) {
  x <- nb_default(cfg$baseline_nb, list())
  base <- dirname(cfg$paths$model_gene)
  out <- list(output_dir = nb_default(x$output_dir, file.path(base, NB_VERSION)),
    batches_gene = as.integer(nb_default(x$batches_gene, 20)),
    batches_transcript = as.integer(nb_default(x$batches_transcript, 100)),
    min_snps = 3L, min_total_counts = 50,
    min_samples = if (cfg$profile == "human_trio") 15L else 4L,
    min_per_cross = 2L, bootstrap = as.integer(nb_default(x$bootstrap, 2000)),
    seed = as.integer(nb_default(x$seed, 1147)),
    bootstrap_chunk = as.integer(nb_default(x$bootstrap_chunk, 25)),
    metadata = nb_default(x$metadata, cfg$paths$covariates),
    sex_column = nb_default(x$sex_column, "sex"),
    ancestry_column = nb_default(x$ancestry_column, "ancestry"),
    max_dense_gb = as.numeric(nb_default(x$max_dense_gb, 4)))
  for (k in c("output_dir", "metadata")) out[[k]] <- reg_path(out[[k]], cfg)
  for (k in c("batches_gene", "batches_transcript", "bootstrap", "bootstrap_chunk"))
    if (length(out[[k]]) != 1L || !is.finite(out[[k]]) || out[[k]] < 1)
      stop("Invalid baseline_nb.", k)
  if (out$bootstrap < 200) stop("Production bootstrap must be >=200; default 2000")
  out
}
nb_path <- function(s, ...) file.path(s$output_dir, ...)
nb_hash <- function(x) digest::digest(x, algo = "sha256")
nb_code_hash <- function(root) {
  files <- c(file.path("R", sort(list.files(file.path(root, "R"), pattern = "\\.R$"))),
    "poe_baseline.R", "bin/export_results.R")
  nb_hash(setNames(vapply(file.path(root, files), reg_sha, ""), files))
}
nb_runtime_key <- function() {
  packages <- c("glmmTMB", "TMB", "Matrix", "lme4", "data.table", "digest",
    "GenomicRanges", "Rsamtools", "Biostrings")
  nb_hash(list(R = as.character(getRversion()), packages = setNames(vapply(packages,
    function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not_installed"),
    ""), packages), contrasts = getOption("contrasts")))
}
nb_save <- function(x, path) reg_atomic_rds(x, path)
nb_write <- function(x, path) reg_write(as.data.table(x), path)
nb_file_key <- function(paths) {
  paths <- unique(paths)
  setNames(vapply(paths, reg_sha, ""), paths)
}
nb_receipt <- function(s, stage) nb_path(s, "receipts", paste0(stage, ".rds"))
nb_current <- function(s, stage, key) {
  p <- nb_receipt(s, stage)
  if (!file.exists(p)) return(FALSE)
  r <- readRDS(p)
  if (!identical(r$key, key))
    stop("Stale ", stage, ": inputs/code/settings changed. Use a NEW baseline_nb.output_dir.")
  if (!all(file.exists(names(r$outputs)))) return(FALSE)
  identical(nb_file_key(names(r$outputs)), r$outputs)
}
nb_finish <- function(s, stage, key, outputs) {
  nb_save(list(version = NB_VERSION, key = key, outputs = nb_file_key(outputs),
    finished = as.character(Sys.time())), nb_receipt(s, stage))
}
nb_require <- function(s, stage, key = NULL) {
  p <- nb_receipt(s, stage)
  if (!file.exists(p)) stop("Run baseline stage first: ", stage)
  r <- readRDS(p)
  if (!identical(r$version, NB_VERSION) ||
      (!is.null(key) && !identical(r$key, key)) ||
      !all(file.exists(names(r$outputs))) ||
      !identical(nb_file_key(names(r$outputs)), r$outputs))
    stop("Stale/incomplete baseline stage: ", stage, ". Never fall back to sentinel outputs.")
  r
}
nb_context_key <- function(cfg, s, root) {
  paths <- unlist(cfg$paths[intersect(c("libraries", "units", "tx_map",
    "event_map", "exons", "events",
    if (cfg$profile == "public_cross") "markers" else "trios"), names(cfg$paths))])
  if (cfg$profile == "human_trio" && isTRUE(cfg$gusto_precomputed$enabled))
    paths <- c(paths, unlist(cfg$gusto_precomputed[c("fastq_id_map", "origin_matrix")]))
  nb_hash(list(version = NB_VERSION, code = nb_code_hash(root), runtime = nb_runtime_key(), settings = s,
    config = nb_file_key(attr(cfg, "file")),
    manifests_and_parent_resources = nb_file_key(paths),
    metadata = if (file.exists(s$metadata)) reg_sha(s$metadata) else "not_present"))
}
nb_validate_raw <- function(d) {
  require_columns(d, c("transcript_id", "gene_id", "locus_id", "maternal_count",
    "paternal_count", "sample_id"))
  assert_unique(d, c("sample_id", "transcript_id", "locus_id"), "Transcript SNP count keys")
  if (any(!grepl("^[^:]+:[1-9][0-9]*:[ACGT]:[ACGT]$", d$locus_id)))
    stop("Transcript rows require full genomic SNP keys CHROM:POS:REF:ALT")
  v <- unlist(d[, .(maternal_count, paternal_count)])
  if (any(!is.na(v) & (!is.finite(v) | v < 0))) stop("Invalid expected counts")
  invisible(TRUE)
}
nb_read_tpms <- function(path, tx) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- sub(",\\s*$", "", lines[nzchar(lines) & !grepl("^[{}]$", lines)])
  at <- regexpr(":[^:]*$", lines)
  keep <- at > 0
  ids <- gsub('^"|"$', "", trimws(substr(lines[keep], 1, at[keep] - 1L)))
  val <- suppressWarnings(as.numeric(trimws(substring(lines[keep], at[keep] + 1L))))
  z <- data.table(transcript_id = ids, tpm = val)
  if (!nrow(z) || anyNA(z) || any(!is.finite(z$tpm) | z$tpm < 0))
    stop("Missing/invalid TPMs: ", path, ". Do not replace an unreadable TPM file with equal weights.")
  assert_unique(z, "transcript_id", "TPM IDs")
  # Exact IDs first. Only trailing numeric versions may be reconciled uniquely.
  missing <- setdiff(tx$transcript_id, z$transcript_id)
  if (length(missing)) {
    norm <- function(v) sub("\\.[0-9]+$", "", v)
    candidate <- z[!transcript_id %in% tx$transcript_id]
    candidate[, normalized := norm(transcript_id)]
    if (anyDuplicated(candidate$normalized)) stop("Ambiguous normalized TPM identifiers")
    add <- candidate[match(norm(missing), normalized)]
    add[, transcript_id := missing]
    z <- rbind(z, add[, .(transcript_id, tpm)])
  }
  w <- merge(tx[, .(transcript_id, gene_id)], z, by = "transcript_id", all.x = TRUE)
  if (anyNA(w$tpm)) stop("Reference transcripts missing from TPM file: ", path)
  w[, weight := if (sum(tpm) > 0) tpm / sum(tpm) else rep(1 / .N, .N), by = gene_id]
  w
}
nb_aggregate_gene <- function(d, weights) {
  nb_validate_raw(d)
  assert_unique(weights, "transcript_id", "Transcript usage weights")
  x <- merge(d, weights[, .(transcript_id, weight)], by = "transcript_id", all.x = TRUE)
  if (anyNA(x$weight)) stop("Allelic transcript absent from TPM weights")
  # Missing parental assignment invalidates that gene/SNP/sample, not the gene's other SNPs.
  x[, .(maternal_count = if (anyNA(maternal_count) || anyNA(paternal_count)) NA_real_ else sum(weight * maternal_count),
    paternal_count = if (anyNA(maternal_count) || anyNA(paternal_count)) NA_real_ else sum(weight * paternal_count)),
    by = .(sample_id, gene_id, locus_id)]
}
nb_metadata <- function(cfg, s, m) {
  lib <- copy(m$libraries)
  require_columns(lib, c("sample_id", "unit_id"))
  assert_unique(lib, "sample_id")
  if (cfg$profile == "human_trio") {
    link <- human_linkage(cfg, m)
    meta <- merge(lib[, .(sample_id, unit_id)], link[, .(sample_id, family_id)], by = "sample_id")
    meta[, cluster_id := family_id]
    meta[, `:=`(cross = NA_character_, independence_verified = TRUE)]
  } else {
    require_columns(m$units, c("unit_id", "maternal_parent", "paternal_parent", "cluster_id"))
    meta <- merge(lib[, .(sample_id, unit_id)], m$units, by = "unit_id")
    if (anyDuplicated(meta$unit_id))
      stop("Multiple configured libraries per biological unit: merge their raw lanes upstream; do not pseudo-replicate.")
    meta[, cross := fifelse(maternal_parent == cfg$parents$focal, "focal_mother",
      fifelse(maternal_parent == cfg$parents$alternate, "alternate_mother", NA_character_))]
    if (anyNA(meta$cross)) stop("Unrecognized reciprocal parents")
    meta[, independence_verified := !is.na(cluster_id) & nzchar(cluster_id)]
    # IDs below are computational sample blocks, explicitly NOT verified biological clusters.
    meta[!independence_verified, cluster_id := paste0("unverified_unit:", unit_id)]
  }
  if (file.exists(s$metadata)) {
    cov <- reg_read(s$metadata)
    require_columns(cov, "sample_id", "Baseline metadata")
    assert_unique(cov, "sample_id", "Baseline metadata")
    for (nm in c("sex", "ancestry")) {
      source <- s[[paste0(nm, "_column")]]
      if (source %in% names(cov)) {
        meta[, (nm) := as.character(cov[[source]][match(sample_id, cov$sample_id)])]
      }
    }
  }
  if (cfg$profile == "human_trio") {
    require_columns(meta, c("sex", "ancestry"), "Literal GUSTO baseline metadata")
    if (anyNA(meta[, .(sex, ancestry)]) || any(!nzchar(meta$sex) | !nzchar(meta$ancestry)))
      stop("Missing GUSTO sex/ancestry: supply baseline_nb.metadata and exact sample_id mapping")
    if (uniqueN(meta$sex) < 2 || uniqueN(meta$ancestry) < 2)
      stop("Literal human baseline requires estimable sex and ancestry factors")
  }
  for (nm in setdiff(c("sex", "ancestry"), names(meta))) meta[, (nm) := NA_character_]
  meta
}
nb_preflight <- function(cfg, s, root) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("Install glmmTMB in configured R library first")
  m <- load_manifests(cfg)
  reference_gate(cfg, root)
  meta <- nb_metadata(cfg, s, m)
  rows <- lapply(seq_len(nrow(meta)), function(i) {
    sample <- meta$sample_id[i]
    l <- m$libraries[sample_id == sample]
    files <- vapply(c("allele_counts", "tpms", "quant_receipt"), function(k)
      cfg_path(cfg, k, sample = sample, mode = "complete", strand = l$strand_mode), "")
    if (cfg$profile == "human_trio") {
      origin_key <- if (isTRUE(cfg$gusto_precomputed$enabled)) "child_tx_map" else "child_query"
      files <- c(files, setNames(cfg_path(cfg, origin_key, sample = sample), origin_key))
    }
    data.table(sample_id = sample, input = names(files), path = files,
      bytes = file.info(files)$size)
  })
  audit <- rbindlist(rows)
  if (anyNA(audit$bytes) || any(audit$bytes <= 0))
    stop("Raw quantification products missing; see paths resolved by cohort config. Existing aggregate tables are not sufficient.")
  invisible(list(manifests = m, metadata = meta, files = audit))
}
nb_prepare <- function(cfg, s, root) {
  checked <- nb_preflight(cfg, s, root)
  m <- checked$manifests; meta <- checked$metadata
  key <- nb_context_key(cfg, s, root)
  products <- checked$files[input != "quant_receipt"]$path
  key <- nb_hash(list(key, raw_products = nb_file_key(products)))
  if (nb_current(s, "prepare", key)) return(invisible(TRUE))
  dir.create(nb_path(s, "private"), recursive = TRUE, showWarnings = FALSE, mode = "0700")
  txmap <- reg_read(cfg$paths$tx_map)
  assert_unique(txmap, "transcript_id")
  audit <- list(); txfiles <- genefiles <- character()
  for (sample in meta$sample_id) {
    l <- m$libraries[sample_id == sample]
    getpath <- function(k) cfg_path(cfg, k, sample = sample, mode = "complete", strand = l$strand_mode)
    qout <- vapply(c("allele_counts", "thetas", "tpms", "quant_log", "resource_log"), getpath, "")
    if (!stage_valid(getpath("quant_receipt"), qout)) stop("Incomplete quantification: ", sample)
    tp <- nb_path(s, "private", paste0(sample, ".transcript.rds"))
    gp <- nb_path(s, "private", paste0(sample, ".gene.rds"))
    samplekey <- nb_hash(list(code = nb_code_hash(root),
      raw = nb_file_key(c(getpath("allele_counts"), getpath("tpms"))),
      parent_resources = key,
      references = nb_file_key(unlist(cfg$paths[c("tx_map", "exons", "event_map", "events")])),
      config = reg_sha(attr(cfg, "file"))))
    recname <- paste0("assign_", sample)
    if (!nb_current(s, recname, samplekey)) {
      x <- assign_sample(cfg, m, sample, return_transcripts = TRUE)
      x <- x[!is.na(locus_id) & grepl("^[^:]+:[1-9][0-9]*:[ACGT]:[ACGT]$", locus_id)]
      x[, sample_id := sample]
      nb_validate_raw(x)
      weights <- nb_read_tpms(getpath("tpms"), txmap)
      gene <- nb_aggregate_gene(x, weights)
      tx <- x[, .(sample_id, transcript_id, gene_id, locus_id, maternal_count, paternal_count)]
      nb_save(tx, tp); nb_save(gene, gp)
      nb_finish(s, recname, samplekey, c(tp, gp))
    }
    tx <- readRDS(tp); gene <- readRDS(gp)
    vals <- c(tx$maternal_count, tx$paternal_count, gene$maternal_count, gene$paternal_count)
    audit[[sample]] <- data.table(sample_id = sample, transcript_rows = nrow(tx),
      gene_rows = nrow(gene), missing_parental_rows = sum(!is.finite(tx$maternal_count) |
        !is.finite(tx$paternal_count)), fractional_values = sum(vals != round(vals), na.rm = TRUE),
      absolute_rounding_change = sum(abs(vals - round(vals)), na.rm = TRUE))
    txfiles <- c(txfiles, tp); genefiles <- c(genefiles, gp)
  }
  # One global, imbalance-blind representative per family. No locus-specific selection.
  if (cfg$profile == "human_trio") {
    depth <- rbindlist(lapply(txfiles, function(p) readRDS(p)[, .(
      coverage = sum(is.finite(maternal_count) & is.finite(paternal_count) &
        maternal_count + paternal_count > 0),
      depth = sum(maternal_count + paternal_count, na.rm = TRUE)), by = sample_id]))
    meta <- merge(meta, depth, by = "sample_id")
    setorder(meta, family_id, -coverage, -depth, sample_id)
    meta[, selected := seq_len(.N) == 1L, by = family_id]
  } else meta[, selected := TRUE]
  nb_save(meta, nb_path(s, "private", "metadata.rds"))
  nb_write(rbindlist(audit), nb_path(s, "private", "input_audit.tsv"))
  nb_write(rbindlist(audit)[, lapply(.SD, sum), .SDcols = setdiff(names(audit[[1]]), "sample_id")],
    nb_path(s, "input_audit_summary.tsv"))
  for (level in c("gene", "transcript")) {
    paths <- if (level == "gene") genefiles else txfiles
    d <- rbindlist(lapply(paths, readRDS))
    entity <- if (level == "gene") "gene_id" else "transcript_id"
    setnames(d, entity, "feature_id")
    d <- d[sample_id %in% meta[selected == TRUE]$sample_id]
    nb_make_batches(d, meta[selected == TRUE], level, s, key)
    rm(d); gc(FALSE)
  }
  nb_write(data.table(setting = c("version", "reference", "effect", "counts", "formula_human",
      "formula_reciprocal", "source_code", "prepared_utc"),
    value = c(NB_VERSION, NB_REFERENCE, "ln_paternal_over_maternal",
      "fractional_expected_round_at_fit", "allele+sex+ancestry+(1|SNPpos)+(1|sample_id)",
      "allele*cross+available_covariates+(1|SNPpos)+(1|sample_id)",
      nb_code_hash(root), as.character(Sys.time()))), nb_path(s, "provenance.tsv"))
  nb_write(data.table(package = rownames(installed.packages()),
    version = installed.packages()[, "Version"]), nb_path(s, "package_versions.tsv"))
  safe_settings <- s[setdiff(names(s), c("output_dir", "metadata"))]
  nb_write(data.table(setting = names(safe_settings),
    value = vapply(safe_settings, as.character, "")), nb_path(s, "baseline_settings.tsv"))
  outputs <- c(nb_path(s, "private", "metadata.rds"), nb_path(s, "input_audit_summary.tsv"),
    nb_path(s, "provenance.tsv"), nb_path(s, "package_versions.tsv"), nb_path(s, "baseline_settings.tsv"),
    unlist(lapply(c("gene", "transcript"), function(lv) nb_path(s,
      paste0(lv, c("_batches.tsv", "_inventory.tsv.gz", "_availability.tsv", "_feature_batches.tsv"))))))
  nb_finish(s, "prepare", key, outputs)
}
# Reference implementation, kept algebraically identical on complete matrices.
nb_size_factors <- function(M, P) {
  lib <- M + P
  geo <- exp(rowMeans(log(lib + 1), na.rm = TRUE)) - 1
  valid <- is.finite(geo) & geo > 0
  if (sum(valid) < 100) {
    keep <- rowSums(lib, na.rm = TRUE) > 0
    geo <- rep(NA_real_, nrow(lib))
    geo[keep] <- exp(rowMeans(log(lib[keep, , drop = FALSE] + .5), na.rm = TRUE)) - .5
    valid <- is.finite(geo) & geo > 0
  }
  ratios <- lib[valid, , drop = FALSE] / geo[valid]
  sf <- apply(ratios, 2, function(v) median(v[is.finite(v) & v > 0], na.rm = TRUE))
  sf <- sf / exp(mean(log(sf[is.finite(sf) & sf > 0]), na.rm = TRUE))
  sf[!is.finite(sf) | sf <= 0] <- 1
  setNames(as.numeric(sf), colnames(M))
}
nb_make_batches <- function(d, meta, level, s, prepare_key) {
  require_columns(d, c("feature_id", "sample_id", "locus_id", "maternal_count", "paternal_count"))
  assert_unique(d, c("feature_id", "sample_id", "locus_id"))
  # Sparse absence = no contributed counts, explicit unknown orientation = NA.
  # Build size factors by chunks of feature-SNP rows, not a duplicated long allele table.
  keys <- unique(d[, .(feature_id, locus_id)])
  setorder(keys, feature_id, locus_id)
  keys[, row := .I]
  d <- merge(d, keys, by = c("feature_id", "locus_id"))
  samples <- meta$sample_id
  matrix_gb <- nrow(keys) * length(samples) * 16 / 1024^3
  message(level, " normalization matrices: ", round(matrix_gb, 3),
    " GiB, excluding tables and fitting overhead")
  if (matrix_gb > s$max_dense_gb) stop("Normalization matrix memory exceeds baseline_nb.max_dense_gb. ",
    "Allocate additional LSF memory, then explicitly increase that limit; do not use a login node.")
  M <- P <- matrix(0, nrow(keys), length(samples), dimnames = list(NULL, samples))
  ij <- cbind(d$row, match(d$sample_id, samples))
  M[ij] <- d$maternal_count; P[ij] <- d$paternal_count
  sf <- nb_size_factors(M, P)
  inventory <- d[, .(n_samples_observed = sum(is.finite(maternal_count) &
      is.finite(paternal_count) & maternal_count + paternal_count > 0),
    total_counts = sum(maternal_count + paternal_count, na.rm = TRUE),
    n_unknown_parental = sum(!is.finite(maternal_count) | !is.finite(paternal_count))),
    by = .(feature_id, locus_id)]
  nb_write(inventory, nb_path(s, paste0(level, "_inventory.tsv.gz")))
  feature_stats <- d[, .(n_samples = uniqueN(sample_id[is.finite(maternal_count) &
      is.finite(paternal_count) & maternal_count + paternal_count > 0]),
    n_snps = uniqueN(locus_id[is.finite(maternal_count) & is.finite(paternal_count) &
      maternal_count + paternal_count > 0]),
    total_counts = sum(maternal_count + paternal_count, na.rm = TRUE)), by = feature_id]
  feature_stats[, prefilter := fifelse(n_snps < s$min_snps, "insufficient_snps",
    fifelse(n_samples < s$min_samples, "insufficient_samples",
      fifelse(total_counts < s$min_total_counts, "insufficient_counts", "pass")))]
  nb_write(feature_stats, nb_path(s, paste0(level, "_availability.tsv")))
  nbat <- s[[paste0("batches_", level)]]
  features <- feature_stats$feature_id
  cost <- keys[, .N, by = feature_id][match(features, feature_id)]$N
  o <- order(cost, decreasing = TRUE)
  assignment <- integer(length(features))
  assignment[o] <- rep(c(seq_len(nbat), rev(seq_len(nbat))), length.out = length(features))
  nb_write(data.table(feature_id = features, batch = assignment),
    nb_path(s, paste0(level, "_feature_batches.tsv")))
  manifests <- list()
  for (b in seq_len(nbat)) {
    ff <- features[assignment == b]
    ii <- which(keys$feature_id %in% ff)
    # All SNPs for each feature enter the same model. No all-isoform marker intersection.
    z <- data.table(feature_id = rep(keys$feature_id[ii], times = length(samples)),
      SNPpos = rep(keys$locus_id[ii], times = length(samples)),
      sample_id = rep(samples, each = length(ii)))
    if (!length(ii)) z <- data.table(feature_id = character(), SNPpos = character(), sample_id = character())
    z[, `:=`(maternal_count = as.vector(M[ii, , drop = FALSE]),
      paternal_count = as.vector(P[ii, , drop = FALSE]))]
    z <- merge(z, meta, by = "sample_id", all.x = TRUE)
    z[, offset := log(sf[sample_id])]
    file <- nb_path(s, "private", "batches", sprintf("%s_%04d.rds", level, b))
    nb_save(list(data = z, key = prepare_key, level = level), file)
    manifests[[b]] <- data.table(batch = b, level = level, n_features = length(ff),
      n_snp_rows = length(ii), n_long_rows = 2 * nrow(z), input_file = file,
      sha256 = reg_sha(file), object_bytes = as.numeric(object.size(z)))
  }
  nb_write(rbindlist(manifests), nb_path(s, paste0(level, "_batches.tsv")))
  nb_save(sf, nb_path(s, "private", paste0(level, "_size_factors.rds")))
}
nb_long <- function(z) {
  z <- z[is.finite(maternal_count) & is.finite(paternal_count)]
  x <- rbindlist(list(z[, .(feature_id, SNPpos, sample_id, allele = "maternal",
    count = maternal_count, offset, sex, ancestry, cross, cluster_id, independence_verified)],
    z[, .(feature_id, SNPpos, sample_id, allele = "paternal",
      count = paternal_count, offset, sex, ancestry, cross, cluster_id, independence_verified)]))
  # Unknown orientation is missing, not zero evidence.
  x <- x[is.finite(count) & count >= 0 & is.finite(offset)]
  x[, `:=`(count = round(count), allele = factor(allele, c("maternal", "paternal")),
    SNPpos = factor(SNPpos), sample_id = factor(sample_id))]
  x
}
nb_linear <- function(beta, V, weights, null = 0) {
  w <- setNames(rep(0, length(beta)), names(beta))
  w[names(weights)] <- weights
  b <- sum(w * beta); v <- as.numeric(t(w) %*% V %*% w)
  se <- if (is.finite(v) && v >= 0) sqrt(v) else NA_real_
  c(beta = b, se = se, p = if (is.finite(se) && se > 0) 2 * pnorm(-abs((b - null) / se)) else NA_real_)
}
nb_cross_flag <- function(p) ifelse(is.finite(p), p < .05, NA)
nb_fit_one <- function(z, cfg, s, bootstrap = FALSE) {
  id <- unique(z$feature_id)
  stopifnot(length(id) == 1L)
  x <- nb_long(z)
  sample_depth <- x[, .(depth = sum(count)), by = sample_id]
  positive <- as.character(sample_depth[depth > 0]$sample_id)
  n_snp <- uniqueN(as.character(x[count > 0]$SNPpos))
  public <- cfg$profile == "public_cross"
  dosage_null <- log(cfg$dosage$paternal / cfg$dosage$maternal)
  out <- data.table(feature_id = id, version = NB_VERSION, status = "filtered",
    n_samples = length(positive), n_snps = n_snp, total_counts = sum(x$count),
    beta = NA_real_, beta_adjusted = NA_real_, se = NA_real_, p = NA_real_,
    p_raw = NA_real_, interaction_beta = NA_real_, interaction_se = NA_real_,
    interaction_p = NA_real_, cross_main_beta = NA_real_, cross_main_p = NA_real_,
    beta_cross0 = NA_real_, beta_cross1 = NA_real_, se_cross0 = NA_real_, se_cross1 = NA_real_,
    beta_reference_cross = NA_real_, not_strictly_POE = NA,
    independence_verified = all(z$independence_verified),
    small_N = public, dosage_null = dosage_null, optimizer_code = NA_integer_,
    pdHess = NA, boundary_warning = FALSE, warnings = "", formula = "",
    interpretation = if (public) "exploratory_small_N" else "reference_NB_no_explicit_cis_adjustment")
  out$n_biological_clusters <- uniqueN(z[sample_id %in% positive]$cluster_id)
  if (length(positive) < s$min_samples) {out$status <- "insufficient_samples"; return(out)}
  if (out$n_biological_clusters < s$min_samples) {out$status <- "insufficient_biological_clusters"; return(out)}
  if (n_snp < s$min_snps) {out$status <- "insufficient_snps"; return(out)}
  if (sum(x$count) < s$min_total_counts) {out$status <- "insufficient_counts"; return(out)}
  terms <- if (public) "allele * cross" else "allele"
  for (nm in c("sex", "ancestry")) {
    values <- unique(na.omit(x[[nm]]))
    if (!public && (length(values) < 2 || anyNA(x[[nm]]))) {
      out$status <- paste0("unidentifiable_or_missing_", nm); return(out)
    }
    if (length(values) > 1) {
      if (anyNA(x[[nm]])) {out$status <- paste0("missing_", nm); return(out)}
      x[, (nm) := factor(get(nm))]
      terms <- c(terms, nm)
    }
  }
  if (public) {
    cross_counts <- unique(z[sample_id %in% positive, .(cluster_id, cross)])[, .N, by = cross]
    if (nrow(cross_counts) != 2L || anyNA(cross_counts$cross) ||
        min(cross_counts$N) < s$min_per_cross) {
      out$status <- "insufficient_reciprocal_units"; return(out)
    }
    x[, cross := factor(cross, c("focal_mother", "alternate_mother"))]
  }
  formula <- as.formula(paste("count ~", paste(c(terms, "(1|SNPpos)", "(1|sample_id)"), collapse = " + ")))
  out$formula <- paste(deparse(formula), collapse = "")
  fixed <- as.formula(paste("~", paste(terms, collapse = " + ")))
  X <- model.matrix(fixed, x)
  if (qr(X)$rank < ncol(X)) {out$status <- "rank_deficient"; return(out)}
  warning_text <- character()
  fit <- tryCatch(withCallingHandlers(glmmTMB::glmmTMB(formula, family = glmmTMB::nbinom2(),
    offset = offset, data = x, control = glmmTMB::glmmTMBControl(optimizer = nlminb,
      parallel = 1L, optCtrl = list(rel.tol = 1e-5, eval.max = 300L, iter.max = 300L))),
    warning = function(w) {warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning")}),
    error = function(e) e)
  out$warnings <- paste(unique(warning_text), collapse = " | ")
  out$boundary_warning <- any(grepl("singular|boundary|non-positive|convergence", warning_text, ignore.case = TRUE))
  if (inherits(fit, "error")) {out$status <- "fit_error"; out$warnings <- conditionMessage(fit); return(out)}
  out$optimizer_code <- fit$fit$convergence
  out$pdHess <- isTRUE(fit$sdr$pdHess)
  if (!out$pdHess || fit$fit$convergence != 0) {out$status <- "failed_convergence_or_Hessian"; return(out)}
  co <- summary(fit)$coefficients$cond
  b <- glmmTMB::fixef(fit)$cond; V <- as.matrix(vcov(fit)$cond)
  if (!all(is.finite(b)) || !all(is.finite(V))) {out$status <- "nonfinite_uncertainty"; return(out)}
  allele <- "allelepaternal"
  if (!allele %in% names(b)) {out$status <- "missing_allele_coefficient"; return(out)}
  raw <- nb_linear(b, V, setNames(1, allele))
  out$beta_reference_cross <- raw[["beta"]]
  if (public) {
    ki <- "allelepaternal:crossalternate_mother"; kc <- "crossalternate_mother"
    if (!all(c(ki, kc) %in% names(b))) {out$status <- "missing_cross_coefficient"; return(out)}
    other <- nb_linear(b, V, setNames(c(1, 1), c(allele, ki)))
    overall <- nb_linear(b, V, setNames(c(1, .5), c(allele, ki)))
    out[, `:=`(interaction_beta = b[[ki]], interaction_se = co[ki, 2],
      interaction_p = co[ki, 4], cross_main_beta = b[[kc]], cross_main_p = co[kc, 4],
      beta_cross0 = raw[["beta"]] - dosage_null,
      beta_cross1 = other[["beta"]] - dosage_null, se_cross0 = raw[["se"]],
      se_cross1 = other[["se"]], not_strictly_POE = nb_cross_flag(co[ki, 4]))]
  } else {overall <- raw; out$not_strictly_POE <- FALSE}
  out[, `:=`(beta = overall[["beta"]], se = overall[["se"]], p_raw = overall[["p"]],
    beta_adjusted = overall[["beta"]] - dosage_null)]
  if (!is.finite(out$se) || out$se <= 0) {out$status <- "nonfinite_uncertainty"; return(out)}
  out[, `:=`(p = 2 * pnorm(-abs(beta_adjusted / se)),
    ci_low = beta_adjusted - qnorm(.975) * se, ci_high = beta_adjusted + qnorm(.975) * se,
    log2_ratio_p_over_m = beta_adjusted / log(2),
    direction = ifelse(beta_adjusted > 0, "paternal", "maternal"),
    status = if (public) "exploratory_NB" else "eligible_NB")]
  if (!out$independence_verified) out$interpretation <- "exploratory_independence_unverified"
  if (isTRUE(out$not_strictly_POE)) out$interpretation <- paste(out$interpretation, "not_strictly_POE", sep = ";")
  out
}
nb_fit_batch <- function(cfg, s, root, level, batch_id) {
  nb_check_prepare(cfg, s, root)
  manifest <- reg_read(nb_path(s, paste0(level, "_batches.tsv")))
  row <- manifest[batch == as.integer(batch_id)]
  if (nrow(row) != 1L) stop("Unknown baseline batch")
  if (!identical(reg_sha(row$input_file), row$sha256)) stop("Batch input changed")
  inp <- readRDS(row$input_file)
  key <- nb_hash(list(inp$key, row$sha256, nb_code_hash(root)))
  stage <- sprintf("fit_%s_%04d", level, batch_id)
  output <- nb_path(s, "fits", sprintf("%s_%04d.tsv", level, batch_id))
  if (nb_current(s, stage, key)) return(invisible(TRUE))
  started <- proc.time()[["elapsed"]]
  results <- lapply(split(inp$data, inp$data$feature_id), nb_fit_one, cfg = cfg, s = s)
  result <- if (length(results)) rbindlist(results, fill = TRUE) else
    data.table(feature_id = character(), status = character(), version = character(),
      p = numeric(), beta_adjusted = numeric(), se = numeric(), interaction_p = numeric())
  nb_write(result, output)
  nb_write(data.table(level = level, batch = batch_id, n_features = nrow(result),
    elapsed_seconds = proc.time()[["elapsed"]] - started, input_bytes = row$object_bytes),
    sub("\\.tsv$", "_timing.tsv", output))
  nb_finish(s, stage, key, c(output, sub("\\.tsv$", "_timing.tsv", output)))
}
nb_check_prepare <- function(cfg, s, root) {
  r <- nb_require(s, "prepare")
  provenance <- reg_read(nb_path(s, "provenance.tsv"))
  if (!identical(provenance[setting == "source_code"]$value, nb_code_hash(root)))
    stop("Source changed after prepare: use a new output_dir")
  manifest <- reg_read(nb_path(s, "gene_batches.tsv"))
  inp <- readRDS(manifest$input_file[1])
  if (!identical(inp$key, r$key)) stop("Mixed prepared inputs")
  # Current raw/config key must match the prepare fingerprint, not just file presence.
  checked <- nb_preflight(cfg, s, root)
  key <- nb_hash(list(nb_context_key(cfg, s, root),
    raw_products = nb_file_key(checked$files[input != "quant_receipt"]$path)))
  if (!identical(r$key, key)) stop("Prepared inputs stale: config, raw files or metadata changed")
  invisible(r)
}
nb_combine <- function(cfg, s, root) {
  nb_check_prepare(cfg, s, root)
  outputs <- character()
  for (level in c("gene", "transcript")) {
    manifest <- reg_read(nb_path(s, paste0(level, "_batches.tsv")))
    for (b in manifest$batch) nb_require(s, sprintf("fit_%s_%04d", level, b))
    files <- nb_path(s, "fits", sprintf("%s_%04d.tsv", level, manifest$batch))
    result <- rbindlist(lapply(files, reg_read), fill = TRUE)
    assert_unique(result, "feature_id", "Combined baseline results")
    availability <- reg_read(nb_path(s, paste0(level, "_availability.tsv")))
    if (!setequal(result$feature_id, availability$feature_id)) stop("Missing/extra features in batch results")
    result[, `:=`(q = p.adjust(p, "BH"), interaction_q = p.adjust(interaction_p, "BH"),
      cohort_id = cfg$cohort_id, species = cfg$species, level = level,
      reference_build = cfg$reference$build, nominal_p_lt_0_05 = !is.na(p) & p < .05)]
    out <- nb_path(s, paste0(level, "_POE.tsv"))
    nb_write(result, out); outputs <- c(outputs, out)
  }
  key <- nb_hash(list(prepare = nb_require(s, "prepare")$key, results = nb_file_key(outputs)))
  nb_finish(s, "combine", key, outputs)
}
