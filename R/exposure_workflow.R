ex_batch_manifest <- function(s, level) {
  m <- ex_read(file.path(s$baseline_dir, paste0(level, "_batches.tsv")))
  a <- ex_read(file.path(s$baseline_dir, paste0(level, "_feature_batches.tsv")))
  ex_need(m, c("batch", "input_file", "object_bytes")); ex_unique(m, "batch")
  ex_need(a, c("feature_id", "batch")); ex_unique(a, "feature_id")
  if (any(!file.exists(m$input_file))) stop("Private baseline batches unavailable; run baseline prepare or restore their declared paths")
  m <- merge(a, m, by = "batch", all.x = TRUE)
  if (anyNA(m$input_file)) stop("Incomplete baseline batch linkage")
  m[, .(level = level, feature_id, source_batch = batch, input_file, object_bytes)]
}
ex_preflight <- function(cfg, s) {
  meta <- ex_validate_crosswalk(ex_read(s$sample_crosswalk), cfg)
  loci <- ex_read(s$te_loci); ex_validate_loci(loci, cfg)
  qc <- ex_read(s$quant_qc); counts <- ex_read(s$te_counts)
  ex_validate_quant(counts, qc, loci, meta, cfg, s)
  annotation <- ex_read(file.path(s$baseline_dir, "transcript_TE_annotation.tsv"))
  ex_unique(annotation, "transcript_id", "baseline transcript annotation")
  ex_need(annotation, c("promoter_te_label", "promoter_start", "promoter_end", "te_tss_window_bp"))
  radius <- cfg$regulatory$te_tss_window_bp %||ex% 3000L
  if (!radius %in% c(0, 1000, 3000) || anyNA(annotation$te_tss_window_bp) ||
      any(annotation$te_tss_window_bp != radius))
    stop("H6/H7 requires a freshly prepared annotation for the configured promoter radius")
  spans <- ex_gene_spans(annotation)
  if (!any(loci$chrom %in% spans$chrom)) stop("TE and transcript contigs do not overlap")
  targets <- ex_match_targets(annotation, loci)
  manifest <- data.table::rbindlist(lapply(c("gene", "transcript"), function(lv) ex_batch_manifest(s, lv)))
  h7_ready <- FALSE; h7_blocker <- "not_GUSTO"
  if (cfg$cohort_id == "GUSTO") {
    ans <- tryCatch({
      dictionary <- ex_read(s$exposure_dictionary); pf <- ex_read(s$pfas)
      ex_join_pfas(meta, pf, dictionary, cfg, s$common_pfas_samples)
    }, error = function(e) e)
    if (inherits(ans, "error")) h7_blocker <- conditionMessage(ans) else {
      meta <- ans; h7_ready <- TRUE; h7_blocker <- ""
    }
  }
  audit <- data.table::data.table(
    cohort_id = cfg$cohort_id, n_crosswalk_samples = nrow(meta), n_quantified_samples = nrow(qc),
    n_unquantified_crosswalk_samples = sum(!meta$sample_id %in% qc$sample_id),
    n_loci = nrow(loci), n_families = data.table::uniqueN(loci$family_id),
    n_transcript_family_pairs = sum(targets$level == "transcript"),
    n_gene_family_pairs = sum(targets$level == "gene"),
    private_count_input_bytes = sum(file.info(unique(manifest$input_file))$size),
    h7_ready = h7_ready, h7_blocker = h7_blocker,
    primary_TE_predictor_available = !s$family_only)
  ex_write(audit, ex_path(s, "preflight.tsv"))
  coverage <- data.table::rbindlist(lapply(if (h7_ready) c("matching_TE", EX_PFAS) else "matching_TE", function(p) {
    complete <- meta$sample_id %in% qc$sample_id
    if (p != "matching_TE") complete <- complete & is.finite(meta[[p]])
    data.table::data.table(exposure = p, n_crosswalk = nrow(meta),
      n_available = sum(complete), n_unavailable = sum(!complete))
  }))
  ex_write(coverage, ex_path(s, "exposure_coverage.tsv"))
  if (cfg$cohort_id == "GUSTO" && !h7_ready)
    message("H6 preparation can proceed. H7 blocked: ", h7_blocker)
  list(metadata = meta, loci = loci, qc = qc, counts = counts, annotation = annotation,
    spans = spans, targets = targets, manifest = manifest, h7_ready = h7_ready,
    h7_blocker = h7_blocker, audit = audit)
}
ex_prepare <- function(cfg, s) {
  obj <- ex_preflight(cfg, s)
  inputs <- c(s$sample_crosswalk, s$te_loci, s$te_counts, s$quant_qc,
    file.path(s$baseline_dir, c("gene_batches.tsv", "transcript_batches.tsv", "transcript_TE_annotation.tsv")),
    unique(obj$manifest$input_file),
    if (obj$h7_ready) c(s$pfas, s$exposure_dictionary),
    if (is.null(s$h5_annotations) && file.exists(file.path(s$baseline_dir, "context", "human_transcript_context.tsv")))
      file.path(s$baseline_dir, "context", "human_transcript_context.tsv"),
    unlist(s[intersect(c("annotation_manifest", "gene_sets", "h5_annotations"), names(s))]))
  if (any(!file.exists(inputs))) stop("Configured optional input is absent")
  ex_run_lock(s, cfg, inputs)
  predictors <- ex_predictors(obj$counts, obj$qc, obj$loci, obj$targets, obj$spans, s)
  selected_modes <- s$predictors %||ex% if (s$family_only) "inclusive_family_only" else "target_excluded"
  if (length(setdiff(selected_modes, unique(predictors$predictor))) && nrow(predictors))
    stop("Requested TE predictor unavailable: ", paste(setdiff(selected_modes, unique(predictors$predictor)), collapse = ", "))
  base <- merge(obj$targets, obj$manifest, by = c("level", "feature_id"), all = FALSE)
  h6 <- data.table::rbindlist(lapply(selected_modes, function(mode) {
    d <- data.table::copy(base); d[, `:=`(module = "h6", exposure = "matching_TE", predictor = mode)]; d
  }))
  jobs <- h6
  if (obj$h7_ready) {
    total <- data.table::rbindlist(lapply(EX_PFAS, function(p) {
      d <- data.table::copy(obj$manifest)
      tx <- obj$annotation$gene_id[match(d$feature_id, obj$annotation$transcript_id)]
      d[, gene_id := ifelse(level == "gene", feature_id, tx)]
      d[, `:=`(module = "h7_total", exposure = p, family_id = NA_character_, predictor = "none")]; d
    }))
    med <- data.table::rbindlist(lapply(EX_PFAS, function(p) {
      d <- data.table::copy(h6); d[, `:=`(module = "h7_mediation", exposure = p)]; d
    }))
    jobs <- data.table::rbindlist(list(jobs, total, med), fill = TRUE)
    if (isTRUE(s$modification)) {
      mod <- data.table::copy(med); mod[, module := "h7_modification"]
      jobs <- data.table::rbindlist(list(jobs, mod), fill = TRUE)
    }
  }
  if (nrow(jobs)) {
    data.table::setorder(jobs, input_file, module, feature_id, family_id, exposure, predictor)
    jobs[, job_id := sprintf("job_%08d", seq_len(.N))]
    jobs[, batch_local := ceiling(seq_len(.N) / s$batch_jobs), by = input_file]
    jobs[, batch := .GRP, by = .(input_file, batch_local)]
    jobs[, batch_local := NULL]
    ex_unique(jobs, "job_id", "planned model jobs")
  } else {
    jobs[, `:=`(job_id = character(), batch = integer())]
  }
  ex_write(jobs, ex_path(s, "private", "jobs.tsv"))
  ex_write(obj$targets, ex_path(s, "target_family_links.tsv"))
  insertion_links <- attr(obj$targets, "insertion_links")
  if (!is.null(insertion_links)) ex_write(insertion_links, ex_path(s, "target_insertion_links.tsv"))
  ex_write(predictors, ex_path(s, "private", "TE_predictors.tsv.gz"))
  scaling <- if (nrow(predictors)) unique(predictors[, .(gene_id, family_id, predictor,
    center, scale, n_expressed, predictor_status)]) else data.table::data.table()
  ex_write(scaling, ex_path(s, "TE_scaling.tsv"))
  # Do not serialize the full TE locus-count input into every model task.
  obj$counts <- NULL; obj$predictors <- predictors; obj$jobs <- jobs
  ex_save(obj, ex_path(s, "private", "prepared.rds"))
  tab <- if (nrow(jobs)) jobs[, .(jobs = .N, batches = data.table::uniqueN(batch)),
    by = .(module, level, predictor)] else data.table::data.table(module = "none", jobs = 0L)
  ex_write(tab, ex_path(s, "job_inventory.tsv"))
  if (!is.null(s$annotation_manifest) && !is.null(s$gene_sets)) {
    ann <- ex_read(s$annotation_manifest); memberships <- ex_read(s$gene_sets)
    ex_annotation_check(ann, memberships, cfg)
    ex_write(ann, ex_path(s, "frozen_annotation_manifest.tsv"))
    ex_save(list(manifest = ann, memberships = memberships,
      frozen_utc = format(Sys.time(), tz = "UTC")),
      ex_path(s, "private", "frozen_annotations.rds"))
  }
  ex_write(data.table::data.table(version = EX_VERSION, base_commit = EX_BASE_COMMIT,
    cohort_id = cfg$cohort_id, species = cfg$species, build = cfg$reference$build,
    no_remaining_confounding = TRUE, covariates = "sex,ancestry",
    temporal_ordering = "assumed", h7_ready = obj$h7_ready),
    ex_path(s, "provenance.tsv"))
  invisible(tab)
}
ex_load_prepared <- function(s) {
  p <- ex_path(s, "private", "prepared.rds")
  if (!file.exists(p)) stop("Run exposure prepare first")
  if ((!is.null(s$annotation_manifest) || !is.null(s$gene_sets)) &&
      !file.exists(ex_path(s, "private", "frozen_annotations.rds")))
    stop("Configured functional annotation was not frozen successfully; complete prepare first")
  obj <- readRDS(p)
  # Validate the prepared exposure run's descriptive metadata.
  lock <- readRDS(ex_path(s, "private", "run_identity.rds"))
  if (!identical(lock$identity$identity_schema, "readable_jobs_v2"))
    stop("Preparation schema unavailable: prepare in a fresh output_dir")
  for (f in lock$identity$files) {
    i <- file.info(f$path)
    if (!file.exists(f$path) || !identical(as.numeric(i$size), as.numeric(f$size)) ||
        !identical(as.numeric(i$mtime), as.numeric(f$mtime))) stop("Prepared input changed; use a new exposure output_dir")
  }
  if (!identical(lock$identity$settings, s)) stop("Settings changed; use a new exposure output_dir")
  obj
}
ex_job_metadata <- function(obj, job) {
  m <- data.table::copy(obj$metadata)
  if (job$module == "h6") {
    p <- obj$predictors[gene_id == job$gene_id & family_id == job$family_id & predictor == job$predictor]
    m[, E := p$TE_abundance[match(sample_id, p$sample_id)]]
  } else {
    if (!isTRUE(obj$h7_ready)) stop("H7 blocked: ", obj$h7_blocker)
    m[, E := get(job$exposure)]
    if (job$module != "h7_total") {
      p <- obj$predictors[gene_id == job$gene_id & family_id == job$family_id & predictor == job$predictor]
      m[, M := p$TE_abundance[match(sample_id, p$sample_id)]]
    }
  }
  m
}
ex_execute_job <- function(job, z, obj, cfg, s) {
  m <- ex_job_metadata(obj, job)
  ans <- tryCatch(switch(job$module,
    h6 = ex_nb(z, m, cfg, s, "h6"),
    h7_total = ex_nb(z, m, cfg, s, "total"),
    h7_mediation = ex_mediation(z, m, cfg, s),
    h7_modification = ex_mediation(z, m, cfg, s, modification = TRUE)),
    error = function(e) ex_failure("input_or_model_error", list(warnings = conditionMessage(e))))
  # Every planned job gets a row; exceptions remain visibly unavailable.
  row <- cbind(job[, .(job_id, batch, module, cohort_id = cfg$cohort_id, species = cfg$species,
    level, feature_id, gene_id, family_id, exposure, predictor)], ans$row)
  for (nm in intersect(c("insertion_ids", "n_insertions", "exact_TSS", "proximal_only",
      "min_distance_bp", "te_tss_window_bp"), names(job))) row[, (nm) := job[[nm]]]
  row
}
ex_fit_batch <- function(cfg, s, batch_id, mediation = FALSE) {
  if (mediation) ex_require_human(cfg)
  obj <- ex_load_prepared(s)
  if (mediation && !obj$h7_ready) stop("H7 blocked: ", obj$h7_blocker)
  jobs <- obj$jobs[obj$jobs$batch == as.integer(batch_id)]
  wanted <- if (mediation) c("h7_mediation", "h7_modification") else c("h6", "h7_total")
  jobs <- jobs[module %in% wanted]
  if (!nrow(jobs)) stop("No jobs of this stage for the requested batch")
  files <- unique(jobs$input_file); out <- list()
  start <- proc.time()[["elapsed"]]
  for (file in files) {
    inp <- readRDS(file)
    z <- if (is.list(inp) && "data" %in% names(inp)) inp$data else inp
    data.table::setDT(z)
    jj <- jobs[input_file == file]
    for (i in seq_len(nrow(jj))) out[[length(out) + 1L]] <- ex_execute_job(jj[i],
      z[feature_id == jj$feature_id[i]], obj, cfg, s)
  }
  results <- data.table::rbindlist(out, fill = TRUE)
  prefix <- if (mediation) "mediation" else "fit"
  dest <- ex_path(s, "private", "fits", sprintf("%s_%06d.tsv", prefix, batch_id))
  ex_write(results, dest)
  ex_write(data.table::data.table(batch = batch_id, stage = prefix, jobs = nrow(jobs),
    elapsed_seconds = proc.time()[["elapsed"]] - start, input_object_bytes = max(jobs$object_bytes)),
    sub(".tsv", "_timing.tsv", dest, fixed = TRUE))
  invisible(results)
}
ex_combine <- function(cfg, s) {
  obj <- ex_load_prepared(s); jobs <- obj$jobs
  if (!nrow(jobs)) {
    d <- jobs[, .(job_id, batch, module, level, feature_id, gene_id, family_id, exposure, predictor)]
    d[, `:=`(cohort_id = character(), status = character(), primary_p = numeric(), primary_q = numeric(),
      discovery = logical())]
    ex_write(d, ex_path(s, "model_results.tsv"))
    ex_write(data.table::data.table(module = "none", planned = 0L), ex_path(s, "model_coverage.tsv"))
    ex_write(data.table::data.table(selection_index = integer(), job_id = character()),
      ex_path(s, "bootstrap_selection.tsv"))
    return(invisible(d))
  }
  expected <- unique(jobs[, .(batch, prefix = ifelse(module %in% c("h7_mediation", "h7_modification"), "mediation", "fit"))])
  files <- ex_path(s, "private", "fits", sprintf("%s_%06d.tsv", expected$prefix, expected$batch))
  if (any(!file.exists(files))) stop("Missing fit stage files; run status for inventory")
  d <- data.table::rbindlist(lapply(files, ex_read), fill = TRUE)
  ex_unique(d, "job_id", "combined model jobs")
  if (!setequal(d$job_id, jobs$job_id)) stop("Missing/extra model job results")
  for (nm in c("slope_p", "indirect_path_p", "c_prime_p", "modification_p", "total_counts", "n_biological_units", "n_snps"))
    if (!nm %in% names(d)) d[, (nm) := NA_real_]
  d <- ex_finish(d)
  ex_write(d, ex_path(s, "model_results.tsv"))
  ex_write(d[, .N, by = .(module, level, predictor, status)], ex_path(s, "model_coverage.tsv"))
  selection <- d[module == "h7_mediation" & status == "eligible" & is.finite(indirect_q)]
  data.table::setorder(selection, indirect_q, indirect_path_p, level, feature_id, family_id, exposure, predictor)
  if (!isTRUE(s$bootstrap_all)) selection <- head(selection, s$bootstrap_top)
  selection[, selection_index := seq_len(.N)]
  ex_write(selection[, .(selection_index, job_id, level, feature_id, family_id, exposure, predictor)],
    ex_path(s, "bootstrap_selection.tsv"))
  invisible(d)
}
ex_bootstrap_job <- function(cfg, s, index) {
  ex_require_human(cfg)
  selection <- ex_read(ex_path(s, "bootstrap_selection.tsv"))
  chosen <- selection[selection_index == as.integer(index)]
  if (nrow(chosen) != 1L) stop("Unknown bootstrap selection index")
  obj <- ex_load_prepared(s); job <- obj$jobs[job_id == chosen$job_id]
  inp <- readRDS(job$input_file); z <- if (is.list(inp) && "data" %in% names(inp)) inp$data else inp
  data.table::setDT(z); z <- z[feature_id == job$feature_id]
  ans <- ex_bootstrap_pair(z, ex_job_metadata(obj, job), cfg, s)
  ex_write(cbind(chosen, ans$summary), ex_path(s, "bootstrap", paste0(job$job_id, ".tsv")))
  ex_write(ans$draws, ex_path(s, "private", "bootstrap_draws", paste0(job$job_id, ".tsv.gz")))
  invisible(ans$summary)
}
ex_pilot <- function(cfg, s, n = 9L) {
  obj <- ex_load_prepared(s); jobs <- obj$jobs
  if (!nrow(jobs)) stop("No model jobs available for pilot")
  data.table::setorder(jobs, object_bytes, level, module)
  idx <- unique(round(seq(1, nrow(jobs), length.out = min(n, nrow(jobs)))))
  out <- lapply(idx, function(i) {
    job <- jobs[i]; inp <- readRDS(job$input_file)
    z <- if (is.list(inp) && "data" %in% names(inp)) inp$data else inp
    data.table::setDT(z); z <- z[feature_id == job$feature_id]
    t <- proc.time()[["elapsed"]]
    ans <- ex_execute_job(job, z, obj, cfg, s)
    ans[, `:=`(elapsed_seconds = proc.time()[["elapsed"]] - t,
      input_object_bytes = as.numeric(object.size(inp)))]
    ans
  })
  d <- data.table::rbindlist(out, fill = TRUE)
  ex_write(d, ex_path(s, "pilot_results.tsv"))
  projection <- d[, .(n_pilot = .N, n_estimable = sum(status == "eligible"),
    median_fit_seconds = stats::median(elapsed_seconds), max_input_bytes = max(input_object_bytes)),
    by = .(module, level)]
  inventory <- jobs[, .(n_planned = .N), by = .(module, level)]
  projection <- merge(inventory, projection, all.x = TRUE)
  projection[, projected_serial_seconds := n_planned * median_fit_seconds]
  ex_write(projection, ex_path(s, "pilot_resources.tsv"))
  invisible(projection)
}
ex_run_enrichment <- function(cfg, s) {
  obj <- ex_load_prepared(s); results <- ex_read(ex_path(s, "model_results.tsv"))
  frozen <- ex_path(s, "private", "frozen_annotations.rds")
  if (!file.exists(frozen)) stop("Freeze annotation_manifest and gene_sets during prepare before fitting")
  ann <- readRDS(frozen)
  scores <- ex_gene_scores(results, obj$spans)
  out <- ex_enrich(scores, ann$memberships, ann$manifest, cfg, s)
  ex_write(scores, ex_path(s, "gene_evidence.tsv"))
  ex_write(out, ex_path(s, "functional_enrichment.tsv"))
  if (cfg$cohort_id == "GUSTO") {
    p <- s$h5_annotations %||ex% file.path(s$baseline_dir, "context", "human_transcript_context.tsv")
    if (file.exists(p)) {
      h5 <- ex_h5_adapter(ex_read(p)); h <- ex_h5_followup(results, h5, cfg)
      if ("p" %in% names(h)) h[, q := stats::p.adjust(p, "BH"), by = .(module, cohort_id, level, predictor)]
      ex_write(h, ex_path(s, "H5_support_enrichment.tsv"))
    } else stop("H5 context unavailable; functional results saved, H5 follow-up blocked")
    gp <- file.path(s$baseline_dir, "context", "gxe_ppBMI_context.tsv.gz")
    if (file.exists(gp)) ex_write(ex_h5_gxe_followup(results, ex_read(gp), cfg),
      ex_path(s, "H5_ppBMI_support_enrichment.tsv"))
  }
  invisible(out)
}
ex_export <- function(s, destination = ex_path(s, "aggregate_export")) {
  allowed <- c("model_results.tsv", "model_coverage.tsv", "job_inventory.tsv", "pilot_resources.tsv",
    "gene_evidence.tsv", "functional_enrichment.tsv", "H5_support_enrichment.tsv", "H5_ppBMI_support_enrichment.tsv",
    "provenance.tsv", "frozen_annotation_manifest.tsv", "bootstrap_selection.tsv", "report_summary.tsv",
    "TE_scaling.tsv", "exposure_coverage.tsv", "target_family_links.tsv", "target_insertion_links.tsv",
    file.path("bootstrap", list.files(ex_path(s, "bootstrap"), pattern = "[.]tsv$")))
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  copied <- character()
  for (nm in allowed) {
    path <- ex_path(s, nm); if (!file.exists(path)) next
    d <- ex_read(path)
    # File allowlist plus removal of private identifiers, raw exposures and errors.
    forbidden <- grepl("sample_id|SubjectID|cluster_id|unit_id|input_file|warning|blocker|^E$|^M$|family_count", names(d))
    if (any(forbidden)) d[, (names(d)[forbidden]) := NULL]
    ex_write(d, file.path(destination, nm)); copied <- c(copied, nm)
  }
  if (!length(copied)) stop("No aggregate output available")
  inventory <- data.table::data.table(file = copied,
    bytes = file.info(file.path(destination, copied))$size)
  ex_write(inventory, file.path(destination, "MANIFEST.tsv")); invisible(inventory)
}
ex_status <- function(s) {
  p <- ex_path(s, "private", "jobs.tsv")
  if (!file.exists(p)) return(data.table::data.table(stage = "prepare", status = "missing"))
  jobs <- ex_read(p)
  if (!nrow(jobs)) return(data.table::data.table(stage = "none", planned_batches = 0L, completed_batches = 0L))
  wanted <- unique(jobs[, .(batch, stage = ifelse(module %in% c("h7_mediation", "h7_modification"), "mediation", "fit"))])
  wanted[, complete := file.exists(ex_path(s, "private", "fits", sprintf("%s_%06d.tsv", stage, batch)))]
  wanted[, .(planned_batches = .N, completed_batches = sum(complete)), by = stage]
}
ex_pool_h6 <- function(paths, output_file) {
  # Optional final cross-cohort correction; never combines biological effect estimates.
  d <- data.table::rbindlist(lapply(paths, ex_read), fill = TRUE)
  ex_need(d, c("module", "cohort_id", "level", "predictor", "exposure", "job_id", "slope_p"))
  d <- d[module == "h6"]
  ex_unique(d, c("cohort_id", "job_id"), "pooled H6 jobs")
  if (data.table::uniqueN(d$cohort_id) < 2L) stop("Pooled-cohort q-values require at least two cohorts")
  out <- ex_finish(d); ex_write(out, output_file); invisible(out)
}
