nb_context_resource_key <- function(cfg, s, root) {
  parent <- nb_te_require(cfg, s, root)$key
  if (cfg$profile != "human_trio") return(nb_hash(list(parent, "nonhuman")))
  z <- reg_settings(cfg)
  qspec <- reg_context_manifest(cfg, z, "qtl_tables")
  beds <- reg_context_manifest(cfg, z, "beds")
  paths <- c(if (!is.null(qspec)) qspec$path, if (!is.null(beds)) beds$path,
    unlist(lapply(z$human_context[intersect(c("qtl_tables", "beds", "phenotype_map"),
      names(z$human_context))], reg_path, cfg = cfg)),
    z$contig_aliases, cfg$paths$genome_fai)
  nb_hash(list(parent, resources = nb_file_key(paths)))
}
nb_context <- function(cfg, s, root) {
  nb_check_prepare(cfg, s, root); nb_te_require(cfg, s, root); nb_require(s, "combine")
  key <- nb_context_resource_key(cfg, s, root)
  if (nb_current(s, "context", key)) return(invisible(TRUE))
  output <- nb_path(s, "context", "human_transcript_context.tsv")
  if (cfg$profile != "human_trio") {
    nb_write(data.table(cohort_id = cfg$cohort_id, status = "not_applicable_nonhuman"), output)
    nb_finish(s, "context", key, output)
    return(invisible(TRUE))
  }
  z <- reg_settings(cfg)
  qspec <- reg_context_manifest(cfg, z, "qtl_tables")
  beds <- reg_context_manifest(cfg, z, "beds")
  # Same gate as the legacy human-context worker: peaks cannot stand in for cis-QTL.
  required <- reg_default(z$human_context$require_qtl, TRUE)
  if (isTRUE(required) && (is.null(qspec) || !nrow(qspec)))
    stop("context blocked: regulatory.human_context.require_qtl is TRUE but no qtl_tables ",
      "are configured; supply qtl_tables and a verified phenotype_map")
  # .fai has no header.
  target <- fread(cfg$paths$genome_fai, header = FALSE, select = 1)[[1]]
  aliases <- if (!is.null(z$contig_aliases)) reg_read(z$contig_aliases) else NULL
  repeat_path <- nb_repeat_path(cfg)
  repeats <- if (file.exists(repeat_path)) reg_read(repeat_path) else NULL
  tx <- reg_read(nb_path(s, "transcript_TE_annotation.tsv"))
  results <- reg_read(nb_path(s, "transcript_POE.tsv"))
  map <- reg_read(cfg$paths$event_map)[window == 50, .(transcript_id, event_id, gene_id)]
  assert_unique(map, "transcript_id", "Transcript-event context map")
  mapped_tx <- merge(tx, map[, .(transcript_id, event_id)], by = "transcript_id")
  events <- unique(map[, .(gene_id, event_id)])
  links <- data.table(event_id = character(), phenotype_id = character(),
    modality = character(), link_level = character())
  coverage <- data.table(status = "QTL_not_configured")
  if (!is.null(qspec) && nrow(qspec)) {
    if (is.null(z$human_context$phenotype_map)) stop("QTL context requires verified phenotype_map")
    q <- rbindlist(lapply(seq_len(nrow(qspec)), function(i) reg_context_qtl(qspec[i])), fill = TRUE)
    mapped <- reg_context_map(q, reg_read(reg_path(z$human_context$phenotype_map, cfg)), events)
    links <- mapped$links; coverage <- mapped$coverage
    # Event maps can group several transcripts. Do not call that exact transcript validation.
    links <- merge(links, unique(map[, .(event_id, transcript_id)]), by = "event_id", allow.cartesian = TRUE)
    links[, exact_transcript_phenotype := sub("\\.[0-9]+$", "", phenotype_id) ==
      sub("\\.[0-9]+$", "", transcript_id)]
    links[, POE_context_scope := fifelse(exact_transcript_phenotype,
      "matching_transcript_ID", "event_or_gene_mapping_not_transcript_validation")]
    # Hypothesis #5 TE term: is a lead/PIP variant inside the TE insertion that
    # overlaps THIS transcript's TSS? Reuse the legacy event-level routine with the
    # transcript as the unit (event_id temporarily carries transcript_id).
    if (!is.null(repeats) && nrow(links)) {
      setnames(links, "event_id", "context_event_id")
      links[, event_id := transcript_id]
      links <- reg_context_variant_links(links,
        mapped_tx[, .(chrom, tss, event_id = transcript_id)], repeats, target, aliases)
      links[, event_id := context_event_id][, context_event_id := NULL]
      if ("variant_at_event_TSS" %in% names(links))
        setnames(links, "variant_at_event_TSS", "variant_at_transcript_TSS")
    }
  }
  for (v in c("variant_at_transcript_TSS", "variant_overlaps_TSS_TE"))
    if (!v %in% names(links)) links[, (v) := as.logical(NA)]
  if (!"pip_all" %in% names(links)) links[, pip_all := NA_real_]
  if (!"variant_id" %in% names(links)) links[, variant_id := NA_character_]
  if (!"input_id" %in% names(links)) links[, input_id := NA_character_]
  peak <- reg_context_overlap(mapped_tx, beds, target, aliases)
  qtl_assessed <- !is.null(qspec) && nrow(qspec) > 0
  summary <- merge(tx, results, by.x = "transcript_id", by.y = "feature_id", all.x = TRUE)
  summary[, `:=`(n_qtl_links = 0L, n_peak_overlaps = 0L,
    n_event_qtl_links = if (qtl_assessed) 0L else NA_integer_,
    n_gene_context_links = if (qtl_assessed) 0L else NA_integer_,
    max_event_link_pip = NA_real_,
    # 0 is meaningful only where the TSS was assessed against repeats and QTL were supplied.
    n_qtl_variants_in_TSS_TE = if (qtl_assessed && !is.null(repeats)) 0L else NA_integer_,
    qtl_assessed = qtl_assessed,
    chromatin_assessed = !is.null(beds) && nrow(beds) > 0,
    TE_insertion_assessed = !is.null(repeats),
    interpretation = "baseline_transcript_POE_with_context_not_causal_or_colocalization")]
  summary[!(tss_te_label %in% c("TE", "nonTE")), n_qtl_variants_in_TSS_TE := NA_integer_]
  if (nrow(links)) {
    counts <- links[, .(N = .N,
      n_event_qtl_links = sum(link_level == "verified_phenotype_event"),
      n_gene_context_links = sum(link_level != "verified_phenotype_event"),
      max_event_link_pip = if (any(link_level == "verified_phenotype_event" & is.finite(pip_all)))
        max(pip_all[link_level == "verified_phenotype_event" & is.finite(pip_all)]) else NA_real_,
      n_qtl_variants_in_TSS_TE = uniqueN(paste(input_id, variant_id)[variant_overlaps_TSS_TE %in% TRUE])),
      by = transcript_id]
    summary[counts, on = "transcript_id", `:=`(n_qtl_links = i.N,
      n_event_qtl_links = i.n_event_qtl_links, n_gene_context_links = i.n_gene_context_links,
      max_event_link_pip = i.max_event_link_pip)]
    summary[counts[transcript_id %in% summary[!is.na(n_qtl_variants_in_TSS_TE)]$transcript_id],
      on = "transcript_id", n_qtl_variants_in_TSS_TE := i.n_qtl_variants_in_TSS_TE]
  }
  if (nrow(peak$hits)) {
    counts <- peak$hits[, .N, by = transcript_id]
    summary[counts, on = "transcript_id", n_peak_overlaps := i.N]
  }
  provenance <- data.table(version = NB_VERSION, build = z$build,
    unit = "transcript", qtl_required = isTRUE(required), qtl_assessed = qtl_assessed,
    chromatin_assessed = !is.null(beds) && nrow(beds) > 0,
    repeats = if (is.null(repeats)) "not_available" else basename(repeat_path),
    overlap_unit = "exact_transcript_TSS_or_same_TE_insertion",
    candidate_filter = "none_no_PIP_or_POE_significance_cutoff",
    pip_note = "max_event_link_pip is NA for lead-only QTL inputs (no PIP column)",
    caveat = "shared_GUSTO_or_annotation_selection_not_independent_replication")
  files <- nb_path(s, "context", c("human_transcript_context.tsv", "human_qtl_links.tsv.gz",
    "human_peak_overlaps.tsv.gz", "human_mapping_coverage.tsv", "human_resource_coverage.tsv",
    "human_context_provenance.tsv"))
  for (i in seq_along(files))
    nb_write(list(summary, links, peak$hits, coverage, peak$coverage, provenance)[[i]], files[i])
  nb_finish(s, "context", key, files)
}
nb_pilot <- function(cfg, s, root, n = 10L, te = FALSE) {
  nb_check_prepare(cfg, s, root)
  records <- list()
  if (te) {
    nb_te_require(cfg, s, root)
    jobs <- reg_read(nb_path(s, "within_gene_jobs.tsv"))
    if (!nrow(jobs)) stop("No eligible mixed gene for TE pilot")
    tx <- nb_estimable(reg_read(nb_path(s, "transcript_TE_POE.tsv"))[gene_id == jobs$gene_id[1]])
    d <- nb_gene_data(s, tx$transcript_id)
    started <- proc.time()[["elapsed"]]
    for (k in seq_len(n)) {
      set.seed(s$seed + k)
      sampled <- nb_bootstrap_units(d, cfg$profile == "public_cross")
      invisible(lapply(split(sampled, sampled$feature_id), nb_fit_one, cfg = cfg, s = s))
    }
    elapsed <- proc.time()[["elapsed"]] - started
    records[[1]] <- data.table(pilot = "TE_joint_refit", features = uniqueN(d$feature_id),
      draws = n, elapsed_seconds = elapsed,
      projected_seconds_this_gene = elapsed / n * s$bootstrap,
      caveat = "one_gene_pilot_not_a_guarantee_for_other_loci")
  } else for (lv in c("gene", "transcript")) {
    manifest <- reg_read(nb_path(s, paste0(lv, "_batches.tsv")))
    chosen <- manifest[n_features > 0][order(n_snp_rows)]
    if (!nrow(chosen)) next
    # Representative input batches across their SNP-row sizes.
    chosen <- chosen[unique(round(seq(1, .N, length.out = min(3L, .N))))]
    done <- 0L; elapsed <- 0; fit_status <- character(); mem <- 0
    for (file in chosen$input_file) {
      d <- readRDS(file)$data
      ids <- unique(d$feature_id)
      ids <- head(ids, max(1L, ceiling(n / nrow(chosen))))
      for (id in ids) {
        started <- proc.time()[["elapsed"]]
        result <- nb_fit_one(d[feature_id == id], cfg, s)
        elapsed <- elapsed + proc.time()[["elapsed"]] - started
        done <- done + 1L; fit_status <- c(fit_status, result$status)
      }
      mem <- max(mem, as.numeric(object.size(d)))
    }
    records[[lv]] <- data.table(pilot = lv, features = done, elapsed_seconds = elapsed,
      max_input_object_bytes = mem, fitted_ok = sum(fit_status %in% c("eligible_NB", "exploratory_NB")),
      projected_all_feature_seconds = elapsed / max(1, done) * sum(manifest$n_features),
      caveat = "input_heterogeneity_and_failed_fits_limit_projection")
  }
  result <- rbindlist(records, fill = TRUE)
  nb_write(result, nb_path(s, if (te) "TE_pilot_resources.tsv" else "pilot_resources.tsv"))
  print(result)
}
nb_status <- function(cfg, s, root) {
  stages <- c("prepare", "combine", "te-inputs", "te-combine", "between-genes", "context", "report-ready")
  verified <- function(expr) tryCatch({force(expr); TRUE}, error = function(e) FALSE)
  prepared <- verified(nb_check_prepare(cfg, s, root))
  te <- prepared && verified(nb_te_require(cfg, s, root))
  context <- te && verified(nb_require(s, "context", nb_context_resource_key(cfg, s, root)))
  rbindlist(lapply(stages, function(stage) {
    parent_ok <- if (stage %in% c("prepare", "combine")) prepared else
      if (stage %in% c("context", "report-ready")) context else te
    ok <- parent_ok && verified(nb_require(s, stage))
    data.table(cohort_id = cfg$cohort_id, stage = stage,
      status = if (ok) "complete_verified" else "missing_or_stale")
  }))
}
nb_report_ready <- function(cfg, s, root) {
  nb_check_prepare(cfg, s, root)
  nb_te_require(cfg, s, root)
  nb_require(s, "context", nb_context_resource_key(cfg, s, root))
  stages <- c("combine", "te-inputs", "te-combine", "between-genes", "context")
  parents <- lapply(stages, function(st) nb_require(s, st))
  nb_transcript_between(cfg, s)
  tx <- reg_read(nb_path(s, "transcript_TE_POE.tsv"))
  candidates <- nb_estimable(tx)[is.finite(p) & p < .05]
  candidates[, selection := "nominal_transcript_POE_not_TE_causality"]
  nb_write(candidates, nb_path(s, "nominal_transcript_candidates.tsv"))
  # No preselected old gene list, no restriction on shared markers or FDR.
  nb_write(data.table(cohort_id = cfg$cohort_id, species = cfg$species,
    version = NB_VERSION, n_gene_fits = nrow(reg_read(nb_path(s, "gene_POE.tsv"))),
    n_transcript_fits = nrow(reg_read(nb_path(s, "transcript_POE.tsv"))),
    n_nominal_transcript_candidates = nrow(candidates),
    marker_requirement = "all_informative_SNPs_per_feature_no_shared_isoform_marker",
    status = "ready_for_reanalysis"), nb_path(s, "report_summary.tsv"))
  nb_finish(s, "report-ready", nb_hash(lapply(parents, `[[`, "key")), nb_path(s,
    c("report_summary.tsv", "nominal_transcript_candidates.tsv",
      "transcript_between_gene_table.tsv", "transcript_between_gene_summary.tsv")))
}
nb_cross_summary <- function(configs, out, root) {
  if (file.exists(out)) stop("Refusing to overwrite cross-cohort summary")
  result <- rbindlist(lapply(configs, function(file) {
    cfg <- read_config_v2(file); s <- nb_settings(cfg)
    nb_check_prepare(cfg, s, root); nb_te_require(cfg, s, root)
    nb_require(s, "context", nb_context_resource_key(cfg, s, root))
    nb_require(s, "report-ready")
    d <- reg_read(nb_path(s, "gene_POE.tsv"))
    m <- load_manifests(cfg)
    data.table(cohort_id = cfg$cohort_id, species = cfg$species, reference_build = cfg$reference$build,
      tissue = paste(unique(m$units$tissue), collapse = ";"),
      stage = paste(unique(m$units$stage), collapse = ";"), n_features = nrow(d),
      n_estimable = nrow(nb_estimable(d)), n_nominal = sum(d$p < .05, na.rm = TRUE),
      n_cross_dependent = sum(d$not_strictly_POE %in% TRUE),
      median_signed = median(d$beta_adjusted, na.rm = TRUE),
      interpretation = "separate_cohort_tissue_stage_not_pooled_cross_species_inference")
  }), fill = TRUE)
  nb_write(result, out)
}
nb_export_files <- function(s) {
  fixed <- c("gene_POE.tsv", "transcript_POE.tsv", "gene_inventory.tsv.gz",
    "transcript_inventory.tsv.gz", "gene_availability.tsv", "transcript_availability.tsv",
    "input_audit_summary.tsv", "provenance.tsv", "package_versions.tsv", "baseline_settings.tsv",
    "transcript_TE_annotation.tsv", "transcript_TE_POE.tsv", "TE_availability.tsv",
    "within_gene_jobs.tsv", "within_gene_contrasts.tsv.gz", "within_gene_cohort_summary.tsv",
    "transcript_between_gene_table.tsv", "transcript_between_gene_summary.tsv",
    "nominal_transcript_candidates.tsv", "report_summary.tsv")
  paths <- nb_path(s, fixed)
  extras <- unlist(lapply(c("context", "between_genes", "within_gene"), function(dir)
    list.files(nb_path(s, dir), pattern = "\\.tsv(\\.gz)?$", full.names = TRUE, recursive = TRUE)))
  unique(c(paths, extras))
}
nb_export <- function(configs, out, root) {
  if (file.exists(out)) stop("Refusing to overwrite archive")
  parsed <- lapply(configs, read_config_v2)
  ids <- vapply(parsed, `[[`, "", "cohort_id")
  if (anyDuplicated(ids)) stop("Duplicate cohort IDs")
  stage <- tempfile("baseline-report-"); dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  top <- file.path(stage, sub("\\.tar\\.gz$", "", basename(out)))
  dir.create(top)
  manifest <- availability <- list()
  banned <- c("sample_id", "family_id", "child_id", "mother_id", "father_id", "unit_id",
    "cluster_id", "run_id", "infant_id", "sample_accession", "library_name")
  for (cfg in parsed) {
    s <- nb_settings(cfg)
    nb_check_prepare(cfg, s, root)
    nb_te_require(cfg, s, root)
    nb_require(s, "context", nb_context_resource_key(cfg, s, root))
    for (part in c("combine", "te-inputs", "te-combine", "between-genes", "context", "report-ready"))
      nb_require(s, part)
    files <- nb_export_files(s)
    for (src in files) {
      if (!file.exists(src) || file.info(src)$size <= 0) stop("Missing required baseline report file: ", src)
      header <- names(reg_read(src))
      if (any(header %in% banned)) stop("Privacy guard rejected identifier columns in ", basename(src))
      relative <- substring(src, nchar(s$output_dir) + 2L)
      rel <- file.path("cohorts", cfg$cohort_id, NB_VERSION, relative)
      dest <- file.path(top, rel); dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(src, dest)) stop("Copy failed: ", src)
      manifest[[length(manifest) + 1L]] <- data.table(cohort_id = cfg$cohort_id,
        part = NB_VERSION, archive_path = rel, bytes = file.info(dest)$size, sha256 = reg_sha(dest))
    }
    # Reuse privacy-safe cohort/settings/QC summaries, not legacy POE endpoints.
    cohort_dir <- file.path(top, "cohorts", cfg$cohort_id)
    if (exists("cohort_summary", mode = "function")) nb_write(cohort_summary(cfg), file.path(cohort_dir, "cohort.tsv"))
    if (exists("cohort_settings", mode = "function")) nb_write(cohort_settings(cfg), file.path(cohort_dir, "settings.tsv"))
    if (exists("qc_summary", mode = "function")) {
      qc <- qc_summary(cfg$paths$metrics, cfg$profile == "human_trio")
      if (!is.null(qc)) nb_write(qc, file.path(cohort_dir, "quantification_QC_summary.tsv"))
    }
    a <- reg_read(nb_path(s, "TE_availability.tsv"))
    availability[[cfg$cohort_id]] <- a
  }
  nb_write(rbindlist(manifest), file.path(top, "MANIFEST.tsv"))
  nb_write(rbindlist(availability), file.path(top, "AVAILABILITY.tsv"))
  nb_cross_summary(configs, file.path(top, "cross_cohort_summary.tsv"), root)
  dir.create(file.path(top, "docs"))
  for (f in c("docs/rerun_baseline_poe.md", "docs/model_comparison.md", "docs/BASELINE_NB_SCHEMA.md")) {
    if (!file.exists(file.path(root, f))) stop("Missing rerun documentation: ", f)
    file.copy(file.path(root, f), file.path(top, "docs", basename(f)))
  }
  writeLines(c("# Baseline NB results for reanalysis", "",
    "All-SNP gene and individual-transcript models replace the previous sentinel-marker analysis.",
    "Positive effects are paternal after parental dosage adjustment.",
    "Reciprocal interaction p < 0.05 is flagged not strictly POE.",
    "Small public cohorts remain exploratory. Missing repeats are unavailable, not negative results.",
    "No individual-level count tables, manifests, model objects or bootstrap records are included.",
    "MANIFEST.tsv records results checksums. AVAILABILITY.tsv records TE coverage.",
    "See docs for model, eligibility, limitations and rerun provenance."), file.path(top, "README.md"))
  old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(stage)
  utils::tar(out, files = basename(top), compression = "gzip", tar = "internal")
  listing <- utils::untar(out, list = TRUE)
  if (!length(listing) || !any(endsWith(listing, "/MANIFEST.tsv"))) stop("Archive validation failed")
  cat("Created ", basename(out), " (", file.info(out)$size, " bytes)\n", sep = "")
}
