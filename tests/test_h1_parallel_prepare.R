#!/usr/bin/env Rscript
# Fixed software fixtures only; no biological inference or model fitting.
args <- commandArgs(TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages(library(data.table))
for (m in c("common", "config", "regulatory_common", "baseline_nb", "baseline_te", "h1_bayes"))
  source(file.path(root, "R", paste0(m, ".R")))

# Frozen preparation reference from the preceding delivered implementation.
h1_serial_reference <- function(cfg, s, root) {
  if (cfg$profile != "public_cross") stop("Bayesian H1 uses reciprocal small-cohort profiles")
  hs <- h1_settings(cfg); nb_te_require(cfg, s, root); nb_check_prepare(cfg, s, root)
  outputs <- nb_path(s, c("h1_transcript_coverage.tsv", "h1_gene_coverage.tsv", "h1_bayes_jobs.tsv"))
  if (nb_current(s, "h1-inputs", NULL)) return(invisible(TRUE))
  a <- reg_read(nb_path(s, "transcript_TE_annotation.tsv"))
  coverage <- rbindlist(lapply(split(a, a$gene_id), function(ann) {
    z <- nb_gene_data(s, ann$transcript_id)
    h1_counts(z, ann, cfg, s)$coverage
  }), fill = TRUE)
  genes <- coverage[, .(n_TE = sum(eligibility == "eligible" & promoter_te_label == "TE"),
    n_nonTE = sum(eligibility == "eligible" & promoter_te_label == "nonTE"),
    n_transcripts = .N, n_count_eligible = sum(eligibility == "eligible")), by = gene_id]
  genes[, eligibility := fifelse(n_TE > 0 & n_nonTE > 0, "eligible", "no_eligible_mixed_classes")]
  jobs <- copy(genes[eligibility == "eligible"]); setorder(jobs, gene_id); jobs[, job_index := .I]
  coverage[, te_tss_window_bp := reg_settings(cfg)$te_tss_window_bp]
  genes[, `:=`(count_provenance = hs$count_provenance,
    assignment_covariance = "shared_random_effects_leave_unknown_read_assignment_covariance_unestimated")]
  for (i in seq_along(outputs)) nb_write(list(coverage, genes, jobs)[[i]], outputs[i])
  nb_finish(s, "h1-inputs", list(settings = hs, profiles = h1_profiles(),
    parent = nb_require(s, "te-inputs")$key), outputs)
}

main <- function() {
  checks <- character()
  check <- function(label, x) {
    if (!isTRUE(x)) stop("FAIL ", label)
    checks <<- c(checks, label); cat("PASS ", label, "\n", sep = "")
  }
  fails <- function(x, pattern = NULL) {
    e <- tryCatch({force(x); NULL}, error = identity)
    inherits(e, "error") && (is.null(pattern) || grepl(pattern, conditionMessage(e)))
  }
  controls <- c("H1_PREP_WORKERS", "H1_PREP_MEMORY_GB", "LSB_DJOB_NUMPROC",
    "LSB_MAX_NUM_PROCESSORS")
  old <- Sys.getenv(controls, unset = NA_character_)
  on.exit(for (i in seq_along(controls)) if (is.na(old[i])) Sys.unsetenv(controls[i]) else
    do.call(Sys.setenv, setNames(list(old[i]), controls[i])), add = TRUE)
  Sys.unsetenv(controls)
  setDTthreads(1L)
  home <- tempfile("h1-parallel-fixture-"); dir.create(home)
  on.exit(unlink(home, recursive = TRUE), add = TRUE)
  events <- file.path(home, "read_events"); dir.create(events)
  # Separate event files also record reads performed by forked workers.
  assign("readRDS", function(file, ...) {
    if (grepl("counts_[0-9]+[.]rds$", file))
      writeLines(file, tempfile(paste0("read_", Sys.getpid(), "_"), tmpdir = events))
    base::readRDS(file, ...)
  }, envir = .GlobalEnv)
  on.exit(rm("readRDS", envir = .GlobalEnv), add = TRUE)
  cfg <- list(profile = "public_cross", reference = list(build = "software"),
    dosage = list(maternal = 2, paternal = 1),
    paths = list(model_gene = file.path(home, "unused.tsv"), scratch = home))
  attr(cfg, "file") <- file.path(home, "software.yaml")
  s0 <- list(min_snps = 3L, min_samples = 4L, min_total_counts = 50, min_per_cross = 2L)
  # Counts and reciprocal layout follow test_h1_bayes_contract.R.
  ann <- CJ(gene_id = sprintf("G%03d", 1:12), tx = 1:4)
  ann[, `:=`(transcript_id = paste0(gene_id, "_T", tx),
    promoter_te_label = c("TE", "nonTE", "TE", "unknown_repeat")[tx])]
  map <- ann[, .(feature_id = transcript_id, batch = tx)]
  z <- CJ(feature_id = ann$transcript_id, sample_id = paste0("S", 1:4), j = 1:3)
  z[, `:=`(SNPpos = paste0(feature_id, ":marker", j), maternal_count = 4.4,
    paternal_count = 2.2, cluster_id = sample_id, independence_verified = TRUE,
    cross = ifelse(sample_id %in% c("S1", "S2"), "focal_mother", "alternate_mother"))]
  z[feature_id == "G002_T1" & sample_id == "S1" & j == 1L, maternal_count := NA_real_]
  z[feature_id == "G002_T1" & sample_id == "S1" & j == 2L,
    `:=`(maternal_count = 0, paternal_count = 0)]
  z[startsWith(feature_id, "G003") & sample_id == "S1", independence_verified := FALSE]
  # Cross-gene sample differences are allowed; do not impose a cohort-wide rule.
  z[startsWith(feature_id, "G004"), cluster_id := paste0("other_", sample_id)]
  z[feature_id == "G005_T3", `:=`(maternal_count = 0, paternal_count = 0)]
  z[feature_id == "G006_T3", maternal_count := NA_real_]
  z <- z[!(feature_id == "G007_T3" & j == 3L)]
  z <- z[!(feature_id == "G008_T3" & sample_id == "S4")]
  z[feature_id == "G009_T3", `:=`(maternal_count = 1, paternal_count = 0)]
  z[startsWith(feature_id, "G010") & sample_id == "S2", cluster_id := "S1"]
  ann <- rbind(ann[, .(transcript_id, gene_id, promoter_te_label)],
    data.table(transcript_id = c("absent", "absent_same_gene"), gene_id = c("Z", "G001"),
      promoter_te_label = c("TE", "nonTE")))
  ann <- ann[rev(seq_len(.N))]
  paths <- file.path(home, paste0("counts_", 1:4, ".rds"))
  pieces <- lapply(1:4, function(b) z[feature_id %in% map[batch == b]$feature_id])
  for (b in 1:4) saveRDS(list(data = pieces[[b]]), paths[b])
  manifest <- data.table(batch = 1:4, input_file = paths,
    object_bytes = vapply(pieces, function(x) as.numeric(object.size(x)), 0))
  stage <- function(label, mp = map, mf = manifest, aa = ann) {
    s <- s0; s$output_dir <- file.path(home, label); dir.create(s$output_dir)
    nb_write(aa, nb_path(s, "transcript_TE_annotation.tsv"))
    nb_write(mp, nb_path(s, "transcript_feature_batches.tsv"))
    nb_write(mf, nb_path(s, "transcript_batches.tsv"))
    nb_finish(s, "prepare", list(software = TRUE), mf$input_file)
    nb_finish(s, "te-inputs", list(software = TRUE), nb_path(s, "transcript_TE_annotation.tsv"))
    s
  }
  tables <- c("h1_transcript_coverage.tsv", "h1_gene_coverage.tsv", "h1_bayes_jobs.tsv")
  raw <- function(p) readBin(p, "raw", n = file.info(p)$size)
  runs <- list()
  run <- function(label, fn, workers) {
    Sys.setenv(H1_PREP_WORKERS = workers, H1_PREP_MEMORY_GB = 4)
    ss <- stage(label)
    before <- list.files(events, full.names = TRUE)
    elapsed <- unname(system.time(fn(cfg, ss, root))[["elapsed"]])
    new <- setdiff(list.files(events, full.names = TRUE), before)
    reads <- vapply(new, function(p) readLines(p)[1L], "")
    runs[[label]] <<- data.table(implementation = label, workers = workers,
      count_file_reads = length(new), elapsed_seconds = elapsed,
      batch_object_bytes = sum(manifest$object_bytes), count_rows = nrow(z))
    list(s = ss, reads = reads)
  }
  ref <- run("old_serial", h1_serial_reference, 1L)
  one <- run("batch_serial", h1_prepare, 1L)
  two <- run("batch_parallel", h1_prepare, 2L)
  for (name in tables) {
    check(paste("serial output bytes match reference", name),
      identical(raw(nb_path(ref$s, name)), raw(nb_path(one$s, name))))
    check(paste("parallel output bytes match reference", name),
      identical(raw(nb_path(ref$s, name)), raw(nb_path(two$s, name))))
    check(paste("parsed serial/parallel values and types match", name),
      identical(as.data.frame(reg_read(nb_path(one$s, name))),
        as.data.frame(reg_read(nb_path(two$s, name)))))
  }
  ref_cov <- rbindlist(lapply(split(ann, ann$gene_id), function(ga)
    h1_counts(nb_gene_data(ref$s, ga$transcript_id), ga, cfg, ref$s)$coverage), fill = TRUE)
  check("in-memory coverage types and missing values match reference",
    identical(as.data.frame(ref_cov),
      as.data.frame(h1_prepare_coverage(ann, cfg, one$s))))
  check("each required count file read once in serial", length(one$reads) == 4L &&
    all(table(one$reads) == 1L))
  check("each required count file read once in parallel", length(two$reads) == 4L &&
    all(table(two$reads) == 1L))
  check("reference count reads repeat across genes", length(ref$reads) == 48L)
  check("fixture covers support and class exclusions", all(c("eligible", "unknown_TE_class",
    "no_informative_counts", "independence_unverified", "insufficient_snps",
    "insufficient_samples", "insufficient_counts", "insufficient_reciprocal_clusters") %in%
    ref_cov$eligibility))
  before <- length(list.files(events))
  h1_prepare(cfg, two$s, root)
  check("completed-stage reuse performs no count reads", before == length(list.files(events)))
  check("success receipt recorded only after all tables",
    nb_current(two$s, "h1-inputs", NULL))
  r <- h1_prep_resources(manifest[1:2], ann, map)
  check("workers capped by available batches", r$workers == 2L)
  Sys.setenv(H1_PREP_WORKERS = 8, LSB_DJOB_NUMPROC = 1)
  check("workers capped by scheduler CPUs", h1_prep_resources(manifest, ann, map)$workers == 1L)
  Sys.unsetenv("LSB_DJOB_NUMPROC")
  cap <- (r$base_bytes + 1.1 * r$worker_bytes) / 1024^3
  Sys.setenv(H1_PREP_MEMORY_GB = cap)
  check("memory estimate reduces concurrency", h1_prep_resources(manifest, ann, map)$workers == 1L)
  Sys.setenv(H1_PREP_MEMORY_GB = .1)
  check("insufficient memory budget fails before batch loading",
    fails(h1_prep_resources(manifest, ann, map), "one-worker"))
  Sys.setenv(H1_PREP_MEMORY_GB = 4)
  fallback <- character()
  withCallingHandlers(r <- h1_prep_resources(manifest[, .(batch, input_file)], ann, map),
    message = function(m) {fallback <<- c(fallback, conditionMessage(m)); invokeRestart("muffleMessage")})
  check("missing object sizes explicitly select serial loading", r$workers == 1L &&
    any(grepl("object_bytes unavailable", fallback, fixed = TRUE)))
  check("missing size metadata still preserves output",
    !fails(h1_prepare(cfg, stage("sizes_absent", mf = manifest[, .(batch, input_file)]), root)))
  Sys.unsetenv("H1_PREP_MEMORY_GB")
  check("parallel request without memory budget selects serial",
    h1_prep_resources(manifest, ann, map)$workers == 1L)
  Sys.setenv(H1_PREP_MEMORY_GB = 4, H1_PREP_WORKERS = "1.5")
  check("noninteger worker control rejected", fails(h1_prep_resources(manifest, ann, map)))
  Sys.setenv(H1_PREP_WORKERS = 2, H1_PREP_MEMORY_GB = "bad")
  check("invalid memory control rejected", fails(h1_prep_resources(manifest, ann, map)))
  Sys.setenv(H1_PREP_MEMORY_GB = 4)
  error_run <- function(label, mp = map, mf = manifest, pattern = NULL) {
    ss <- stage(label, mp, mf)
    check(paste(label, "fails"), fails(h1_prepare(cfg, ss, root), pattern))
    check(paste(label, "writes neither tables nor success receipt"),
      !any(file.exists(nb_path(ss, tables))) && !file.exists(nb_receipt(ss, "h1-inputs")))
  }
  error_run("duplicate_feature_index", rbind(map, map[1L]))
  error_run("missing_batch_index", transform(map, batch = 99L))
  error_run("duplicate_manifest", mf = rbind(manifest, manifest[1L]))
  bad_mf <- copy(manifest); bad_mf[2L, input_file := paths[1L]]
  error_run("duplicate_batch_file", mf = bad_mf)
  bad_mf <- copy(manifest); bad_mf[2L, input_file := file.path(home, "missing.rds")]
  error_run("unreadable_batch", mf = bad_mf)
  # Worker-side errors and cross-batch consistency on count-ineligible features.
  mutate_batch <- function(label, change, pattern, reference_error = TRUE) {
    saveRDS(list(data = change(copy(pieces[[4L]]))), paths[4L])
    if (reference_error) check(paste(label, "also rejected by original preparation"),
      fails(h1_serial_reference(cfg, stage(paste0(label, "_reference")), root)))
    error_run(label, pattern = pattern)
    saveRDS(list(data = pieces[[4L]]), paths[4L])
  }
  mutate_batch("negative_counts", function(x) {x[1L, maternal_count := -1]; x}, "Negative")
  mutate_batch("duplicate_count_row", function(x) rbind(x, x[1L]), "unique|Duplicate|duplicate")
  mutate_batch("wrong_batch_features", function(x) {x[1L, feature_id := "not_indexed"]; x},
    "features disagree", reference_error = FALSE)
  mutate_batch("cross_batch_sample_conflict", function(x) {
    x[feature_id == "G001_T4" & sample_id == "S1", cluster_id := "different"]; x
  }, "Inconsistent biological clusters")
  mutate_batch("cross_batch_cluster_conflict", function(x) {
    x[feature_id == "G001_T4", `:=`(sample_id = paste0(sample_id, "_different"),
      cross = ifelse(cross == "focal_mother", "alternate_mother", "focal_mother"))]; x
  }, "Inconsistent biological clusters")
  mutate_batch("cross_batch_independence_conflict", function(x) {
    x[feature_id == "G001_T4" & sample_id == "S1", independence_verified := FALSE]; x
  }, "Inconsistent biological clusters")
  # Metadata on excluded zero/missing parental rows must not create new conflicts.
  saveRDS(list(data = {
    zz <- copy(pieces[[4L]])
    zz[feature_id == "G001_T4", `:=`(maternal_count = 0, paternal_count = 0,
      cluster_id = "excluded_conflict")]; zz
  }), paths[4L])
  check("zero-count metadata remain excluded from consistency checks",
    !fails(h1_prepare(cfg, stage("excluded_zero_metadata"), root)))
  saveRDS(list(data = pieces[[4L]]), paths[4L])
  saveRDS(list(data = {
    zz <- copy(pieces[[4L]])
    zz[feature_id == "G001_T4", `:=`(maternal_count = NA_real_,
      cluster_id = "excluded_conflict")]; zz
  }), paths[4L])
  check("missing-parental metadata remain excluded from consistency checks",
    !fails(h1_prepare(cfg, stage("excluded_missing_metadata"), root)))
  saveRDS(list(data = pieces[[4L]]), paths[4L])
  writeLines("deliberately unreadable software RDS fixture", paths[4L])
  error_run("invalid_rds")
  saveRDS(list(data = pieces[[4L]]), paths[4L])
  before <- length(list.files(events))
  h1_prepare_coverage(ann[transcript_id == "G001_T1"], cfg, stage("one_required_batch"))
  check("unneeded count batches are not read", length(list.files(events)) - before == 1L)
  missing_only <- ann[transcript_id == "absent"]
  check("entirely absent gene retains reference coverage",
    identical(as.data.frame(h1_counts(data.table(), missing_only, cfg, s0)$coverage),
      as.data.frame(h1_prepare_coverage(missing_only, cfg, stage("all_absent")))))
  benchmark <- rbindlist(runs)
  print(benchmark)
  if (length(args) > 1L) {
    dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
    fwrite(data.table(check = checks, status = "pass"),
      file.path(args[2], "h1_parallel_prepare_checks.tsv"), sep = "\t")
    fwrite(benchmark, file.path(args[2], "h1_parallel_fixture_benchmark.tsv"), sep = "\t")
  }
  cat("Completed ", length(checks), " H1 parallel preparation software checks.\n", sep = "")
}
main()
