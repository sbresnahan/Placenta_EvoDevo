#!/usr/bin/env Rscript
# Software-only backend and orchestration checks; no sampler is invoked.
a <- commandArgs(TRUE)
root <- normalizePath(if (length(a)) a[1L] else ".")
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages(library(data.table))
for (m in c("common", "config", "regulatory_common", "baseline_nb", "baseline_te", "h1_bayes"))
  source(file.path(root, "R", paste0(m, ".R")))
main <- function() {
  checks <- character()
  check <- function(label, x) {
    if (!isTRUE(x)) stop("FAIL ", label)
    checks <<- c(checks, label); cat("PASS ", label, "\n", sep = "")
  }
  fails <- function(x, text) {
    e <- tryCatch({force(x); NULL}, error = identity)
    inherits(e, "error") && grepl(text, conditionMessage(e), fixed = TRUE)
  }
  home <- tempfile("backend-recovery-"); dir.create(home)
  on.exit(unlink(home, recursive = TRUE), add = TRUE)
  calls <- character()
  api <- list(load = function(p) {calls <<- c(calls, p); TRUE},
    classes = function() TRUE, version = function(p) "software",
    path = function(p) "/software/package")
  check("backend preflight checks all six packages", nrow(poe_h1_backend_check(api)) == 6L &&
    identical(calls, c("rstan", "StanHeaders", "RcppParallel", "brms", "cmdstanr", "posterior")))
  bad <- api; bad$load <- function(p) stop("undefined TBB symbol software fixture")
  check("backend load error is actionable and stops before sampling",
    fails(poe_h1_backend_check(bad), "No new sampling started"))
  bad <- api; bad$classes <- function() stop("software class failure")
  check("backend validates RStan classes explicitly",
    fails(poe_h1_backend_check(bad), "RStan class validation failed"))
  e <- new.env(parent = globalenv()); e$root <- root
  sys.source(file.path(root, "R/h1_bayes.R"), e)
  e$h1_encode <- function(d) list(data = d, map = "software-map")
  runtime_calls <- 0L
  e$h1_cmdstan_runtime <- function(...) {runtime_calls <<- runtime_calls + 1L; stop("NO_SAMPLER")}
  e$poe_h1_backend_check <- function(...) stop("BACKEND_STOP")
  check("direct fitting stops before tool initialization when backend fails",
    fails(e$h1_fit_once(NULL, NULL, list(), list(), 1L, file.path(home, "missing.rds")),
      "BACKEND_STOP") && runtime_calls == 0L)
  e$poe_h1_backend_check <- function(...) TRUE
  check("direct fitting reaches tool initialization only after backend passes",
    fails(e$h1_fit_once(NULL, NULL, list(), list(), 1L, file.path(home, "missing.rds")),
      "NO_SAMPLER") && runtime_calls == 1L)
  saved <- file.path(home, "saved.rds")
  saveRDS(list(model_spec = e$H1_MODEL_SPEC, settings = list(), profile = list(),
    seed = 1L, map = "software-map", derived = NULL, fit = NULL,
    diagnostics = list(status = "diagnostic_failure", summary = "software")), saved)
  e$h1_summarize <- function(x) data.table(variable = "software")
  e$h1_residuals <- function(...) "software residuals"
  reused <- e$h1_fit_once(NULL, NULL, list(), list(), 1L, saved)
  check("saved diagnostic failure remains reusable without new tool initialization",
    runtime_calls == 1L && identical(reused$summary$status, "diagnostic_failure"))
  check("saved-task settings remain strictly checked",
    fails(e$h1_fit_once(NULL, NULL, list(chains = 1), list(), 1L, saved), "settings differ"))
  profile <- h1_profiles()[1L]
  map <- data.table(feature_id = c("T1", "T2"), transcript = c("tx000001", "tx000002"))
  e$h1_encode <- function(d) list(data = d, map = map)
  x <- readRDS(saved); x$profile <- profile; x$map <- map; saveRDS(x, saved)
  check("real table profiles and mappings survive serialization for reuse",
    !fails(e$h1_fit_once(NULL, NULL, list(), profile, 1L, saved), "settings differ"))
  changed <- copy(profile); changed[, effect_scale := 2]
  check("changed profile values are still rejected",
    fails(e$h1_fit_once(NULL, NULL, list(), changed, 1L, saved), "settings differ"))
  map <- map[2:1]
  check("changed transcript mapping order is still rejected",
    fails(e$h1_fit_once(NULL, NULL, list(), profile, 1L, saved), "settings differ"))
  check("table comparison retains column types",
    !h1_task_table_equal(data.table(x = 1L), data.table(x = 1)))

  # End-to-end gene scheduling with a fixed count-format fixture and fake fit API.
  cfg <- list(profile = "public_cross", dosage = list(maternal = 2, paternal = 1),
    paths = list(model_gene = file.path(home, "unused.tsv"), scratch = home),
    reference = list(build = "software"))
  attr(cfg, "file") <- file.path(home, "cfg.yaml")
  s <- list(output_dir = file.path(home, "gene"), min_snps = 3L, min_samples = 4L,
    min_total_counts = 50, min_per_cross = 2L)
  dir.create(s$output_dir)
  ann <- data.table(transcript_id = c("T1", "T2"), gene_id = "G",
    promoter_te_label = c("TE", "nonTE"))
  z <- CJ(feature_id = ann$transcript_id, sample_id = paste0("S", 1:4), j = 1:3)
  z[, `:=`(SNPpos = paste0(feature_id, ":", j), maternal_count = 4.4,
    paternal_count = 2.2, cluster_id = sample_id, independence_verified = TRUE,
    cross = ifelse(sample_id %in% c("S1", "S2"), "focal_mother", "alternate_mother"))]
  input <- file.path(home, "counts.rds"); saveRDS(list(data = z), input)
  nb_write(ann, nb_path(s, "transcript_TE_annotation.tsv"))
  nb_write(data.table(feature_id = ann$transcript_id, batch = 1L),
    nb_path(s, "transcript_feature_batches.tsv"))
  nb_write(data.table(batch = 1L, input_file = input), nb_path(s, "transcript_batches.tsv"))
  nb_write(data.table(gene_id = "G", job_index = 1L), nb_path(s, "h1_bayes_jobs.tsv"))
  nb_finish(s, "h1-inputs", list(software = TRUE), nb_path(s, "h1_bayes_jobs.tsv"))
  e <- new.env(parent = globalenv()); e$root <- root
  sys.source(file.path(root, "R/h1_bayes.R"), e)
  new_tasks <- reuse_tasks <- 0L; fail_one <- TRUE
  e$h1_cmdstan_runtime <- function(...) TRUE
  e$poe_h1_backend_check <- function(...) stop("BACKEND_STOP")
  e$h1_fit_once <- function(d, ann, hs, profile, seed, file) {
    if (file.exists(file)) {
      x <- readRDS(file); reuse_tasks <<- reuse_tasks + 1L
      stopifnot(identical(x$settings, hs), h1_task_table_equal(x$profile, profile), identical(x$seed, seed))
    } else {
      new_tasks <<- new_tasks + 1L
      if (profile$profile == "half_effect_scales" && fail_one) stop("deliberate task failure")
      status <- if (profile$profile == "primary") "diagnostic_failure" else "ready"
      x <- list(settings = hs, profile = profile, seed = seed, status = status)
      nb_save(x, file)
    }
    list(summary = data.table(variable = "software", status = x$status),
      diagnostics = data.table(status = x$status), residuals = data.table(note = "no real fit"))
  }
  receipt <- nb_receipt(s, "h1_000001")
  di_path <- nb_path(s, "h1_bayes/gene_000001_diagnostics.tsv")
  check("gene backend failure stops before task loop",
    fails(e$h1_gene(cfg, s, root, 1L), "BACKEND_STOP") &&
      new_tasks == 0L && !file.exists(receipt) && !file.exists(di_path))
  e$poe_h1_backend_check <- function(...) TRUE
  e$h1_gene(cfg, s, root, 1L)
  di <- reg_read(di_path)
  check("runtime failures and diagnostic failures remain distinct",
    sum(di$status == "fit_failure") == 1L && sum(di$status == "diagnostic_failure") == 1L)
  check("support-invalid omission tasks remain unavailable",
    sum(di$status == "unavailable_LOCO_support") == 4L && new_tasks == 5L)
  check("failed completed gene requires explicit retry",
    fails(e$h1_gene(cfg, s, root, 1L), "h1-retry-failed 1") && new_tasks == 5L)
  old_di <- readBin(di_path, "raw", file.info(di_path)$size)
  old_receipt <- readBin(receipt, "raw", file.info(receipt)$size)
  e$poe_h1_backend_check <- function(...) stop("BACKEND_STOP")
  check("failed retry preflight preserves existing receipt and output",
    fails(e$h1_gene(cfg, s, root, 1L, retry_failed = TRUE), "BACKEND_STOP") &&
      identical(old_di, readBin(di_path, "raw", file.info(di_path)$size)) &&
      identical(old_receipt, readBin(receipt, "raw", file.info(receipt)$size)))
  e$poe_h1_backend_check <- function(...) TRUE
  fail_one <- FALSE
  task_files <- list.files(nb_path(s, "private/h1_bayes"), full.names = TRUE)
  task_bytes <- lapply(task_files, function(p) readBin(p, "raw", file.info(p)$size))
  e$h1_gene(cfg, s, root, 1L, retry_failed = TRUE)
  check("explicit retry samples only the previously failed task",
    new_tasks == 6L && reuse_tasks == 4L && !any(reg_read(di_path)$status == "fit_failure"))
  check("valid saved task bytes are preserved by retry", all(vapply(seq_along(task_files),
    function(i) identical(task_bytes[[i]], readBin(task_files[i], "raw", file.info(task_files[i])$size)), TRUE)))
  backups <- list.files(nb_path(s, "private/h1_recovery"), full.names = TRUE)
  check("prior failed diagnostics and receipt are preserved",
    length(backups) == 1L &&
      identical(old_di, readBin(file.path(backups, basename(di_path)), "raw", length(old_di))) &&
      identical(old_receipt, readBin(file.path(backups, basename(receipt)), "raw", length(old_receipt))))
  e$h1_gene(cfg, s, root, 1L)
  check("completed diagnostic failure does not trigger automatic refitting",
    new_tasks == 6L && reuse_tasks == 4L)
  # Explicitly create an interrupted software state by moving its receipt/task.
  file.rename(receipt, file.path(home, "interrupted_receipt.rds"))
  first_task <- list.files(nb_path(s, "private/h1_bayes"), full.names = TRUE)[1L]
  file.rename(first_task, file.path(home, "interrupted_task.rds"))
  e$h1_gene(cfg, s, root, 1L)
  check("normal restart resumes interrupted state without refitting saved tasks",
    new_tasks == 7L && reuse_tasks == 8L)
  check("retry CLI and scheduler action are wired",
    any(grepl("`h1-retry-failed` =", readLines(file.path(root, "poe_baseline.R")), fixed = TRUE)) &&
      any(grepl("|h1-retry-failed|", readLines(file.path(root, "bin/submit_baseline_nb.sh")), fixed = TRUE)))
  if (length(a) > 1L) {
    dir.create(a[2L], recursive = TRUE, showWarnings = FALSE)
    fwrite(data.table(check = checks, status = "pass"),
      file.path(a[2L], "h1_backend_recovery_checks.tsv"), sep = "\t")
  }
  cat("Completed ", length(checks), " backend/recovery software checks. No sampling.\n", sep = "")
}
main()
