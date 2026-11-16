#!/usr/bin/env Rscript
# Hand-specified aggregate software records; no biological estimates or fitting.
a <- commandArgs(TRUE)
root <- normalizePath(if (length(a)) a[1] else ".")
suppressPackageStartupMessages(library(data.table))
source(file.path(root, "R/exposure_archive.R"))
checks <- character()
check <- function(label, expr) {
  if (!isTRUE(expr)) stop("FAIL ", label)
  checks <<- c(checks, label); cat("PASS ", label, "\n", sep = "")
}
fails <- function(expr) inherits(tryCatch({force(expr); NULL}, error = function(e) e), "error")
home <- tempfile("archive-contract-"); dir.create(home)
fixture <- function(id = "GUSTO", empty = FALSE) {
  run <- tempfile("run-", tmpdir = home); dir.create(run)
  snap <- file.path(run, "aggregate_export"); dir.create(snap)
  d <- data.table(job_id = c("job_00000001", "job_00000002"), module = "h6",
    cohort_id = id, level = "gene", feature_id = c("software_A", "software_B"),
    family_id = "native_TE_family", status = c("unavailable", "failed_convergence_or_Hessian"))
  if (empty) d <- d[0]
  fwrite(d, file.path(snap, "model_results.tsv"), sep = "\t")
  fwrite(data.table(module = "h6", planned = nrow(d)), file.path(snap, "model_coverage.tsv"), sep = "\t")
  fwrite(data.table(module = "h6", jobs = nrow(d)), file.path(snap, "job_inventory.tsv"), sep = "\t")
  fwrite(data.table(version = "exposure_nb_v1", cohort_id = id, species = "software_fixture",
    build = "software_fixture"), file.path(snap, "provenance.tsv"), sep = "\t")
  config <- file.path(run, "config.yaml")
  # Deliberately not a full schema-2 config and no baseline/private paths exist.
  yaml::write_yaml(list(cohort_id = id, exposure_nb = list(output_dir = ".")), config)
  list(id = id, snapshot = snap, config = config)
}
refresh <- function(x) {
  fs <- list.files(x$snapshot, recursive = TRUE)
  fs <- setdiff(fs, "MANIFEST.tsv")
  fwrite(data.table(file = fs, bytes = file.info(file.path(x$snapshot, fs))$size),
    file.path(x$snapshot, "MANIFEST.tsv"), sep = "\t")
}
reject <- function(label, change) {
  x <- fixture(); refresh(x); change(x)
  check(label, fails(ex_archive_snapshot(x)))
}
x <- fixture(); refresh(x)
out <- file.path(home, "single.tar.gz")
ex_archive(x$config, out, root)
check("single cohort archive created", file.info(out)$size > 0)
unpack <- file.path(home, "extract"); dir.create(unpack); utils::untar(out, exdir = unpack)
base <- file.path(unpack, "single", "cohorts", "GUSTO", "exposure_nb_v1")
result <- fread(file.path(base, "model_results.tsv"))
check("native TE family retained", all(result$family_id == "native_TE_family"))
check("failed and unavailable rows retained", nrow(result) == 2L && all(result$status != "eligible"))
check("inventories use path and size only", identical(names(fread(file.path(base, "MANIFEST.tsv"))), c("file", "bytes")))
av <- fread(file.path(unpack, "single", "AVAILABILITY.tsv"))
check("optional absences explicit", all(av$state == "absent_from_snapshot"))
paths <- utils::untar(out, list = TRUE)
check("private products excluded", !any(grepl("private|prepared[.]rds|config[.]yaml|initial_report", paths)))
check("overwrite refused", fails(ex_archive(x$config, out, root)))
y <- fixture("MOUSE", TRUE); refresh(y)
ex_archive(c(x$config, y$config), file.path(home, "multiple.tar.gz"), root)
check("multiple cohorts archived", sum(grepl("/provenance.tsv$",
  utils::untar(file.path(home, "multiple.tar.gz"), list = TRUE))) == 2L)
check("header only model table accepted", nrow(ex_archive_snapshot(y)) == 4L)
check("duplicate cohorts rejected", fails(ex_archive(c(x$config, x$config), file.path(home, "dup.tar.gz"), root)))
# Unlisted stale files must never be collected.
fwrite(data.table(sample_id = "never_share"), file.path(x$snapshot, "stale.tsv"), sep = "\t")
ex_archive(x$config, file.path(home, "stale_safe.tar.gz"), root)
check("unlisted stale files excluded", !any(grepl("stale.tsv", utils::untar(file.path(home, "stale_safe.tar.gz"), list = TRUE), fixed = TRUE)))
# Extra old metadata is not interpreted or re-exported.
m <- fread(file.path(y$snapshot, "MANIFEST.tsv")); m[, obsolete_metadata := "ignored"]
fwrite(m, file.path(y$snapshot, "MANIFEST.tsv"), sep = "\t")
check("extra manifest metadata ignored", identical(names(ex_archive_snapshot(y)), c("file", "bytes")))
z <- fixture(); dir.create(file.path(z$snapshot, "bootstrap"))
fwrite(data.table(selection_index = 1L, bootstrap_status = "unavailable"),
  file.path(z$snapshot, "bootstrap", "job_00000001.tsv"), sep = "\t")
refresh(z)
ex_archive(z$config, file.path(home, "optional.tar.gz"), root)
check("aggregate bootstrap summary accepted", any(grepl("/bootstrap/job_00000001.tsv",
  utils::untar(file.path(home, "optional.tar.gz"), list = TRUE), fixed = TRUE)))
reject("missing core file rejected", function(x) unlink(file.path(x$snapshot, "provenance.tsv")))
reject("missing manifest rejected", function(x) unlink(file.path(x$snapshot, "MANIFEST.tsv")))
reject("wrong declared size rejected", function(x) {
  m <- fread(file.path(x$snapshot, "MANIFEST.tsv")); m$bytes[1] <- m$bytes[1] + 1
  fwrite(m, file.path(x$snapshot, "MANIFEST.tsv"), sep = "\t")
})
reject("duplicate paths rejected", function(x) {
  m <- fread(file.path(x$snapshot, "MANIFEST.tsv"))
  fwrite(rbind(m, m[1]), file.path(x$snapshot, "MANIFEST.tsv"), sep = "\t")
})
for (bad in c("../private.tsv", "/absolute.tsv", "bootstrap/../private.tsv", "bootstrap\\evil.tsv")) {
  reject(paste("unsafe path rejected", bad), function(x) {
    m <- fread(file.path(x$snapshot, "MANIFEST.tsv")); m$file[1] <- bad
    fwrite(m, file.path(x$snapshot, "MANIFEST.tsv"), sep = "\t")
  })
}
reject("symlinked file rejected", function(x) {
  p <- file.path(x$snapshot, "provenance.tsv"); backup <- paste0(p, ".outside")
  file.rename(p, backup); stopifnot(file.symlink(backup, p))
})
reject("symlinked parent rejected", function(x) {
  p <- x$snapshot; backup <- paste0(p, "_outside")
  file.rename(p, backup); stopifnot(file.symlink(backup, p))
})
for (column in c("SubjectID", "SAMPLE_ID", "cluster_id", "maternal_count", "PFOS_mat", "E", "input_file")) {
  reject(paste("private column rejected", column), function(x) {
    p <- file.path(x$snapshot, "model_results.tsv"); d <- fread(p); d[, (column) := "private"]
    fwrite(d, p, sep = "\t"); refresh(x)
  })
}
reject("unknown file rejected", function(x) {
  fwrite(data.table(a = 1L), file.path(x$snapshot, "unapproved.tsv"), sep = "\t"); refresh(x)
})
reject("missing required table columns rejected", function(x) {
  fwrite(data.table(a = 1L), file.path(x$snapshot, "job_inventory.tsv"), sep = "\t"); refresh(x)
})
reject("cohort mismatch rejected", function(x) {
  p <- file.path(x$snapshot, "provenance.tsv"); d <- fread(p); d[, cohort_id := "OTHER"]
  fwrite(d, p, sep = "\t"); refresh(x)
})
reject("duplicate result jobs rejected", function(x) {
  p <- file.path(x$snapshot, "model_results.tsv"); d <- fread(p)
  fwrite(rbind(d, d[1]), p, sep = "\t"); refresh(x)
})
invoke <- function(args) {
  text <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(file.path(root, "bin/export_results.R")), shQuote(args)), stdout = TRUE, stderr = TRUE))
  list(status = if (is.null(attr(text, "status"))) 0L else attr(text, "status"), text = text)
}
r <- invoke(c("--exposure-nb", "--out", file.path(home, "cli.tar.gz"), y$config))
if (r$status != 0L) stop(paste(r$text, collapse = "\n"))
check("CLI independent mode works without private inputs", r$status == 0L)
for (flags in list(c("--exposure-nb", "--baseline-nb"), c("--exposure-nb", "--force"),
                  c("--exposure-nb", "--partial", "MOUSE"), c("--exposure-nb", "--partial", ""))) {
  check(paste("CLI rejects incompatible flags", paste(flags, collapse = " ")),
    invoke(c(flags, "--out", file.path(home, "bad_cli.tar.gz"), y$config))$status != 0L)
}
if (length(a) > 1L) {
  dir.create(a[2], recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(check = checks, status = "pass"), file.path(a[2], "exposure_archive_tests.tsv"), sep = "\t")
}
unlink(home, recursive = TRUE)
cat("Completed ", length(checks), " archive software checks.\n", sep = "")
