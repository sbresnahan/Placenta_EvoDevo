#!/usr/bin/env Rscript
# Fake package trees, local source archives and executable mocks only.
a <- commandArgs(TRUE)
root <- normalizePath(if (length(a)) a[1L] else ".")
source(file.path(root, "R/library_paths.R"))
main <- function() {
  checks <- character()
  check <- function(label, x) {
    if (!isTRUE(x)) stop("FAIL ", label)
    checks <<- c(checks, label); cat("PASS ", label, "\n", sep = "")
  }
  fails <- function(x, text) {
    e <- tryCatch({force(x); NULL}, error = identity)
    if (!inherits(e, "error")) return(FALSE)
    if (!grepl(text, conditionMessage(e), fixed = TRUE)) message("Unexpected error: ", conditionMessage(e))
    grepl(text, conditionMessage(e), fixed = TRUE)
  }
  home <- tempfile("rstan-maintenance-"); dir.create(home)
  on.exit(unlink(home, recursive = TRUE), add = TRUE)
  e <- new.env(parent = globalenv())
  sys.source(file.path(root, "R/rstan_repair.R"), e)
  e$poe_h1_backend_check <- function(...) data.frame(status = "software-only")
  versions <- c(e$poe_rstan_versions(), brms = "2.23.0", cmdstanr = "0.99.0")
  make_package <- function(lib, p, version, label) {
    dir.create(file.path(lib, p, "libs"), recursive = TRUE, showWarnings = FALSE)
    write.dcf(data.frame(Package = p, Version = version), file.path(lib, p, "DESCRIPTION"))
    writeLines(label, file.path(lib, p, "libs", paste0(p, ".so")))
  }
  fixture <- function(name) {
    dir <- file.path(home, name); dir.create(dir)
    lib <- file.path(dir, "library"); dir.create(lib)
    for (p in names(versions)) make_package(lib, p, versions[[p]], paste("original", p))
    make_package(lib, "untouched", "1.0", "unchanged software sentinel")
    list(dir = dir, lib = lib, work = file.path(dir, "repair"),
      original = lapply(names(versions), function(p)
        readBin(file.path(lib, p, "libs", paste0(p, ".so")), "raw", 1000)))
  }
  calls <- new.env(); calls$build <- character(); calls$validate <- character()
  calls$fail <- ""; calls$bad_version <- FALSE; calls$space <- Inf
  api <- list(guard = function(library, allocated) TRUE,
    space = function(path) calls$space,
    download = function(package, version, dest) {
      src <- tempfile("archive-source-", tmpdir = home); dir.create(src)
      ver <- if (calls$bad_version) "0.0" else version
      make_package(src, package, ver, "not an executable or real native library")
      old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(src)
      utils::tar(dest, files = package, compression = "gzip", tar = "internal")
    },
    build = function(source, library, shared, makevars, log) {
      d <- read.dcf(file.path(source, "DESCRIPTION"))
      p <- unname(d[1L, "Package"])
      calls$build <- c(calls$build, p)
      if (calls$fail == paste0("build_", p)) stop("deliberate build failure")
      stopifnot(file.exists(makevars), dir.exists(shared))
      make_package(library, p, d[1L, "Version"], paste("rebuilt", p))
      writeLines("software build callback; no installation", log)
    },
    validate = function(library, shared, work, label) {
      calls$validate <- c(calls$validate, label)
      writeLines("software validation callback; no native loading", file.path(work, paste0(label, ".log")))
      for (p in names(e$poe_rstan_versions()))
        e$poe_rstan_description(p, library, e$poe_rstan_versions()[[p]])
      if (calls$fail == label) stop(paste("deliberate validation failure", label))
      TRUE
    })
  invoke <- function(f, action, confirmed = FALSE)
    e$poe_rstan_repair(root, action, f$work, confirmed, library = f$lib, .api = api)
  unchanged <- function(f) all(vapply(seq_along(versions), function(i) {
    p <- names(versions)[i]
    identical(f$original[[i]], readBin(file.path(f$lib, p, "libs", paste0(p, ".so")), "raw", 1000))
  }, TRUE))
  f <- fixture("normal")
  check("read-only check does not create a repair directory",
    nrow(invoke(f, "check")) == 1L && !file.exists(f$work) && unchanged(f))
  calls$space <- 1
  check("insufficient disk stops before work-directory writes",
    fails(invoke(f, "stage"), "Insufficient work-space") && !file.exists(f$work))
  calls$space <- Inf
  state <- invoke(f, "stage")
  check("pinned packages build in dependency order",
    identical(calls$build, names(e$poe_rstan_versions())))
  check("staging leaves all shared packages unchanged",
    state$phase == "staged" && unchanged(f))
  check("staged validation occurs before promotion",
    identical(calls$validate, "validation_staged"))
  check("staging refuses to replace existing work and logs",
    fails(invoke(f, "stage"), "new empty repair"))
  check("apply requires maintenance confirmation",
    fails(invoke(f, "apply"), "--maintenance-confirmed") && unchanged(f))
  state <- invoke(f, "apply", TRUE)
  check("promotion validates again using normal shared library",
    state$phase == "applied" && identical(tail(calls$validate, 2L),
      c("validation_before_apply", "validation_shared")))
  check("all original package directories are backed up",
    all(dir.exists(file.path(state$backup, "originals", names(e$poe_rstan_versions())))))
  check("promoted packages equal validated staged files",
    all(vapply(names(e$poe_rstan_versions()), function(p)
      e$poe_rstan_same(file.path(f$lib, p), file.path(state$staging, p)), TRUE)))
  check("unrelated package remains unchanged",
    readLines(file.path(f$lib, "untouched/libs/untouched.so")) == "unchanged software sentinel")
  check("repeat apply is refused", fails(invoke(f, "apply", TRUE), "successfully staged"))
  state <- invoke(f, "rollback", TRUE)
  check("explicit rollback restores all original package bytes",
    state$phase == "rolled_back" && unchanged(f))
  check("rollback preserves rejected rebuilt packages",
    length(list.files(file.path(state$backup, "rejected"))) == 3L)
  check("maintenance lock is released after completion",
    !dir.exists(file.path(f$lib, ".poe_rstan_repair.lock")))

  f <- fixture("bad_source"); calls$bad_version <- TRUE
  check("wrong source version fails before any shared package change",
    fails(invoke(f, "stage"), "Unexpected installed version") && unchanged(f))
  calls$bad_version <- FALSE
  f <- fixture("build_error"); calls$fail <- "build_StanHeaders"
  check("build failure preserves originals and logs",
    fails(invoke(f, "stage"), "deliberate build failure") && unchanged(f) &&
      file.exists(file.path(f$work, "build_RcppParallel.log")))
  check("incomplete builds cannot be applied",
    fails(invoke(f, "apply", TRUE), "successfully staged") && unchanged(f))
  f <- fixture("stage_error"); calls$fail <- "validation_staged"
  check("staged validation failure does not modify shared library",
    fails(invoke(f, "stage"), "validation_staged") && unchanged(f))
  calls$fail <- ""
  f <- fixture("post_promotion_error"); invoke(f, "stage")
  calls$fail <- "validation_shared"
  check("post-promotion validation failure triggers automatic rollback",
    fails(invoke(f, "apply", TRUE), "validation_shared") && unchanged(f) &&
      readRDS(file.path(f$work, "state.rds"))$phase == "rolled_back")
  calls$fail <- ""
  f <- fixture("space_apply"); invoke(f, "stage"); calls$space <- 1
  check("incoming package space checked before original directories move",
    fails(invoke(f, "apply", TRUE), "Insufficient library space") && unchanged(f))
  calls$space <- Inf
  f <- fixture("state_guard"); invoke(f, "stage")
  state <- readRDS(file.path(f$work, "state.rds")); state$library <- "/not/selected"
  saveRDS(state, file.path(f$work, "state.rds"))
  check("repair state cannot redirect another library",
    fails(invoke(f, "apply", TRUE), "does not match") && unchanged(f))
  f <- fixture("changed_after_apply"); invoke(f, "stage"); state <- invoke(f, "apply", TRUE)
  writeLines("later external package change", file.path(f$lib, "rstan/libs/rstan.so"))
  check("rollback will not overwrite subsequent external package changes",
    fails(invoke(f, "rollback", TRUE), "changed after promotion") &&
      readLines(file.path(f$lib, "rstan/libs/rstan.so")) == "later external package change")
  check("rollback refusal leaves every other promoted package untouched",
    e$poe_rstan_same(file.path(f$lib, "RcppParallel"), file.path(state$staging, "RcppParallel")) &&
      e$poe_rstan_same(file.path(f$lib, "StanHeaders"), file.path(state$staging, "StanHeaders")) &&
      length(list.files(file.path(state$backup, "originals"))) == 3L)

  # Emulate a hard interruption between same-filesystem directory renames.
  f <- fixture("interrupted_promotion"); state <- invoke(f, "stage")
  state$backup <- file.path(f$lib, ".poe_rstan_repair_interrupted")
  for (d in c("originals", "incoming", "rejected"))
    dir.create(file.path(state$backup, d), recursive = TRUE)
  file.rename(file.path(f$lib, "RcppParallel"), file.path(state$backup, "originals/RcppParallel"))
  e$poe_rstan_copy(file.path(state$staging, "RcppParallel"), file.path(f$lib, "RcppParallel"))
  file.rename(file.path(f$lib, "StanHeaders"), file.path(state$backup, "originals/StanHeaders"))
  state$phase <- "promoting"; saveRDS(state, file.path(f$work, "state.rds"))
  lock <- file.path(f$lib, ".poe_rstan_repair.lock"); dir.create(lock)
  saveRDS(list(work = state$work), file.path(lock, "owner.rds"))
  check("rollback recovers an interrupted promotion with a missing active package",
    invoke(f, "rollback", TRUE)$phase == "rolled_back" && unchanged(f))
  check("interrupted maintenance lock is cleared only for explicit matching rollback",
    !dir.exists(lock))
  f <- fixture("foreign_lock"); invoke(f, "stage")
  lock <- file.path(f$lib, ".poe_rstan_repair.lock"); dir.create(lock)
  saveRDS(list(work = "/another/repair"), file.path(lock, "owner.rds"))
  check("another repair's lock is not removed",
    fails(invoke(f, "apply", TRUE), "Repair lock exists") &&
      fails(invoke(f, "rollback", TRUE), "Repair lock exists") && dir.exists(lock))
  env <- e$poe_rstan_clean_env("/stage", "/shared", "/private/Makevars")
  check("build environment neutralizes external TBB and compiler include paths",
    all(c("TBB_LIB", "TBB_INC", "TBB_LINK_LIB", "CPATH", "CPLUS_INCLUDE_PATH",
      "LIBRARY_PATH", "LD_LIBRARY_PATH", "LD_PRELOAD") %in% env))
  check("build environment restricts compilation and native threads",
    all(c("MAKEFLAGS=-j1", "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1") %in% env))
  check("build environment keeps staging before shared dependencies",
    all(c("R_LIBS_USER=/stage", "R_LIBS=/shared", "R_MAKEVARS_USER=/private/Makevars",
      "R_MAKEVARS_SITE=/dev/null") %in% env))
  check("default production guard requires the target R version",
    as.character(getRversion()) == "4.3.1" ||
      fails(e$poe_rstan_guard(home), "R 4.3.1"))
  check("container wrapper forwards allocation identity",
    any(grepl("LSB_MAX_NUM_PROCESSORS LSB_JOBID", readLines(file.path(root, "bin/Rscript_sif")), fixed = TRUE)))
  source <- file.path(home, "dependency-source"); dir.create(source)
  write.dcf(data.frame(Package = "software", Version = "1.0",
    Imports = "SoftwareMissingDependency"), file.path(source, "DESCRIPTION"))
  check("missing prerequisites stop without automatic installation",
    fails(e$poe_rstan_prerequisites(source, f$lib), "Missing prerequisite"))
  write.dcf(data.frame(Package = "software", Version = "1.0",
    Imports = "rstan (>= 99.0)"), file.path(source, "DESCRIPTION"))
  check("prerequisite version requirements remain enforced",
    fails(e$poe_rstan_prerequisites(source, f$lib), "Unsatisfied prerequisite"))
  before <- Sys.getenv("TBB_LIB", unset = NA_character_)
  Sys.setenv(TBB_LIB = "/external/tbb/software-fixture")
  log <- file.path(home, "clean_subprocess.log")
  e$poe_rstan_run("/bin/sh", c("-c",
    'test -z "${TBB_LIB+x}" && test "$MAKEFLAGS" = "-j1" && test "$R_LIBS_USER" = "/stage" && printf "CLEAN_SUBPROCESS_OK\\n"'),
    env, log)
  if (is.na(before)) Sys.unsetenv("TBB_LIB") else Sys.setenv(TBB_LIB = before)
  check("real clean subprocess applies environment guards",
    identical(readLines(log), "CLEAN_SUBPROCESS_OK"))
  # Login-node fetching uses base R only; no allocation, compiler or native load.
  f <- fixture("fetch_only")
  calls$download <- character()
  fetch_api <- list(download = function(package, version, dest) {
    calls$download <- c(calls$download, package)
    api$download(package, version, dest)
  }, guard = function(...) stop("guard must not run during fetch"),
    build = function(...) stop("build must not run during fetch"),
    validate = function(...) stop("native validation must not run during fetch"))
  e$poe_h1_backend_check <- function(...) stop("backend must not load during fetch")
  old_job <- Sys.getenv("LSB_JOBID", unset = NA_character_)
  Sys.unsetenv("LSB_JOBID")
  fetched <- e$poe_rstan_repair(root, "fetch", f$work,
    library = "/not/an/installed/library", .api = fetch_api)
  if (is.na(old_job)) Sys.unsetenv("LSB_JOBID") else Sys.setenv(LSB_JOBID = old_job)
  e$poe_h1_backend_check <- function(...) data.frame(status = "software-only")
  source_cache <- f$work
  expected_archives <- paste0(names(e$poe_rstan_versions()), "_",
    e$poe_rstan_versions(), ".tar.gz")
  check("fetch needs neither allocation nor installed packages nor backend loading",
    nrow(fetched) == 3L && unchanged(f))
  check("fetch retrieves exactly the three pinned sources",
    identical(calls$download, names(e$poe_rstan_versions())) &&
      identical(fetched$version, unname(e$poe_rstan_versions())) &&
      all(file.exists(file.path(source_cache, expected_archives))))
  check("fetch validates sources without creating a staging library",
    file.exists(file.path(source_cache, "sources.tsv")) &&
      !dir.exists(file.path(source_cache, "stage_library")))
  no_network <- fetch_api
  no_network$download <- function(...) stop("unexpected network request")
  fetched_again <- e$poe_rstan_repair(root, "fetch", source_cache, .api = no_network)
  check("fetch reuses validated archives without another network call",
    all(fetched_again$status == "reused"))
  check("fetch rejects relative source directories",
    fails(e$poe_rstan_repair(root, "fetch", "relative-cache", .api = no_network),
      "absolute source directory"))
  check("sources option is rejected outside staging",
    fails(e$poe_rstan_repair(root, "fetch", source_cache, .api = no_network,
      sources = source_cache), "only valid for stage"))
  # Exercise the actual CLI with cached software archives and no scheduler job.
  cli_log <- file.path(home, "fetch_cli.log")
  cli_rc <- system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c("--vanilla", file.path(root, "bin/repair_rstan_stack.R"), "fetch", source_cache)),
    stdout = cli_log, stderr = cli_log,
    env = c(paste0("POE_R_LIBRARY=", shQuote(.libPaths()[1L])), "LSB_JOBID="))
  check("real fetch CLI works without allocation or native backend packages",
    cli_rc == 0L &&
      any(grepl("No compilation or installed-package changes", readLines(cli_log), fixed = TRUE)))
  cli_rc <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c("--vanilla", file.path(root, "bin/repair_rstan_stack.R"),
      "stage", file.path(home, "unused-cli-work"), "--sources")),
    stdout = cli_log, stderr = cli_log,
    env = paste0("POE_R_LIBRARY=", shQuote(.libPaths()[1L]))))
  check("CLI rejects an incomplete sources option before maintenance",
    cli_rc != 0L && any(grepl("Usage:", readLines(cli_log), fixed = TRUE)) &&
      !dir.exists(file.path(home, "unused-cli-work")))
  f <- fixture("offline_stage")
  offline_api <- api
  offline_api$download <- no_network$download
  calls$build <- character(); calls$validate <- character()
  state <- e$poe_rstan_repair(root, "stage", f$work, library = f$lib,
    .api = offline_api, sources = source_cache)
  check("offline staging builds in order without any downloader call",
    state$phase == "staged" && identical(calls$build, names(e$poe_rstan_versions())))
  check("offline staging preserves shared packages and validates staged backend",
    unchanged(f) && identical(calls$validate, "validation_staged"))
  check("offline archives are copied into retained work sources",
    all(vapply(expected_archives, function(name) {
      from <- file.path(source_cache, name); to <- file.path(f$work, "sources", name)
      identical(readBin(from, "raw", file.info(from)$size),
        readBin(to, "raw", file.info(to)$size))
    }, TRUE)))
  state <- e$poe_rstan_repair(root, "apply", f$work, TRUE, library = f$lib, .api = offline_api)
  check("offline staged repair retains guarded promotion and normal-library validation",
    state$phase == "applied" && tail(calls$validate, 1L) == "validation_shared")
  state <- e$poe_rstan_repair(root, "rollback", f$work, TRUE, library = f$lib, .api = offline_api)
  check("offline staged repair retains original-package rollback",
    state$phase == "rolled_back" && unchanged(f))
  for (kind in c("missing_last", "wrong_last", "invalid_last", "empty_last")) {
    f <- fixture(paste0("offline_", kind))
    cache <- file.path(f$dir, "cache"); dir.create(cache)
    stopifnot(all(file.copy(file.path(source_cache, expected_archives), cache)))
    last <- file.path(cache, tail(expected_archives, 1L))
    if (kind == "missing_last") unlink(last)
    if (kind == "wrong_last") {
      calls$bad_version <- TRUE
      api$download("rstan", e$poe_rstan_versions()[["rstan"]], last)
      calls$bad_version <- FALSE
    }
    if (kind == "invalid_last") writeLines("not a source archive", last)
    if (kind == "empty_last") writeBin(raw(), last)
    calls$build <- character(); calls$validate <- character()
    failed <- suppressWarnings(tryCatch({
      e$poe_rstan_repair(root, "stage", f$work, library = f$lib,
        .api = offline_api, sources = cache); NULL
    }, error = identity))
    check(paste("offline", kind, "source stops before the first build or shared change"),
      inherits(failed, "error") && !length(calls$build) &&
        !length(calls$validate) && unchanged(f))
  }
  f <- fixture("fetch_wrong_version"); calls$bad_version <- TRUE
  check("fetch rejects wrong-version downloads before publishing an archive",
    fails(e$poe_rstan_repair(root, "fetch", f$work, .api = fetch_api),
      "Unexpected installed version") &&
      !file.exists(file.path(f$work, expected_archives[1L])) && unchanged(f))
  calls$bad_version <- FALSE
  f <- fixture("fetch_invalid")
  invalid_api <- fetch_api
  invalid_api$download <- function(package, version, dest) writeLines("incomplete archive", dest)
  check("fetch rejects invalid download content before publishing an archive",
    suppressWarnings(fails(e$poe_rstan_repair(root, "fetch", f$work, .api = invalid_api),
      "Unexpected source archive layout")) &&
      !file.exists(file.path(f$work, expected_archives[1L])) && unchanged(f))
  # A failed low-level transfer never leaves a completed destination archive.
  dest <- file.path(home, "download-result.tar.gz")
  attempts <- 0L
  check("nonzero download status never publishes partial data",
    fails(e$poe_rstan_download("software", "1.0", dest, .download = function(url, path, ...) {
      attempts <<- attempts + 1L; writeLines("partial data", path); 1L
    }), "check network access") && attempts == 2L && !file.exists(dest) &&
      !length(list.files(home, pattern = "[.]partial$")))
  attempts <- 0L
  check("download exceptions preserve network cause and remove partial data",
    fails(e$poe_rstan_download("software", "1.0", dest, .download = function(url, path, ...) {
      attempts <<- attempts + 1L; writeLines("partial data", path); stop("network timeout fixture")
    }), "network timeout fixture") && attempts == 2L && !file.exists(dest))
  attempts <- 0L
  e$poe_rstan_download("software", "1.0", dest, .download = function(url, path, ...) {
    attempts <<- attempts + 1L
    writeLines(if (attempts == 1L) "partial" else "completed transfer", path)
    if (attempts == 1L) 1L else 0L
  })
  check("successful fallback publishes only the successful transfer",
    attempts == 2L && identical(readLines(dest), "completed transfer"))
  check("download refuses to overwrite an existing destination",
    fails(e$poe_rstan_download("software", "1.0", dest, .download = function(url, path, ...) {
      writeLines("replacement", path); 0L
    }), "without replacing an existing file") &&
      identical(readLines(dest), "completed transfer"))
  if (length(a) > 1L) {
    dir.create(a[2L], recursive = TRUE, showWarnings = FALSE)
    write.table(data.frame(check = checks, status = "pass"),
      file.path(a[2L], "rstan_repair_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  cat("Completed ", length(checks), " repair software checks. No packages installed.\n", sep = "")
}
main()
