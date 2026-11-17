# Explicit operator maintenance only. No package operations occur on source().
poe_rstan_versions <- function() c(RcppParallel = "5.1.11-2",
  StanHeaders = "2.32.10", rstan = "2.32.7")
poe_rstan_save <- function(x, path) {
  tmp <- tempfile("repair-state-", tmpdir = dirname(path))
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Cannot record repair state: ", path)
  invisible(x)
}
poe_rstan_files <- function(path) {
  z <- list.files(path, recursive = TRUE, full.names = FALSE, all.files = TRUE, no.. = TRUE)
  sort(z[!dir.exists(file.path(path, z))])
}
poe_rstan_bytes <- function(path) sum(file.info(file.path(path, poe_rstan_files(path)))$size)
poe_rstan_entries <- function(path) {
  c(".", sort(list.files(path, recursive = TRUE, all.files = TRUE,
    no.. = TRUE, include.dirs = TRUE)))
}
poe_rstan_same <- function(a, b) {
  if (!dir.exists(a) || !dir.exists(b)) return(FALSE)
  ea <- poe_rstan_entries(a); eb <- poe_rstan_entries(b)
  if (!identical(ea, eb)) return(FALSE)
  pa <- file.path(a, ea); pb <- file.path(b, eb)
  if (!identical(dir.exists(pa), dir.exists(pb)) ||
      !identical(bitwAnd(as.integer(file.info(pa)$mode), 511L),
        bitwAnd(as.integer(file.info(pb)$mode), 511L))) return(FALSE)
  fa <- ea[!dir.exists(pa)]
  all(vapply(fa, function(f) {
    pa <- file.path(a, f); pb <- file.path(b, f)
    n <- file.info(pa)$size
    identical(n, file.info(pb)$size) &&
      identical(readBin(pa, "raw", n), readBin(pb, "raw", n))
  }, TRUE))
}
poe_rstan_copy <- function(from, to) {
  if (file.exists(to)) stop("Refusing to replace an existing maintenance directory: ", to)
  if (!dir.create(to, recursive = TRUE, mode = "0700")) stop("Cannot create ", to)
  fs <- list.files(from, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (!all(file.copy(fs, to, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)))
    stop("Package content copy failed: ", from)
  # file.copy applies the process umask even with copy.mode=TRUE.
  # Restore exact source permissions inside the private backup parent, including
  # the package root that will later become visible in the shared library.
  entries <- poe_rstan_entries(from)
  modes <- file.info(file.path(from, entries))$mode
  if (!all(Sys.chmod(file.path(to, entries), modes, use_umask = FALSE)))
    stop("Cannot preserve copied package permissions: ", from)
  if (!poe_rstan_same(from, to)) stop("Package copy verification failed: ", from)
  invisible(to)
}
poe_rstan_description <- function(package, library, expected) {
  path <- file.path(library, package, "DESCRIPTION")
  if (!file.exists(path)) stop("Missing installed package: ", package, " in ", library)
  d <- read.dcf(path)
  if (nrow(d) != 1L || !identical(unname(d[1L, "Package"]), package) ||
      utils::compareVersion(d[1L, "Version"], expected) != 0L)
    stop("Unexpected installed version of ", package, "; expected ", expected)
  invisible(d)
}
poe_rstan_prerequisites <- function(source, libraries) {
  d <- read.dcf(file.path(source, "DESCRIPTION"))
  fields <- intersect(c("Depends", "Imports", "LinkingTo"), colnames(d))
  if (!length(fields)) return(invisible(TRUE))
  deps <- trimws(unlist(strsplit(paste(d[1L, fields], collapse = ","), ",", fixed = TRUE)))
  for (dep in deps[nzchar(deps)]) {
    package <- sub("[[:space:](].*$", "", dep)
    if (package == "R") actual <- as.character(getRversion()) else {
      paths <- file.path(unique(libraries), package, "DESCRIPTION")
      paths <- paths[file.exists(paths)]
      if (!length(paths)) stop("Missing prerequisite ", package, "; no dependencies will be installed automatically")
      actual <- read.dcf(paths[1L])[1L, "Version"]
    }
    constraint <- regmatches(dep, regexec("\\((>=|<=|==|!=|>|<)[[:space:]]*([^)]*)\\)", dep))[[1L]]
    if (length(constraint)) {
      cmp <- utils::compareVersion(actual, trimws(constraint[3L]))
      ok <- switch(constraint[2L], ">=" = cmp >= 0, "<=" = cmp <= 0,
        "==" = cmp == 0, "!=" = cmp != 0, ">" = cmp > 0, "<" = cmp < 0)
      if (!ok) stop("Unsatisfied prerequisite ", dep, "; installed version is ", actual,
        "; no automatic dependency upgrade is permitted")
    }
  }
  invisible(TRUE)
}
poe_rstan_guard <- function(library, allocated = FALSE) {
  if (as.character(getRversion()) != "4.3.1") stop("This repair requires the existing R 4.3.1 container")
  if (!dir.exists(library) || file.access(library, 7L) != 0L)
    stop("Selected library must be readable, writable and searchable: ", library)
  if (allocated && !nzchar(Sys.getenv("LSB_JOBID")))
    stop("Run maintenance on an allocated Seadragon node (LSB_JOBID is required)")
  for (tool in c("/usr/bin/gcc", "/usr/bin/g++"))
    if (file.access(tool, 1L) != 0L) stop("Missing container compiler: ", tool)
  message("Container compiler: ", system2("/usr/bin/g++", "--version", stdout = TRUE)[1L])
  locked <- list.files(library, "^00LOCK", all.files = TRUE)
  if (length(locked)) stop("Existing package installation lock; coordinate maintenance: ",
    paste(locked, collapse = ", "))
  invisible(TRUE)
}
poe_rstan_clean_env <- function(library, shared, makevars = "/dev/null") {
  clear <- c("LD_LIBRARY_PATH", "LD_PRELOAD", "LIBRARY_PATH", "CPATH", "C_INCLUDE_PATH",
    "CPLUS_INCLUDE_PATH", "OBJC_INCLUDE_PATH", "TBB_LIB", "TBB_INC", "TBB_LINK_LIB",
    "TBB_ROOT", "TBBROOT", "RCPP_PARALLEL_BACKEND", "CONDA_PREFIX", "CONDA_DEFAULT_ENV",
    "CMAKE_PREFIX_PATH", "PKG_CONFIG_PATH", "CC", "CXX", "CXX11", "CXX14", "CXX17",
    "CPPFLAGS", "CFLAGS", "CXXFLAGS", "LDFLAGS", "MAKE", "MAKEFLAGS",
    "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "R_MAKEVARS_USER", "R_MAKEVARS_SITE")
  c(as.vector(rbind("-u", clear)),
    "PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
    paste0("R_LIBS_USER=", library), paste0("R_LIBS=", shared),
    paste0("POE_R_LIBRARY=", library), paste0("R_MAKEVARS_USER=", makevars),
    "R_MAKEVARS_SITE=/dev/null", "R_ENVIRON_USER=/dev/null", "R_ENVIRON=/dev/null",
    "R_PROFILE_USER=/dev/null", "R_PROFILE=/dev/null", "MAKE=/usr/bin/make",
    "MAKEFLAGS=-j1", "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1")
}
poe_rstan_run <- function(command, args, env, log) {
  rc <- system2("env", shQuote(c(env, command, args)), stdout = log, stderr = log)
  if (rc != 0L) stop("Maintenance subprocess failed; inspect ", log)
  invisible(TRUE)
}
poe_rstan_space <- function(path) {
  while (!dir.exists(path)) path <- dirname(path)
  out <- system2("df", c("-Pk", shQuote(path)), stdout = TRUE)
  last <- strsplit(trimws(tail(out, 1L)), "[[:space:]]+")[[1L]]
  n <- suppressWarnings(as.numeric(last[4L]) * 1024)
  if (!is.finite(n)) stop("Cannot determine available disk space for ", path)
  n
}
poe_rstan_download <- function(package, version, dest, .download = utils::download.file) {
  name <- paste0(package, "_", version, ".tar.gz")
  urls <- c(paste0("https://cran.r-project.org/src/contrib/Archive/", package, "/", name),
    paste0("https://cran.r-project.org/src/contrib/", name))
  partial <- tempfile("source-download-", tmpdir = dirname(dest), fileext = ".partial")
  on.exit(unlink(partial), add = TRUE)
  errors <- character()
  for (url in urls) {
    unlink(partial)
    ok <- tryCatch({
      rc <- .download(url, partial, mode = "wb", quiet = FALSE)
      identical(as.integer(rc), 0L) && file.exists(partial) && file.info(partial)$size > 0
    }, error = function(e) {
      errors <<- c(errors, conditionMessage(e)); FALSE
    })
    if (isTRUE(ok)) {
      if (file.exists(dest) || !file.rename(partial, dest))
        stop("Cannot preserve downloaded source without replacing an existing file: ", dest)
      return(invisible(dest))
    }
  }
  stop("Unable to retrieve exact source: ", name,
    "; check network access and source availability; no version substitution is permitted",
    if (length(errors)) paste0("\n", paste(unique(errors), collapse = "\n")) else "")
}
poe_rstan_unpack <- function(tarfile, package, version, dest) {
  if (!file.exists(tarfile) || dir.exists(tarfile) || file.info(tarfile)$size == 0 ||
      file.access(tarfile, 4L) != 0L) stop("Missing or unreadable source archive: ", tarfile)
  entries <- utils::untar(tarfile, list = TRUE)
  if (!length(entries) || any(!startsWith(entries, paste0(package, "/"))) ||
      any(grepl("(^|/)\\.\\.(/|$)", entries)) || any(startsWith(entries, "/")))
    stop("Unexpected source archive layout: ", package)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  if (utils::untar(tarfile, exdir = dest) != 0L) stop("Source archive extraction failed: ", package)
  source <- file.path(dest, package)
  poe_rstan_description(package, dest, version)
  source
}
poe_rstan_fetch <- function(directory, .download = poe_rstan_download) {
  if (is.null(directory) || length(directory) != 1L || is.na(directory) ||
      !startsWith(directory, "/")) stop("Supply one absolute source directory")
  if (!dir.exists(directory) &&
      !dir.create(directory, recursive = TRUE, mode = "0700"))
    stop("Cannot create source directory: ", directory)
  directory <- normalizePath(directory, mustWork = TRUE)
  if (file.access(directory, 7L) != 0L) stop("Source directory is not writable: ", directory)
  scratch <- tempfile(".fetch-", tmpdir = directory)
  if (!dir.create(scratch, mode = "0700")) stop("Cannot create source validation directory")
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  versions <- poe_rstan_versions()
  rows <- lapply(names(versions), function(p) {
    name <- paste0(p, "_", versions[[p]], ".tar.gz")
    target <- file.path(directory, name)
    existing <- file.exists(target)
    candidate <- if (existing) target else file.path(scratch, name)
    message(if (existing) "Checking cached source " else "Fetching source ", name)
    if (!existing) .download(p, versions[[p]], candidate)
    poe_rstan_unpack(candidate, p, versions[[p]], file.path(scratch, "check"))
    if (!existing && (file.exists(target) || !file.rename(candidate, target)))
      stop("Cannot preserve validated source: ", target)
    data.frame(package = p, version = versions[[p]], archive = name,
      bytes = file.info(target)$size, status = if (existing) "reused" else "downloaded")
  })
  result <- do.call(rbind, rows)
  write.table(result, file.path(directory, "sources.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)
  message("All three exact source archives validated. No compilation or installed-package changes.")
  invisible(result)
}
poe_rstan_validate <- function(root, library, shared, report) {
  versions <- poe_rstan_versions()
  for (p in names(versions)) poe_rstan_description(p, library, versions[[p]])
  for (p in c("brms", "cmdstanr")) poe_rstan_description(p, shared,
    c(brms = "2.23.0", cmdstanr = "0.99.0")[[p]])
  x <- poe_h1_backend_check()
  for (p in names(versions)) if (!identical(normalizePath(find.package(p)),
      normalizePath(file.path(library, p)))) stop("Validation selected the wrong package path: ", p)
  print(x)
  links <- unlist(lapply(c("rstan", "RcppParallel"), function(p) {
    so <- system.file("libs", paste0(p, .Platform$dynlib.ext), package = p)
    if (!nzchar(so)) stop("Missing native package library: ", p)
    c(paste("PACKAGE", p, so), system2("ldd", shQuote(so), stdout = TRUE, stderr = TRUE))
  }))
  cat(paste(links, collapse = "\n"), "\n")
  if (any(grepl("leafcutter", links, ignore.case = TRUE)))
    stop("Rebuilt package linkage still resolves through LeafCutter; stop and inspect build logs")
  print(RcppParallel::tbbLibraryPath())
  write.table(x, report, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(TRUE)
}
poe_rstan_api <- function(root) list(
  guard = poe_rstan_guard,
  space = poe_rstan_space,
  download = poe_rstan_download,
  build = function(source, library, shared, makevars, log) {
    poe_rstan_run(file.path(R.home("bin"), "R"),
      c("CMD", "INSTALL", "--preclean", "--no-multiarch", paste0("--library=", library), source),
      poe_rstan_clean_env(library, shared, makevars), log)
  },
  validate = function(library, shared, work, label) {
    code <- paste(
      'a <- commandArgs(TRUE); root <- a[1]; source(file.path(root,"R/library_paths.R"));',
      'source(file.path(root,"R/h1_backend.R")); source(file.path(root,"R/rstan_repair.R"));',
      'poe_rstan_validate(root, a[2], a[3], a[4])')
    poe_rstan_run(file.path(R.home("bin"), "Rscript"),
      c("--vanilla", "-e", code, root, library, shared, file.path(work, paste0(label, ".tsv"))),
      poe_rstan_clean_env(library, shared), file.path(work, paste0(label, ".log")))
  })
poe_rstan_state <- function(work, library) {
  state <- readRDS(file.path(work, "state.rds"))
  if (!identical(state$format, "poe_rstan_repair_v1") ||
      !identical(state$library, normalizePath(library)) ||
      !identical(state$work, normalizePath(work)) ||
      !identical(state$versions, poe_rstan_versions()) ||
      !identical(state$staging, file.path(state$work, "stage_library")))
    stop("Repair state does not match this work directory, selected library or pinned versions")
  if (nzchar(state$backup) &&
      (!identical(dirname(state$backup), state$library) ||
       !startsWith(basename(state$backup), ".poe_rstan_repair_")))
    stop("Unexpected backup location in repair state")
  state
}
poe_rstan_restore <- function(state) {
  if (!nzchar(state$backup) || !dir.exists(state$backup)) stop("No package backup is available")
  # Check every selected directory before the first rollback rename.
  for (p in names(state$versions)) {
    old <- file.path(state$backup, "originals", p)
    target <- file.path(state$library, p)
    if (dir.exists(old) && file.exists(target) &&
        !poe_rstan_same(target, file.path(state$staging, p)))
      stop("Selected package changed after promotion; refusing to replace it automatically: ", p)
    if (!dir.exists(old) && !dir.exists(target))
      stop("Neither original nor active package exists: ", p)
  }
  for (p in names(state$versions)) {
    old <- file.path(state$backup, "originals", p)
    target <- file.path(state$library, p)
    if (dir.exists(old)) {
      if (file.exists(target)) {
        rejected <- tempfile(paste0(p, "_"), tmpdir = file.path(state$backup, "rejected"))
        if (!file.rename(target, rejected)) stop("Cannot preserve promoted package during rollback: ", p)
      }
      if (!file.rename(old, target)) stop("Cannot restore original package: ", p)
    } else if (!dir.exists(target)) stop("Neither original nor active package exists: ", p)
  }
  state$phase <- "rolled_back"
  poe_rstan_save(state, file.path(state$work, "state.rds"))
  message("Original package directories restored. This restores the prior state, not a healthy backend.")
  invisible(state)
}
poe_rstan_repair <- function(root, action = c("check", "stage", "apply", "rollback", "fetch"),
                             work = NULL, confirmed = FALSE, library = .libPaths()[1L],
                             .api = NULL, sources = NULL) {
  action <- match.arg(action)
  if (is.null(.api)) .api <- poe_rstan_api(root)
  if (!is.null(sources) && action != "stage") stop("--sources is only valid for stage")
  if (action == "fetch") return(poe_rstan_fetch(work, .api$download))
  library <- normalizePath(library, mustWork = TRUE)
  .api$guard(library, allocated = action != "check")
  versions <- poe_rstan_versions()
  if (action %in% c("check", "stage"))
    for (p in names(versions)) poe_rstan_description(p, library, versions[[p]])
  for (p in c("brms", "cmdstanr")) poe_rstan_description(p, library,
    c(brms = "2.23.0", cmdstanr = "0.99.0")[[p]])
  message("Selected shared library: ", library, "; R: ", getRversion())
  if (action == "check") return(poe_h1_backend_check())
  if (is.null(work) || length(work) != 1L || !startsWith(work, "/"))
    stop("Supply one absolute private repair work directory")
  if (action == "stage") {
    if (file.exists(work) && (!dir.exists(work) || length(list.files(work, all.files = TRUE, no.. = TRUE))))
      stop("Use a new empty repair work directory; existing logs and state are retained")
    required <- max(4 * 1024^3, 6 * sum(vapply(names(versions),
      function(p) poe_rstan_bytes(file.path(library, p)), 0)))
    available <- .api$space(dirname(work))
    message(sprintf("Build-space estimate %.2f GB; available %.2f GB; one compiler process",
      required / 1024^3, available / 1024^3))
    if (available < required) stop("Insufficient work-space estimate before staging")
    dir.create(work, recursive = TRUE, mode = "0700", showWarnings = FALSE)
    work <- normalizePath(work, mustWork = TRUE)
    state <- list(format = "poe_rstan_repair_v1", phase = "building",
      work = work, library = library, staging = file.path(work, "stage_library"),
      versions = versions, backup = "")
    dir.create(state$staging, mode = "0700")
    for (d in c("sources", "unpacked")) dir.create(file.path(work, d), mode = "0700")
    poe_rstan_save(state, file.path(work, "state.rds"))
    makevars <- file.path(work, "Makevars.clean")
    writeLines(c("CC=/usr/bin/gcc", paste0(c("CXX", "CXX11", "CXX14", "CXX17"), "=/usr/bin/g++"),
      paste0(c("CFLAGS", "CXXFLAGS", "CXX11FLAGS", "CXX14FLAGS", "CXX17FLAGS"), "=-O2 -g0"),
      "CPPFLAGS=", "LDFLAGS="), makevars)
    if (!is.null(sources)) {
      if (length(sources) != 1L || is.na(sources) || !startsWith(sources, "/") ||
          !dir.exists(sources)) stop("Supply one absolute existing --sources directory")
      sources <- normalizePath(sources, mustWork = TRUE)
      message("Offline source mode: using local archives only; no download fallback")
    }
    unpacked <- setNames(character(length(versions)), names(versions))
    # Prepare every source before the first compiler process.
    for (p in names(versions)) {
      message("Preparing source ", p, " ", versions[[p]])
      tar <- file.path(work, "sources", paste0(p, "_", versions[[p]], ".tar.gz"))
      if (is.null(sources)) .api$download(p, versions[[p]], tar) else {
        from <- file.path(sources, basename(tar))
        if (!file.exists(from) || dir.exists(from) || file.access(from, 4L) != 0L)
          stop("Missing or unreadable source archive: ", from, "; no download fallback")
        n <- file.info(from)$size
        if (!file.copy(from, tar, copy.mode = TRUE) ||
            !identical(file.info(tar)$size, n) ||
            !identical(readBin(from, "raw", n), readBin(tar, "raw", n)))
          stop("Local source copy verification failed: ", from)
      }
      unpacked[[p]] <- poe_rstan_unpack(tar, p, versions[[p]], file.path(work, "unpacked"))
    }
    for (p in names(versions)) {
      message("Staging ", p, " ", versions[[p]])
      source <- unpacked[[p]]
      poe_rstan_prerequisites(source, c(state$staging, library, .libPaths()))
      .api$build(source, state$staging, library, makevars, file.path(work, paste0("build_", p, ".log")))
      poe_rstan_description(p, state$staging, versions[[p]])
    }
    .api$validate(state$staging, library, work, "validation_staged")
    state$phase <- "staged"
    poe_rstan_save(state, file.path(work, "state.rds"))
    message("Staged stack validated; shared packages are unchanged. Apply requires maintenance confirmation.")
    return(invisible(state))
  }
  if (!isTRUE(confirmed)) stop("Apply/rollback require --maintenance-confirmed; stop affected sessions first")
  state <- poe_rstan_state(work, library)
  lock <- file.path(library, ".poe_rstan_repair.lock")
  if (dir.exists(lock)) {
    owner <- tryCatch(readRDS(file.path(lock, "owner.rds")), error = function(e) NULL)
    if (action != "rollback" || !identical(owner$work, state$work))
      stop("Repair lock exists; do not run concurrent maintenance: ", lock)
    # Explicit rollback confirmation asserts the prior repair process is stopped.
    unlink(lock, recursive = TRUE)
  }
  if (!dir.create(lock, mode = "0700")) stop("Cannot acquire maintenance lock")
  saveRDS(list(work = state$work, pid = Sys.getpid(), host = Sys.info()[["nodename"]]),
    file.path(lock, "owner.rds"))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (action == "rollback") return(poe_rstan_restore(state))
  if (!identical(state$phase, "staged")) stop("Apply requires a successfully staged, unapplied repair")
  for (p in names(versions)) poe_rstan_description(p, library, versions[[p]])
  .api$validate(state$staging, library, state$work, "validation_before_apply")
  needed <- sum(vapply(names(versions), function(p)
    poe_rstan_bytes(file.path(state$staging, p)), 0)) + 256 * 1024^2
  if (.api$space(library) < needed) stop("Insufficient library space for verified incoming packages")
  state$backup <- file.path(library, paste0(".poe_rstan_repair_",
    format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid()))
  if (!dir.create(state$backup, mode = "0700")) stop("Cannot create private backup directory")
  for (d in c("originals", "incoming", "rejected"))
    dir.create(file.path(state$backup, d), mode = "0700")
  for (p in names(versions)) poe_rstan_copy(file.path(state$staging, p),
    file.path(state$backup, "incoming", p))
  state$phase <- "promoting"
  poe_rstan_save(state, file.path(state$work, "state.rds"))
  committed <- FALSE
  on.exit(if (!committed) {
    tryCatch(poe_rstan_restore(state), error = function(e)
      message("Automatic rollback needs operator attention: ", conditionMessage(e),
        "; retained state: ", file.path(state$work, "state.rds")))
  }, add = TRUE, after = FALSE)
  for (p in names(versions)) {
    target <- file.path(library, p)
    if (!file.rename(target, file.path(state$backup, "originals", p)))
      stop("Cannot preserve original package: ", p)
    if (!file.rename(file.path(state$backup, "incoming", p), target))
      stop("Cannot promote staged package: ", p)
  }
  .api$validate(library, library, state$work, "validation_shared")
  state$phase <- "applied"
  poe_rstan_save(state, file.path(state$work, "state.rds"))
  committed <- TRUE
  message("Shared stack validated. Backups retained at ", state$backup)
  invisible(state)
}
