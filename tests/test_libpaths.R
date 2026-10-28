#!/usr/bin/env Rscript
# Startup/argument fixtures only; no real container, scheduler, or cohort data.
f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
f <- gsub("~+~", " ", f, fixed = TRUE)
root <- normalizePath(file.path(dirname(f), ".."), mustWork = TRUE)

main <- function() {
  n <- 0L
  check <- function(x) { stopifnot(isTRUE(x)); n <<- n + 1L }
  rs <- file.path(R.home("bin"), "Rscript")
  d <- tempfile("R library startup fixtures ")
  dir.create(d)
  keys <- c("SCRIPTS_DIR", "POE_SITE_FILE", "POE_R_LIBRARY", "POE_R_BOOTSTRAP",
            "POE_NATIVE_RECORD", "LD_LIBRARY_PATH", "PYTHONPATH", "PYTHONHOME")
  old <- Sys.getenv(keys, unset = NA_character_)
  on.exit({
    for (k in keys) {
      if (is.na(old[[k]])) Sys.unsetenv(k) else
        do.call(Sys.setenv, setNames(list(old[[k]]), k))
    }
    unlink(d, recursive = TRUE)
  }, add = TRUE)
  put <- function(path, text, executable = FALSE) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(text, path)
    if (executable) Sys.chmod(path, "0755")
    path
  }
  run <- function(command, args) {
    out <- suppressWarnings(system2(command, vapply(args, shQuote, character(1)),
                                   stdout = TRUE, stderr = TRUE))
    rc <- attr(out, "status")
    list(rc = if (is.null(rc)) 0L else rc, out = out)
  }
  # A spaced script root and library catch shell/R interpolation mistakes.
  fixture_root <- file.path(d, "repository with spaces")
  dir.create(file.path(fixture_root, "R"), recursive = TRUE)
  dir.create(file.path(fixture_root, "bin"))
  stopifnot(file.copy(file.path(root, "R/library_paths.R"),
                     file.path(fixture_root, "R/library_paths.R")))
  stopifnot(file.copy(file.path(root, "bin/Rscript_sif"),
                     file.path(fixture_root, "bin/Rscript_sif")))
  lib <- file.path(d, "library with spaces")
  dir.create(lib)
  put(file.path(fixture_root, "config/r_library.path"), lib)
  image <- put(file.path(d, "fixture.sif"), "NOT_A_REAL_CONTAINER")
  mock <- put(file.path(d, "mock singularity"), c(
    "#!/bin/bash", "set -euo pipefail",
    '[[ -z "${LD_LIBRARY_PATH:-}" && -z "${PYTHONPATH:-}" && -z "${PYTHONHOME:-}" ]] || exit 71',
    '[[ "$1" == exec && "$2" == --cleanenv ]] || exit 72; shift 2',
    'pass=(); binds=()',
    'while [[ $# -gt 0 ]]; do case "$1" in',
    '  --bind) binds+=("$2"); shift 2 ;;',
    '  --env) pass+=("$2"); shift 2 ;;',
    '  *) shift; break ;;',
    'esac; done',
    '[[ "${binds[*]}" == "/rsrch5 /rsrch9 /tmp" && "$1" == Rscript ]] || exit 73; shift',
    '# Intentionally discard R_LIBS_USER: .libPaths must work without it.',
    paste('exec env -i "PATH=$PATH" "${pass[@]}" R_LIBS_USER=',
          shQuote(rs), '"$@"')
  ), TRUE)
  site <- put(file.path(d, "site.sh"), c(
    paste0("export POE_R_SIF=", shQuote(image)),
    paste0("export POE_R_LIBRARY=", shQuote(lib)),
    paste0("export POE_SINGULARITY=", shQuote(mock))
  ))
  Sys.setenv(SCRIPTS_DIR = fixture_root, POE_SITE_FILE = site,
             POE_NATIVE_RECORD = "fixture record with spaces",
             LD_LIBRARY_PATH = "/fixture/not/a/native/library",
             PYTHONPATH = "/fixture/python", PYTHONHOME = "/fixture/python/home")
  wrapper <- file.path(fixture_root, "bin/Rscript_sif")
  check(run("bash", c("-n", wrapper))$rc == 0L)
  expr <- paste0(
    'stopifnot(identical(.libPaths()[1L], normalizePath(Sys.getenv("POE_R_LIBRARY"))), ',
    '.Library %in% .libPaths(), ',
    'identical(Sys.getenv("POE_R_BOOTSTRAP"), file.path(',
    deparse(fixture_root), ', "R/library_paths.R")), ',
    'identical(Sys.getenv("POE_NATIVE_RECORD"), "fixture record with spaces")); ',
    'cat("EXPLICIT_LIBPATH_OK\\n")')
  for (options in list(character(), "--vanilla", "--no-init-file",
                       "--default-packages=NULL")) {
    z <- run("bash", c(wrapper, options, "-e", expr))
    check(z$rc == 0L && "EXPLICIT_LIBPATH_OK" %in% z$out)
  }
  # Leading options, several expressions, and trailing values are unchanged.
  args <- c("one space", "-e", "literal ' quote $HOME ;")
  code <- paste0("stopifnot(x == 7L, identical(commandArgs(TRUE), ",
                 paste(deparse(c("--args", args)), collapse = ""), ")); cat('ARGS_OK\\n')")
  z <- run("bash", c(wrapper, "--vanilla", "-e", "x <- 7L",
                     "-e", code, "--args", args))
  check(z$rc == 0L && "ARGS_OK" %in% z$out)
  z <- run("bash", c(wrapper, "-e", "quit(status=42L)"))
  check(z$rc == 42L)
  # Script mode retains --file= and does not interpret a script argument as -e.
  probe <- file.path(d, "script with spaces.R")
  put(probe, c(
    'source(Sys.getenv("POE_R_BOOTSTRAP"))',
    'f <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1L])',
    'f <- gsub("~+~", " ", f, fixed=TRUE)',
    paste0("stopifnot(identical(normalizePath(f), ", deparse(probe), "),"),
    paste0("identical(commandArgs(TRUE), ", paste(deparse(args), collapse = ""), "),"),
    'identical(.libPaths()[1L], normalizePath(Sys.getenv("POE_R_LIBRARY"))))',
    'cat("SCRIPT_ARGS_OK\\n")'
  ))
  z <- run("bash", c(wrapper, "--vanilla", probe, args))
  if (z$rc != 0L) stop(paste(z$out, collapse = "\n"))
  check(z$rc == 0L && "SCRIPT_ARGS_OK" %in% z$out)
  put(probe, c('source(Sys.getenv("POE_R_BOOTSTRAP"))',
               "quit(status=23L)"))
  check(run("bash", c(wrapper, probe))$rc == 23L)
  # Direct Rscript runs do not go through the host wrapper.
  Sys.unsetenv(c("LD_LIBRARY_PATH", "PYTHONPATH", "PYTHONHOME"))
  Sys.setenv(POE_R_LIBRARY = lib)
  direct <- file.path(d, "direct probe.R")
  put(direct, c(
    'before <- .libPaths()',
    paste0("source(", deparse(file.path(root, "R/library_paths.R")), ")"),
    'stopifnot(identical(.libPaths()[1L], normalizePath(Sys.getenv("POE_R_LIBRARY"))),',
    'all(before %in% .libPaths()))',
    'cat("DIRECT_OK\\n")'
  ))
  z <- run(rs, c("--vanilla", direct))
  check(z$rc == 0L && "DIRECT_OK" %in% z$out)
  # Unset environment selects the explicit config default, not the cwd.
  Sys.unsetenv("POE_R_LIBRARY")
  default_probe <- file.path(d, "default probe.R")
  put(default_probe, c(
    paste0("root <- ", deparse(fixture_root)),
    'source(file.path(root, "R/library_paths.R"))',
    paste0("stopifnot(identical(.libPaths()[1L], normalizePath(", deparse(lib), ")))"),
    'cat("CONFIG_DEFAULT_OK\\n")'
  ))
  z <- run(rs, c("--vanilla", default_probe))
  check(z$rc == 0L && "CONFIG_DEFAULT_OK" %in% z$out)
  # Verify initialization before actual YAML lookup, without relying on R_LIBS_USER.
  yaml_lib <- dirname(find.package("yaml"))
  Sys.setenv(POE_R_LIBRARY = yaml_lib)
  yaml_input <- put(file.path(d, "input.yaml"), "startup: verified")
  put(direct, c(
    ".libPaths(.Library)",
    paste0("source(", deparse(file.path(root, "R/library_paths.R")), ")"),
    paste0("stopifnot(identical(.libPaths()[1L], normalizePath(",
           deparse(yaml_lib), ")),"),
    paste0('identical(yaml::read_yaml(', deparse(yaml_input), ')$startup, "verified"))'),
    'cat("YAML_AFTER_INIT_OK\\n")'
  ))
  z <- run(rs, c("--vanilla", direct))
  check(z$rc == 0L && "YAML_AFTER_INIT_OK" %in% z$out)
  entries <- c("poe.R", "run_pipeline.R", "runtime.R", "00_setup/preflight.R",
    "01_retrieval/retrieve_public.R", "02_reference/prepare_reference.R",
    "03_genotypes/prepare_child.R", "03_genotypes/prepare_markers.R",
    "05_parental_expression/assign_parents.R", "06_promoter_models/run_models.R",
    "07_exports/prepare_candidate.R", "07_exports/release_aggregates.R")
  for (entry in entries) {
    text <- readLines(file.path(root, entry))
    init <- grep('source(file.path(root,"R/library_paths.R"))', text, fixed = TRUE)
    config <- grep('source(file.path(root,"R/config.R"))', text, fixed = TRUE)
    check(length(init) == 1L && length(config) == 1L && init < config)
    Sys.setenv(POE_R_LIBRARY = lib)
    z <- run(rs, c("--vanilla", file.path(root, entry)))
    check(z$rc != 0L && any(grepl("Usage:|One explicit private CONFIG", z$out)))
    Sys.setenv(POE_R_LIBRARY = file.path(d, "missing library"))
    z <- run(rs, c("--vanilla", file.path(root, entry)))
    check(z$rc != 0L && any(grepl("R library is missing", z$out, fixed = TRUE)))
  }
  Sys.setenv(POE_R_LIBRARY = "")
  z <- run(rs, c("--vanilla", direct))
  check(z$rc != 0L && any(grepl("R library is missing", z$out, fixed = TRUE)))
  Sys.setenv(POE_R_LIBRARY = lib)
  Sys.chmod(lib, "0000")
  if (file.access(lib, 5L) != 0L) {
    z <- run(rs, c("--vanilla", direct))
    check(z$rc != 0L && any(grepl("not readable/searchable", z$out, fixed = TRUE)))
  } else {
    cat("SKIP unreadable-directory check: current user bypasses permission bits.\n")
  }
  Sys.chmod(lib, "0700")
  # The existing default is retained explicitly, not guessed from the cwd/YAML.
  check(any(grepl("/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1",
                 readLines(file.path(root, "config/r_library.path")), fixed = TRUE)))
  cat("PASS:", n, "explicit-library startup/argument assertions (mock container).\n")
}
main()
