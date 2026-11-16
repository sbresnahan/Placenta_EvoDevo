#!/usr/bin/env Rscript
a <- commandArgs(TRUE)
root <- normalizePath(if (length(a)) a[1] else ".")
source(file.path(root, "R/library_paths.R"))
source(file.path(root, "R/cmdstan_runtime.R"))
main <- function() {
  checks <- character()
  check <- function(label, value) {
    if (!isTRUE(value)) stop("FAIL ", label)
    checks <<- c(checks, label); cat("PASS ", label, "\n", sep = "")
  }
  fails <- function(expr, text) {
    e <- tryCatch({force(expr); NULL}, error = identity)
    ok <- inherits(e, "error") && grepl(text, conditionMessage(e), fixed = TRUE)
    if (!ok) message("Expected error containing: ", text, "\nActual: ",
      if (inherits(e, "error")) conditionMessage(e) else "no error")
    ok
  }
  vars <- c("POE_CMDSTAN_ENV", "CONDA_PREFIX", "CMDSTAN", "MAKE", "CXX", "PATH",
    "LD_LIBRARY_PATH", "PYTHONPATH", "PYTHONHOME", "SCRIPTS_DIR", "POE_SITE_FILE")
  old <- Sys.getenv(vars, unset = NA_character_)
  on.exit(for (v in vars) {
    if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, setNames(list(old[[v]]), v))
  }, add = TRUE)
  Sys.unsetenv(c("POE_CMDSTAN_ENV", "CMDSTAN"))
  home <- tempfile("cmdstan-runtime-"); dir.create(home)
  on.exit(unlink(home, recursive = TRUE), add = TRUE)
  prefix <- file.path(home, "tools")
  stan <- file.path(prefix, "bin/cmdstan")
  dir.create(file.path(stan, "bin"), recursive = TRUE)
  dir.create(file.path(stan, "make"))
  writeLines("# software fixture", file.path(stan, "makefile"))
  executable <- function(path, body = "exit 0") {
    writeLines(c("#!/bin/bash", body), path); Sys.chmod(path, "0755")
  }
  executable(file.path(stan, "bin/stanc"))
  compiler <- file.path(prefix, "bin/x86_64-conda-linux-gnu-c++")
  executable(compiler, "printf 'software-test-compiler\\n'")
  stopifnot(file.copy(Sys.which("make"), file.path(prefix, "bin/make")))
  Sys.chmod(file.path(prefix, "bin/make"), "0755")
  configroot <- file.path(home, "source")
  dir.create(file.path(configroot, "config"), recursive = TRUE)
  writeLines(prefix, file.path(configroot, "config/cmdstan_env.path"))
  preserved <- c(paste0("TBB_LIB=", prefix, "/lib/"),
    paste0("TBB_INC=", prefix, "/include/"), "LDFLAGS+=-pthread", "PRECOMPILED_HEADERS=false")
  local <- file.path(stan, "make/local")
  writeLines(preserved, local)
  calls <- new.env(); calls$set <- calls$rebuild <- calls$check <- calls$model <- 0L
  api <- list(
    set_path = function(p) {stopifnot(identical(p, stan)); calls$set <- calls$set + 1L},
    rebuild = function() calls$rebuild <- calls$rebuild + 1L,
    check = function() calls$check <- calls$check + 1L,
    model = function(path) {
      stopifnot(file.exists(path)); calls$model <- calls$model + 1L
      exe <- sub("[.]stan$", "", path); executable(exe)
      list(exe_file = function() exe)
    })
  before_library <- .libPaths()
  before_ld <- Sys.getenv("LD_LIBRARY_PATH", unset = NA_character_)
  check("missing prefix fails without installing tools",
    fails(poe_cmdstan_runtime(root, prefix = file.path(home, "missing"), .api = api),
      "Missing CmdStan tool environment"))
  check("makefile CXX is checked instead of trusting the R environment",
    fails(poe_cmdstan_runtime(root, prefix = prefix, .api = api),
      "make/local does not select"))
  check("read-only failure preserves make/local and skips rebuild",
    identical(readLines(local), preserved) && calls$rebuild == 0L)
  Sys.unsetenv("CMDSTAN")
  result <- poe_cmdstan_runtime(configroot, "configure", prefix = "", .api = api)
  check("configured-default prefix resolves tools and registers CmdStan",
    identical(result$prefix, prefix) && identical(result$compiler, compiler) &&
      identical(Sys.getenv("MAKE"), file.path(prefix, "bin/make")) && calls$set == 1L)
  check("explicit setup preserves TBB/linker settings and rebuilds once",
    identical(readLines(local)[seq_along(preserved)], preserved) && calls$rebuild == 1L &&
      length(list.files(dirname(local), "^local[.]poe-backup-")) == 1L)
  check("R libraries and linker environment remain unchanged",
    identical(before_library, .libPaths()) &&
      identical(before_ld, Sys.getenv("LD_LIBRARY_PATH", unset = NA_character_)))
  snapshot <- readBin(local, "raw", n = file.info(local)$size)
  checked <- poe_cmdstan_runtime(root, prefix = prefix, .api = api)
  check("fit-time check does not edit or rebuild the shared installation",
    identical(snapshot, readBin(local, "raw", n = file.info(local)$size)) &&
      calls$rebuild == 1L && calls$model == 0L && calls$check == 0L)
  check("tool bin takes PATH precedence without duplicate entries",
    identical(strsplit(Sys.getenv("PATH"), ":", fixed = TRUE)[[1]][1], file.path(prefix, "bin")) &&
      sum(strsplit(Sys.getenv("PATH"), ":", fixed = TRUE)[[1]] == file.path(prefix, "bin")) == 1L)
  poe_cmdstan_runtime(root, "compile", prefix = prefix, .api = api)
  check("compile preflight invokes checker and model compilation without sampling",
    calls$check == 1L && calls$model == 1L && calls$rebuild == 1L)
  Sys.setenv(CMDSTAN = home)
  check("conflicting CmdStan installation is rejected",
    fails(poe_cmdstan_runtime(root, prefix = prefix, .api = api), "CMDSTAN conflicts"))
  Sys.unsetenv("CMDSTAN")
  Sys.chmod(compiler, "0644")
  check("non-executable compiler is rejected",
    fails(poe_cmdstan_runtime(root, prefix = prefix, .api = api), "one Conda Linux C++ compiler"))
  Sys.chmod(compiler, "0755")
  check("tool prefix rejects make and shell metacharacters",
    fails(poe_cmdstan_runtime(root, prefix = "/bad tool", .api = api), "absolute tool prefix"))
  # Exercise direct fitting before any brms invocation; no model is fitted.
  e <- new.env(parent = globalenv()); e$root <- root
  sys.source(file.path(root, "R/h1_bayes.R"), envir = e)
  e$poe_h1_backend_check <- function(...) TRUE
  e$h1_encode <- function(d) list(data = d, map = "map")
  e$h1_cmdstan_runtime <- function(...) stop("RUNTIME_SENTINEL")
  check("direct h1_fit_once initializes before calling brms",
    fails(e$h1_fit_once(NULL, NULL, list(), list(), 1L, file.path(home, "absent.rds")),
      "RUNTIME_SENTINEL"))
  suppressPackageStartupMessages(library(data.table))
  saved <- file.path(home, "saved.rds")
  saveRDS(list(model_spec = e$H1_MODEL_SPEC, settings = list(), profile = list(),
    seed = 1L, map = "map", derived = NULL,
    diagnostics = list(status = "ready", summary = "summary"), fit = NULL), saved)
  e$h1_summarize <- function(x) data.table(value = 1)
  e$h1_residuals <- function(...) "residuals"
  resumed <- e$h1_fit_once(NULL, NULL, list(), list(), 1L, saved)
  check("compatible saved task resumes without toolchain initialization",
    identical(resumed$residuals, "residuals") && resumed$summary$status == "ready")
  gene_calls <- all.names(body(e$h1_gene), functions = TRUE)
  check("gene-level path validates runtime before task loop",
    match("h1_cmdstan_runtime", gene_calls) < match("h1_profiles", gene_calls))
  # Actual wrapper, with a clean-environment executable double for the container.
  mock <- file.path(home, "container")
  executable(mock, c(
    '[[ -z "${LD_LIBRARY_PATH+x}" && -z "${PYTHONPATH+x}" && -z "${PYTHONHOME+x}" ]] || exit 81',
    '[[ "$1" == exec ]] || exit 82', 'shift', 'forward=()', 'clean=0',
    'while [[ "$1" == --* ]]; do',
    '  case "$1" in',
    '    --cleanenv) clean=1; shift ;;',
    '    --bind) shift 2 ;;',
    '    --env) forward+=("$2"); shift 2 ;;',
    '    *) exit 83 ;;',
    '  esac',
    'done', '[[ "$clean" == 1 ]] || exit 84', 'shift',
    '[[ "$1" == Rscript ]] || exit 85', 'shift',
    paste0('exec /usr/bin/env -i PATH=/usr/bin:/bin HOME="$HOME" "${forward[@]}" ',
      shQuote(file.path(R.home("bin"), "Rscript")), ' "$@"')))
  image <- file.path(home, "R.sif"); file.create(image)
  site <- file.path(home, "site.sh")
  writeLines(c(paste0("POE_R_SIF=", shQuote(image)),
    paste0("POE_R_LIBRARY=", shQuote(.libPaths()[1])),
    paste0("POE_SINGULARITY=", shQuote(mock))), site)
  Sys.setenv(SCRIPTS_DIR = root, POE_SITE_FILE = site, POE_CMDSTAN_ENV = prefix,
    CMDSTAN = stan, CONDA_PREFIX = "/host/unrelated", LD_LIBRARY_PATH = "/host/lib",
    PYTHONPATH = "/host/python", PYTHONHOME = "/host/pythonhome")
  expression <- paste0(
    "stopifnot(Sys.getenv('POE_CMDSTAN_ENV') == ", deparse(prefix), ", ",
    "Sys.getenv('CMDSTAN') == ", deparse(stan), ", ",
    "Sys.getenv('CONDA_PREFIX') == '', Sys.getenv('LD_LIBRARY_PATH') == '', ",
    "Sys.getenv('PYTHONPATH') == '', Sys.getenv('PYTHONHOME') == '', ",
    "normalizePath(.libPaths()[1]) == ", deparse(normalizePath(.libPaths()[1])), "); ",
    "cat('WRAPPER_OK\\n')")
  output <- suppressWarnings(system2(file.path(root, "bin/Rscript_sif"),
    c("-e", shQuote(expression)), stdout = TRUE, stderr = TRUE))
  check("wrapper forwards only selected tool settings into a clean R session",
    is.null(attr(output, "status")) && "WRAPPER_OK" %in% output)
  if (!is.null(attr(output, "status"))) stop(paste(output, collapse = "\n"))
  if (length(a) >= 2L) {
    dir.create(a[2], recursive = TRUE, showWarnings = FALSE)
    write.table(data.frame(check = checks, status = "pass", evidence = "local software test doubles"),
      file.path(a[2], "cmdstan_runtime_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  cat("CmdStan runtime checks: ", length(checks), "; no real compilation or sampling\n", sep = "")
}
main()
