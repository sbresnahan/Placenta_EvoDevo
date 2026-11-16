#!/usr/bin/env Rscript
# Container and scheduler argument mocks only.
args <- commandArgs(TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
source(file.path(root, "R/library_paths.R"))
main <- function() {
  checks <- character()
  check <- function(label, x) {
    if (!isTRUE(x)) stop("FAIL ", label)
    checks <<- c(checks, label); cat("PASS ", label, "\n", sep = "")
  }
  vars <- c("H1_PREP_WORKERS", "H1_PREP_MEMORY_GB", "LSB_DJOB_NUMPROC",
    "LSB_MAX_NUM_PROCESSORS", "SCRIPTS_DIR", "POE_SITE_FILE", "NB_CPUS",
    "NB_MEMORY_GB", "NB_ARRAY_END", "NB_LOG_DIR", "H1_TEST_CAPTURE")
  old <- Sys.getenv(vars, unset = NA_character_)
  on.exit(for (i in seq_along(vars)) if (is.na(old[i])) Sys.unsetenv(vars[i]) else
    do.call(Sys.setenv, setNames(list(old[i]), vars[i])), add = TRUE)
  Sys.unsetenv(vars)
  home <- tempfile("h1 deployment "); dir.create(home)
  on.exit(unlink(home, recursive = TRUE), add = TRUE)
  bin <- file.path(home, "bin"); dir.create(bin)
  capture <- file.path(home, "arguments.txt")
  executable <- function(path, lines) {writeLines(c("#!/bin/bash", lines), path); Sys.chmod(path, "0755")}
  mock_container <- file.path(bin, "container")
  executable(mock_container, paste0("printf '%s\\n' \"$@\" > ", shQuote(capture)))
  sif <- file.path(home, "R.sif"); file.create(sif)
  site <- file.path(home, "site.sh")
  writeLines(c(paste0("POE_R_SIF=", shQuote(sif)),
    paste0("POE_R_LIBRARY=", shQuote(.libPaths()[1L])),
    paste0("POE_SINGULARITY=", shQuote(mock_container))), site)
  Sys.setenv(SCRIPTS_DIR = root, POE_SITE_FILE = site, H1_PREP_WORKERS = 2,
    H1_PREP_MEMORY_GB = 4, LSB_DJOB_NUMPROC = 2)
  invoke <- function(path, arguments) {
    out <- suppressWarnings(system2("bash", shQuote(c(path, arguments)), stdout = TRUE, stderr = TRUE))
    list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), text = out)
  }
  rc <- invoke(file.path(root, "bin/Rscript_sif"), c("poe_baseline.R", "cfg.yaml", "h1-inputs"))
  ca <- readLines(capture)
  check("clean container accepts bounded H1 controls", rc$status == 0L &&
    all(c("--cleanenv", "H1_PREP_WORKERS=2", "H1_PREP_MEMORY_GB=4",
      "LSB_DJOB_NUMPROC=2") %in% ca))
  check("preparation container sets one native thread per worker",
    all(paste0(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"), "=1") %in% ca))
  Sys.setenv(H1_PREP_WORKERS = "2;false")
  check("container rejects nonnumeric resource values",
    invoke(file.path(root, "bin/Rscript_sif"), "poe_baseline.R")$status != 0L)
  Sys.unsetenv(c("H1_PREP_WORKERS", "H1_PREP_MEMORY_GB", "LSB_DJOB_NUMPROC"))
  executable(file.path(bin, "bsub_site"), 'printf "%s\\n" "$@" > "$H1_TEST_CAPTURE"')
  worker_env <- file.path(home, "worker_controls.txt")
  executable(file.path(bin, "Rscript_sif"),
    paste0('printf "%s\\n" "$H1_PREP_WORKERS" "$H1_PREP_MEMORY_GB" > ', shQuote(worker_env)))
  file.copy(file.path(root, "bin/submit_baseline_nb.sh"), file.path(bin, "submit_baseline_nb.sh"))
  config <- file.path(home, "cfg.yaml"); writeLines("software: true", config)
  Sys.setenv(SCRIPTS_DIR = home, NB_CPUS = 4, NB_MEMORY_GB = 16,
    NB_LOG_DIR = file.path(home, "logs"), H1_TEST_CAPTURE = capture)
  submit <- function() invoke(file.path(bin, "submit_baseline_nb.sh"), c(config, "h1-inputs"))
  rc <- submit(); ca <- readLines(capture)
  check("preparation scheduler mock succeeds", rc$status == 0L)
  check("preparation processes are allocated on one host", "span[hosts=1]" %in% ca)
  check("preparation scheduler CPU and memory requests preserved",
    ca[match("-n", ca) + 1L] == "4" && ca[match("-M", ca) + 1L] == "16")
  script <- tail(ca, 1L)
  rc <- invoke(script, character())
  check("preparation worker defaults derive from requested resources",
    rc$status == 0L && identical(readLines(worker_env), c("4", "16")))
  Sys.setenv(H1_PREP_WORKERS = 2, H1_PREP_MEMORY_GB = 8)
  rc <- submit(); rc2 <- invoke(tail(readLines(capture), 1L), character())
  check("explicit smaller preparation resources reach worker",
    rc$status == 0L && rc2$status == 0L && identical(readLines(worker_env), c("2", "8")))
  Sys.setenv(H1_PREP_WORKERS = 5)
  check("worker count above requested CPUs is rejected", submit()$status != 0L)
  Sys.setenv(H1_PREP_WORKERS = 2, H1_PREP_MEMORY_GB = 17)
  check("memory budget above request is rejected", submit()$status != 0L)
  Sys.setenv(H1_PREP_MEMORY_GB = 0)
  check("zero memory budget is rejected", submit()$status != 0L)
  Sys.setenv(H1_PREP_MEMORY_GB = 16, NB_ARRAY_END = 2)
  check("preparation arrays are rejected to prevent duplicate writers", submit()$status != 0L)
  if (length(args) > 1L) {
    dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
    write.table(data.frame(check = checks, status = "pass"),
      file.path(args[2], "h1_prepare_deployment_checks.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE)
  }
  cat("Completed ", length(checks), " H1 deployment software checks. No site jobs submitted.\n", sep = "")
}
main()
