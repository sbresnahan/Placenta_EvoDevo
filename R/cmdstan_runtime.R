# Tool-only Conda integration. Loading this module has no runtime side effects.
# Call after R/library_paths.R; never activate Conda or change R libraries here.
poe_cmdstan_runtime <- function(root, mode = c("check", "configure", "compile"),
                                prefix = Sys.getenv("POE_CMDSTAN_ENV"), .api = NULL) {
  mode <- match.arg(mode)
  if (!nzchar(prefix)) {
    config <- file.path(root, "config/cmdstan_env.path")
    if (!file.exists(config)) stop("Missing CmdStan environment configuration: ", config)
    prefix <- trimws(readLines(config, warn = FALSE))
  }
  if (length(prefix) != 1L || !nzchar(prefix) ||
      !grepl("^/[A-Za-z0-9_./+-]+$", prefix))
    stop("POE_CMDSTAN_ENV must be one absolute tool prefix without whitespace or shell/make metacharacters")
  if (!dir.exists(prefix)) stop("Missing CmdStan tool environment: ", prefix,
    "; install tools and run bin/check_cmdstan_runtime.R configure")
  prefix <- normalizePath(prefix, mustWork = TRUE)
  bin <- file.path(prefix, "bin")
  stan <- file.path(bin, "cmdstan")
  if (!file.exists(file.path(stan, "makefile")) ||
      file.access(file.path(stan, "bin/stanc"), 1L) != 0L)
    stop("Expected a Conda CmdStan installation at ", stan)
  stan <- normalizePath(stan, mustWork = TRUE)
  override <- Sys.getenv("CMDSTAN")
  if (nzchar(override) && (!dir.exists(override) ||
      !identical(normalizePath(override, mustWork = TRUE), stan)))
    stop("CMDSTAN conflicts with POE_CMDSTAN_ENV; unset CMDSTAN or select the matching tool prefix")
  make <- file.path(bin, "make")
  candidates <- list.files(bin, "-conda-linux-gnu-c[+][+]$", full.names = TRUE)
  if (!length(candidates))
    candidates <- list.files(bin, "-conda-linux-gnu-g[+][+]$", full.names = TRUE)
  if (file.access(make, 1L) != 0L || length(candidates) != 1L ||
      any(file.access(candidates, 1L) != 0L))
    stop("Tool prefix needs executable make and one Conda Linux C++ compiler: ", bin)
  compiler <- candidates[1L]
  # Keep R and its package library unchanged. In particular, no LD_LIBRARY_PATH.
  path <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1L]]
  Sys.setenv(POE_CMDSTAN_ENV = prefix, CONDA_PREFIX = prefix, CMDSTAN = stan,
    MAKE = make, CXX = compiler,
    PATH = paste(c(bin, path[path != bin]), collapse = .Platform$path.sep))
  run <- function(exe, args) {
    out <- suppressWarnings(system2(exe, shQuote(args), stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0L)
      stop("CmdStan tool command failed: ", exe, "\n", paste(out, collapse = "\n"))
    out
  }
  compiler_version <- run(compiler, "--version")
  local <- file.path(stan, "make/local")
  if (mode == "configure") {
    before <- if (file.exists(local)) readLines(local, warn = FALSE) else character()
    begin <- "# BEGIN POE TOOL-ONLY COMPILER"
    end <- "# END POE TOOL-ONLY COMPILER"
    first <- which(before == begin); last <- which(before == end)
    if (length(first) || length(last)) {
      if (length(first) != 1L || length(last) != 1L ||
          last != length(before) || first >= last)
        stop("Malformed or non-final POE compiler block in make/local; inspect manually")
      before <- if (first > 1L) before[seq_len(first - 1L)] else character()
    }
    if (file.exists(local)) {
      backup <- paste0(local, ".poe-backup-", format(Sys.time(), "%Y%m%d%H%M%S"),
        "-", Sys.getpid())
      if (!file.copy(local, backup, overwrite = FALSE)) stop("Cannot back up make/local")
    }
    tmp <- tempfile("poe-local-", tmpdir = dirname(local))
    writeLines(c(before, begin, paste("CXX :=", compiler), end), tmp)
    if (file.exists(local)) Sys.chmod(tmp, file.info(local)$mode)
    if (!file.rename(tmp, local)) stop("Cannot install compiler configuration")
  }
  # Query make's effective value, not Sys.getenv("CXX"): makefile assignments can
  # override environment variables. The probe neither builds nor edits CmdStan.
  if (!file.exists(local)) stop("Missing make/local; run check_cmdstan_runtime.R configure once")
  probe <- tempfile(fileext = ".mk")
  on.exit(unlink(probe), add = TRUE)
  writeLines(c("CXX :=", paste("include", local), ".PHONY: poe-print-cxx",
    "poe-print-cxx:", "\t@printf '%s\\n' 'POE_CXX=$(CXX)'"), probe)
  effective <- run(make, c("--no-print-directory", "-s", "-f", probe, "poe-print-cxx"))
  effective <- sub("^POE_CXX=", "", effective[startsWith(effective, "POE_CXX=")])
  if (length(effective) != 1L || !identical(effective, compiler))
    stop("make/local does not select the tool environment C++ compiler; run configure once")
  if (is.null(.api)) {
    if (!requireNamespace("cmdstanr", quietly = TRUE))
      stop("cmdstanr is absent from the selected R library; do not install Conda R")
    .api <- list(set_path = cmdstanr::set_cmdstan_path,
      check = function() cmdstanr::check_cmdstan_toolchain(fix = FALSE, quiet = FALSE),
      rebuild = function() cmdstanr::rebuild_cmdstan(cores = 1, quiet = FALSE),
      model = function(path) cmdstanr::cmdstan_model(path, force_recompile = TRUE,
        quiet = FALSE))
  }
  .api$set_path(stan)
  if (mode == "configure") {
    # Explicit operator action only; never called by a fitting task.
    .api$rebuild()
  }
  if (mode == "compile") {
    .api$check()
    smoke <- tempfile("poe-cmdstan-")
    dir.create(smoke, mode = "0700")
    on.exit(unlink(smoke, recursive = TRUE), add = TRUE)
    model <- file.path(smoke, "toolchain_smoke.stan")
    writeLines(c("parameters { real y; }", "model { y ~ normal(0, 1); }"), model)
    compiled <- .api$model(model)
    if (!is.function(compiled$exe_file) ||
        file.access(compiled$exe_file(), 1L) != 0L)
      stop("Compile preflight did not produce an executable")
    # --help loads dynamic libraries without sampling.
    run(compiled$exe_file(), "help")
  }
  invisible(list(prefix = prefix, cmdstan = stan, make = make, compiler = compiler,
    compiler_version = compiler_version[1L], mode = mode))
}
