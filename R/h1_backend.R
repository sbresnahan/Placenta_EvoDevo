# brms uses RStan classes even when CmdStanR performs sampling.
h1_task_table_equal <- function(x, y) {
  if (!is.data.frame(x) || !is.data.frame(y)) return(identical(x, y))
  identical(names(x), names(y)) && identical(nrow(x), nrow(y)) &&
    all(vapply(names(x), function(n) identical(x[[n]], y[[n]]), TRUE))
}
poe_h1_backend_classes <- function() {
  dso <- methods::new("cxxdso", sig = list(character()), dso_saved = FALSE,
    dso_filename = character(), modulename = character(), system = R.version$system,
    cxxflags = character(), .CXXDSOMISC = new.env(parent = emptyenv()))
  model <- methods::new("stanmodel", model_name = "poe_backend_check",
    model_code = character(), model_cpp = list(), dso = dso)
  fit <- methods::new("stanfit", model_name = "poe_backend_check",
    model_pars = character(), par_dims = list(), mode = 0L, sim = list(),
    inits = list(), stan_args = list(), stanmodel = model, date = as.character(Sys.time()),
    .MISC = new.env(parent = emptyenv()))
  stopifnot(methods::is(fit, "stanfit"))
  invisible(TRUE)
}
poe_h1_backend_check <- function(.api = NULL) {
  if (is.null(.api)) .api <- list(
    load = function(p) loadNamespace(p),
    classes = poe_h1_backend_classes,
    version = function(p) as.character(utils::packageVersion(p)),
    path = function(p) find.package(p))
  packages <- c("rstan", "StanHeaders", "RcppParallel", "brms", "cmdstanr", "posterior")
  for (p in packages) tryCatch(.api$load(p), error = function(e)
    stop("H1 backend cannot load ", p, ": ", conditionMessage(e),
      "\nNo new sampling started. Run bin/repair_rstan_stack.R check through bin/Rscript_sif; ",
      "repair the selected R library before fitting. Fitting never installs packages.", call. = FALSE))
  tryCatch(.api$classes(), error = function(e)
    stop("H1 RStan class validation failed: ", conditionMessage(e),
      "\nNo new sampling started; validate the RStan stack before fitting.", call. = FALSE))
  data.frame(package = packages, version = vapply(packages, .api$version, ""),
    path = vapply(packages, .api$path, ""), stringsAsFactors = FALSE, row.names = NULL)
}
h1_archive_gene_failure <- function(s, job_index, paths) {
  dest <- nb_path(s, "private", "h1_recovery",
    sprintf("gene_%06d_%s_%d", job_index, format(Sys.time(), "%Y%m%dT%H%M%S"), Sys.getpid()))
  if (!dir.create(dest, recursive = TRUE, mode = "0700")) stop("Cannot create H1 recovery directory")
  for (p in paths[file.exists(paths)]) {
    to <- file.path(dest, basename(p))
    if (!file.copy(p, to, copy.mode = TRUE) ||
        !identical(readBin(p, "raw", n = file.info(p)$size),
          readBin(to, "raw", n = file.info(to)$size))) stop("Cannot preserve prior H1 output: ", p)
  }
  receipt <- nb_receipt(s, sprintf("h1_%06d", job_index))
  if (file.exists(receipt) && !file.rename(receipt, file.path(dest, basename(receipt))))
    stop("Cannot preserve prior failed H1 receipt")
  nb_save(list(job_index = job_index, action = "explicit_failed_task_retry",
    prior_files = basename(paths[file.exists(paths)])), file.path(dest, "recovery.rds"))
  message("Prior gene failure outputs preserved in ", dest)
  invisible(dest)
}
