#!/usr/bin/env Rscript
# Standalone entry-point startup with an isolated package search path.
a <- commandArgs(TRUE)
root <- normalizePath(if (length(a)) a[1] else ".")
source(file.path(root, "R/library_paths.R"))
checks <- character()
check <- function(label, value) {
  if (!isTRUE(value)) stop("FAIL ", label)
  checks <<- c(checks, label)
  cat("PASS ", label, "\n", sep = "")
}
entries <- list.files(root, "[.]R$", recursive = TRUE)
entries <- entries[!startsWith(entries, "R/")]
startup_rows <- lapply(entries, function(f) {
  lines <- readLines(file.path(root, f), warn = FALSE)
  init <- grep("^[[:space:]]*source\\(.*R/library_paths[.]R", lines)
  tokens <- getParseData(parse(text = lines, keep.source = TRUE))
  calls <- tokens[tokens$token == "SYMBOL_FUNCTION_CALL", ]
  loads <- calls$line1[calls$text %in% c("library", "require", "requireNamespace")]
  sources <- calls$line1[calls$text == "source"]
  first_load <- if (length(loads)) min(loads) else Inf
  first_source <- if (length(sources)) min(sources) else Inf
  data.frame(entry_point = f, initializes_library = length(init) == 1L,
    before_package_calls = length(init) == 1L && init < first_load,
    first_source_is_initializer = length(init) == 1L && init == first_source)
})
startup <- do.call(rbind, startup_rows)
check("all standalone R entry points initialize their library",
  nrow(startup) > 0L && all(startup$initializes_library))
check("initialization precedes package calls and module sourcing",
  all(startup$before_package_calls) && all(startup$first_source_is_initializer))

home <- tempfile("R startup fixture ")
dir.create(home)
fixture <- file.path(home, "source with spaces")
for (p in c("bin", "R", "config")) dir.create(file.path(fixture, p), recursive = TRUE)
for (p in c("bin/prepare_promoter_configs.R", "R/library_paths.R", "R/h1_annotations.R"))
  stopifnot(file.copy(file.path(root, p), file.path(fixture, p)))
selected <- file.path(home, "selected library")
hidden <- file.path(home, "empty default library")
dir.create(selected); dir.create(hidden)
stopifnot(file.copy(find.package("yaml"), selected, recursive = TRUE))
writeLines(selected, file.path(fixture, "config/r_library.path"))

profile <- file.path(home, "isolated_profile.R")
literal <- function(x) paste(deparse(x), collapse = "")
writeLines(c(
  "for (p in c('utils', 'stats', 'graphics', 'grDevices', 'datasets', 'methods'))",
  "  suppressPackageStartupMessages(library(p, character.only = TRUE))",
  paste0("assign('.lib.loc', ", literal(hidden), ", envir = environment(.libPaths))"),
  "stopifnot(!requireNamespace('yaml', quietly = TRUE))",
  "cat('ISOLATION: yaml unavailable before entry point\\n')",
  ".Last <- function() {",
  paste0("  stopifnot(identical(normalizePath(.libPaths()[1]), ", literal(normalizePath(selected)), "))"),
  paste0("  stopifnot(identical(normalizePath(find.package('yaml')), ", literal(normalizePath(file.path(selected, "yaml"))), "))"),
  "  cat('ISOLATION: configured library supplies yaml\\n')",
  "}"), profile)

project <- file.path(home, "project")
cohorts <- c("GUSTO", "mouse_Babr_CAST", "mouse_FVB_CAST_E12",
  "maize_B73_Ki11", "maize_B73_W22", "fish_Pprolifica_AB")
for (id in cohorts) {
  cfg <- yaml::read_yaml(file.path(root, "config", paste0(id, ".example.yaml")))
  cfg$baseline_nb <- list(output_dir = file.path(project, "data", id, "results", "baseline_nb_v1_run1"))
  path <- file.path(project, "data", id, "private", paste0(id, ".baseline_nb.yaml"))
  dir.create(dirname(path), recursive = TRUE)
  yaml::write_yaml(cfg, path)
}
operator <- file.path(home, "operator.yaml")
yaml::write_yaml(list(h1_bayes = list(count_provenance = "assignment_covariance_unknown",
  producing_pipeline = "software startup fixture declaration")), operator)

invoke <- function(script, args, library, command = file.path(R.home("bin"), "Rscript"),
                   wrapper = FALSE, extra_env = character()) {
  vars <- c("POE_R_LIBRARY", "R_PROFILE_USER", "R_LIBS_USER", "R_LIBS_SITE",
    "SCRIPTS_DIR", "POE_SITE_FILE")
  old <- Sys.getenv(vars, unset = NA_character_)
  on.exit(for (v in vars) {
    if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, setNames(list(old[[v]]), v))
  }, add = TRUE)
  Sys.setenv(R_PROFILE_USER = profile, R_LIBS_USER = hidden, R_LIBS_SITE = hidden)
  if (is.null(library)) Sys.unsetenv("POE_R_LIBRARY") else Sys.setenv(POE_R_LIBRARY = library)
  if (wrapper) {
    Sys.setenv(SCRIPTS_DIR = fixture, POE_SITE_FILE = file.path(home, "site.sh"))
    command_args <- c(script, args)
  } else command_args <- c("--no-environ", "--no-site-file", script, args)
  out <- suppressWarnings(system2(command, shQuote(command_args), env = extra_env,
    stdout = TRUE, stderr = TRUE))
  list(code = if (is.null(attr(out, "status"))) 0L else attr(out, "status"),
    text = paste(out, collapse = "\n"))
}
success <- function(x) x$code == 0L &&
  grepl("ISOLATION: yaml unavailable before entry point", x$text, fixed = TRUE) &&
  grepl("ISOLATION: configured library supplies yaml", x$text, fixed = TRUE)
generator <- file.path(fixture, "bin/prepare_promoter_configs.R")
negative <- file.path(home, "package_without_initializer.R")
writeLines("suppressPackageStartupMessages(library(yaml))", negative)
negative_result <- invoke(negative, character(), selected)
check("isolated dependency is unavailable without initialization",
  negative_result$code != 0L && grepl("there is no package called", negative_result$text, fixed = TRUE))

explicit_dest <- file.path(home, "explicit_configs")
explicit <- invoke(generator, c(project, operator, explicit_dest), selected)
if (!success(explicit)) stop(explicit$text)
check("generator uses explicit library with hidden default dependency", success(explicit))
default_dest <- file.path(home, "default_configs")
default <- invoke(generator, c(project, operator, default_dest), NULL)
if (!success(default)) stop(default$text)
check("generator uses configured default library", success(default))
manifest <- read.delim(file.path(explicit_dest, "CONFIGS.tsv"))
check("isolated startup still generates eighteen independent configurations",
  nrow(manifest) == 18L && length(unique(manifest$cohort_id)) == 6L &&
    setequal(manifest$radius, c(0L, 1000L, 3000L)) &&
    length(unique(manifest$destination_config)) == 18L)
check("explicit and configured-default generation agree",
  identical(read.delim(file.path(default_dest, "CONFIGS.tsv"))[c("cohort_id", "radius", "poe_input_config")],
    manifest[c("cohort_id", "radius", "poe_input_config")]))

bad_library <- file.path(home, "missing library")
invalid <- invoke(generator, c(project, operator, file.path(home, "invalid_configs")), bad_library)
check("invalid explicit library fails before package loading or writing",
  invalid$code != 0L && grepl("R library is missing or not readable/searchable", invalid$text, fixed = TRUE) &&
    !dir.exists(file.path(home, "invalid_configs")))
sampler <- invoke(file.path(root, "tests/test_h1_sampler.R"),
  c(root, file.path(home, "unused_sampler"), file.path(home, "absent_fixture.rds")), bad_library)
check("sampler entry point initializes before reading fixture or fitting",
  sampler$code != 0L && grepl("R library is missing or not readable/searchable", sampler$text, fixed = TRUE) &&
    !dir.exists(file.path(home, "unused_sampler")))

mock <- file.path(home, "container_mock")
writeLines(c("#!/bin/bash", "set -euo pipefail", "shift",
  'while [[ $# -gt 0 ]]; do',
  '  case "$1" in',
  '    --cleanenv) shift ;;',
  '    --bind) shift 2 ;;',
  '    --env) export "$2"; shift 2 ;;',
  '    *) shift; break ;;',
  '  esac',
  'done',
  'exec "$@"'), mock)
Sys.chmod(mock, "0755")
image <- file.path(home, "fixture.sif")
file.create(image)
writeLines(c(paste0("POE_R_SIF=", shQuote(image)),
  paste0("POE_R_LIBRARY=", shQuote(selected)),
  paste0("POE_SINGULARITY=", shQuote(mock))), file.path(home, "site.sh"))
wrapped <- invoke(generator, c(project, operator, file.path(home, "wrapped_configs")),
  selected, command = file.path(root, "bin/Rscript_sif"), wrapper = TRUE)
if (!success(wrapped)) stop(wrapped$text)
check("container wrapper mock executes generator with configured library", success(wrapped))

if (length(a) >= 2L) {
  dir.create(a[2], recursive = TRUE, showWarnings = FALSE)
  write.table(startup, file.path(a[2], "R_entry_point_startup.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(check = checks, status = "pass"),
    file.path(a[2], "R_library_startup_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
unlink(home, recursive = TRUE)
cat("R library startup checks: ", length(checks), "; entry points: ", nrow(startup), "\n", sep = "")
