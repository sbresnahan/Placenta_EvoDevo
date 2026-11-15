#!/usr/bin/env Rscript
# Software-only CLI and metadata fixtures. No site jobs are submitted.
a <- commandArgs(TRUE)
root <- if (length(a)) normalizePath(a[1]) else normalizePath(".")
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages(library(data.table))
source(file.path(root, "R/config.R"))
source(file.path(root, "R/common.R"))
source(file.path(root, "R/regulatory_common.R"))
source(file.path(root, "R/baseline_nb.R"))
checks <- 0L
check <- function(name, expr) {stopifnot(isTRUE(expr)); checks <<- checks + 1L; cat("PASS ", name, "\n", sep = "")}
fails <- function(expr) inherits(tryCatch(force(expr), error = function(e) e), "error")
home <- tempfile("nb-deployment-"); dir.create(home)
cfg <- yaml::read_yaml(file.path(root, "config/GUSTO.example.yaml"))
cfg$gusto_precomputed$fastq_id_map <- file.path(home, "link.tsv")
cfg$paths$libraries <- file.path(home, "libraries.tsv")
cfg$baseline_nb <- list(metadata = file.path(home, "metadata.tsv"))
fwrite(data.table(subject = c("fixture_infantA", "fixture_infantB"), fastq = c("001", "002")),
  cfg$gusto_precomputed$fastq_id_map, sep = "\t", col.names = FALSE)
fwrite(data.table(sample_id = c("J002", "J001"), unit_id = c("unitB", "unitA")),
  cfg$paths$libraries, sep = "\t")
cov <- data.table(SubjectID = c("fixture_infantA", "fixture_infantB"), sex = c("F", "M"),
  mother_ethnicity = c("fixture_A", "fixture_B"))
csv <- file.path(home, "cov.csv"); fwrite(cov, csv)
configfile <- file.path(home, "config.yaml"); yaml::write_yaml(cfg, configfile)
parsed <- read_config_v2(configfile)
check("full schema config accepts baseline extension", parsed$baseline_nb$metadata == cfg$baseline_nb$metadata)
old <- Sys.getenv("POE_R_LIBRARY", unset = NA_character_)
Sys.setenv(POE_R_LIBRARY = .libPaths()[1])
invoke <- function(args) {
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(file.path(root, "bin/prepare_baseline_gusto_metadata.R")), shQuote(args)),
    stdout = TRUE, stderr = TRUE))
  list(code = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), text = out)
}
output <- cfg$baseline_nb$metadata
r <- invoke(c(configfile, csv, output))
if (r$code != 0L) stop(paste(r$text, collapse = "\n"))
d <- fread(output, colClasses = "character")
check("metadata CLI exact J FASTQ linkage", d[sample_id == "J001"]$ancestry == "fixture_A" &&
  d[sample_id == "J002"]$sex == "M")
check("metadata CLI refuses overwrite", invoke(c(configfile, csv, output))$code != 0L)
fwrite(cov[1], csv)
check("metadata CLI rejects missing infant", invoke(c(configfile, csv, paste0(output, ".missing")))$code != 0L)
fwrite(rbind(cov, cov[1]), csv)
check("metadata CLI rejects duplicate infant", invoke(c(configfile, csv, paste0(output, ".dup")))$code != 0L)
fwrite(data.table(sample_id = c("J001", "J002")), output, sep = "\t")
human_linkage <- function(cfg, m) data.table(sample_id = c("J001", "J002"), family_id = c("famA", "famB"))
m <- list(libraries = fread(cfg$paths$libraries))
s <- nb_settings(parsed)
check("human model rejects missing sex and ancestry", fails(nb_metadata(parsed, s, m)))
fwrite(data.table(sample_id = c("J001", "J002"), sex = "F", ancestry = "fixture_A"), output, sep = "\t")
check("human model rejects unidentifiable covariates", fails(nb_metadata(parsed, s, m)))
if (is.na(old)) Sys.unsetenv("POE_R_LIBRARY") else Sys.setenv(POE_R_LIBRARY = old)

# Exercise the real submission helper with mocked site wrappers.
mock <- file.path(home, "mock_site"); dir.create(file.path(mock, "bin"), recursive = TRUE)
writeLines(c("#!/bin/bash", 'printf "%s\\n" "$@" > "$NB_TEST_CAPTURE"'),
  file.path(mock, "bin/bsub_site"))
writeLines(c("#!/bin/bash", "exit 0"), file.path(mock, "bin/Rscript_sif"))
file.copy(file.path(root, "bin/submit_baseline_nb.sh"), file.path(mock, "bin/submit_baseline_nb.sh"))
Sys.chmod(list.files(file.path(mock, "bin"), full.names = TRUE), "0755")
capture <- file.path(home, "scheduler_args.txt")
env <- c(paste0("SCRIPTS_DIR=", shQuote(mock)),
  paste0("NB_LOG_DIR=", shQuote(file.path(home, "logs"))),
  paste0("NB_TEST_CAPTURE=", shQuote(capture)), "NB_ARRAY_END=3", "NB_MAX_CONCURRENT=2")
out <- system2("bash", c(shQuote(file.path(mock, "bin/submit_baseline_nb.sh")),
  shQuote(configfile), "fit", "transcript", "ARRAY_INDEX"), env = env, stdout = TRUE, stderr = TRUE)
check("submission helper returns mocked scheduler success", is.null(attr(out, "status")))
args <- readLines(capture)
check("submission helper uses bounded job array", "nb_fit[1-3]%2" %in% args)
script <- tail(args, 1)
check("worker script preserves scheduler index expansion",
  any(grepl('${LSB_JOBINDEX:?Missing array index}', readLines(script), fixed = TRUE)))
check("worker script uses container R entry point",
  any(grepl("Rscript_sif", readLines(script), fixed = TRUE)))
cat("Completed ", checks, " deployment software checks. No site jobs submitted.\n", sep = "")
unlink(home, recursive = TRUE)
