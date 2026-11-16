#!/usr/bin/env Rscript
# Repository text and archive dispatch contracts; no fitting.
a <- commandArgs(TRUE)
root <- normalizePath(if (length(a)) a[1] else ".")
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages(library(data.table))
for (m in c("common", "config", "regulatory_common", "regulatory_resources",
            "regulatory_fish", "baseline_nb", "baseline_handoff"))
  source(file.path(root, "R", paste0(m, ".R")))
checks <- character()
check <- function(label, expr) {
  if (!isTRUE(expr)) stop("FAIL ", label)
  checks <<- c(checks, label); cat("PASS ", label, "\n", sep = "")
}
home <- tempfile("repository-contract-"); dir.create(home)
files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
files <- files[!grepl("/[.]git/|/__pycache__/", files) & !dir.exists(files)]
files <- files[grepl("[.](R|r|py|sh|md|yaml|yml|tsv|csv|txt)$", files)]
# Construct tokens from parts so this guard does not match its own source.
tokens <- c(paste0("check", "sum"), paste0("sha", "[-_]?(1|224|256|384|512)"),
  paste0("m", "d5"), paste0("di", "gest"), paste0("ha", "sh"), paste0("finger", "print"),
  paste0("ck", "sum"), paste0("cr", "c32"))
hits <- files[vapply(files, function(p)
  any(grepl(paste(tokens, collapse = "|"), readLines(p, warn = FALSE), ignore.case = TRUE)), logical(1))]
if (length(hits)) stop("Unexpected content identity operations/fields in: ",
  paste(substring(hits, nchar(root) + 2L), collapse = ", "))
check("maintained text uses descriptive identity only", !length(hits))
check("runtime provenance remains structured", is.list(nb_runtime_key()))
fasta <- file.path(home, "fixture.fa")
writeLines(c(">software_contig", "ACGTNN"), fasta)
inventory <- reg_fasta_inventory(fasta)
check("FASTA inventory reports names and lengths only",
  identical(names(inventory), c("chrom", "length")) && inventory$length == 6L)
z <- list(build = "GCF_027474105.1", output_dir = file.path(home, "regulatory"),
  cache_dir = file.path(home, "cache"))
dir.create(z$output_dir)
resources <- data.table(resource_id = c("fish_genome", "fish_report"),
  path = file.path(home, c("resource_A.txt", "resource_B.txt")), format = "text",
  build = z$build, url = NA_character_)
for (p in resources$path) writeLines("software contract fixture", p)
manifest <- file.path(home, "resources.tsv"); fwrite(resources, manifest, sep = "\t")
cfg <- list(species = "fish")
reg_finalize_resources(cfg, z, root, manifest)
rd <- fread(file.path(z$output_dir, "resource_manifest.tsv"))
check("regulatory resource receipt completes with descriptive fields",
  all(c("bytes", "retrieved_utc") %in% names(rd)) && reg_valid(z, "regulatory-resources"))
local <- file.path(home, "local.txt"); writeLines("local format fixture", local)
status <- system2("bash", c(shQuote(file.path(root, "bin/fetch_regulatory.sh")),
  shQuote("https://ftp.ncbi.nlm.nih.gov/not-requested"), shQuote(local), "text"))
check("local regulatory transfer uses three arguments without network", status == 0L)
invoke <- function(args) {
  text <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(file.path(root, "bin/export_results.R")), shQuote(args)), stdout = TRUE, stderr = TRUE))
  list(status = if (is.null(attr(text, "status"))) 0L else attr(text, "status"), text = text)
}
cfg <- yaml::read_yaml(file.path(root, "config/maize_B73_Ki11.example.yaml"))
oldbase <- dirname(cfg$paths$work)
cfg$paths <- lapply(cfg$paths, function(p) sub(oldbase, home, p, fixed = TRUE))
cfg$templates <- lapply(cfg$templates, function(p) sub(oldbase, home, p, fixed = TRUE))
for (k in c("sources", "libraries", "units", "marker_sources")) {
  dir.create(dirname(cfg$paths[[k]]), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file.path(root, "manifests", basename(cfg$paths[[k]])), cfg$paths[[k]]))
}
cfg$baseline_nb <- list(output_dir = file.path(home, "baseline"))
config <- file.path(home, "cohort.yaml"); yaml::write_yaml(cfg, config)
parsed <- read_config_v2(config); s <- nb_settings(parsed)
for (f in NB_POE_LAYER_FILES) nb_write(data.table(note = "software_contract_only"), nb_path(s, f))
gene <- data.table(feature_id = character(), status = character(), beta_adjusted = numeric(),
  se = numeric(), p = numeric(), not_strictly_POE = logical())
for (f in c("gene_POE.tsv", "transcript_POE.tsv")) nb_write(gene, nb_path(s, f))
nb_finish(s, "prepare", list(software = TRUE), nb_path(s, NB_POE_LAYER_FILES))
nb_finish(s, "combine", list(software = TRUE), nb_path(s, c("gene_POE.tsv", "transcript_POE.tsv")))
r <- invoke(c("--baseline-nb", "--partial", cfg$cohort_id, "--out", file.path(home, "baseline.tar.gz"), config))
if (r$status != 0L) stop(paste(r$text, collapse = "\n"))
check("baseline CLI partial archive remains independently runnable", r$status == 0L)
r <- invoke(c("--force", "--out", file.path(home, "promoter.tar.gz"), config))
if (r$status != 0L) stop(paste(r$text, collapse = "\n"))
check("promoter CLI archive remains runnable", r$status == 0L)
for (name in c("baseline", "promoter")) {
  dest <- file.path(home, paste0("unpack_", name)); dir.create(dest)
  utils::untar(file.path(home, paste0(name, ".tar.gz")), exdir = dest)
  m <- fread(file.path(dest, name, "MANIFEST.tsv"))
  if ("state" %in% names(m)) m <- m[state == "included"]
  check(paste(name, "archive inventory sizes match"),
    all(file.info(file.path(dest, name, m$archive_path))$size == m$bytes))
}
if (length(a) > 1L) {
  dir.create(a[2], recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(check = checks, status = "pass"), file.path(a[2], "repository_contract_tests.tsv"), sep = "\t")
}
unlink(home, recursive = TRUE)
cat("Completed ", length(checks), " repository contract checks.\n", sep = "")
