#!/usr/bin/env Rscript
# End-to-end SOFTWARE harness, not biological POE data.
# Observed Salamanders example counts are reused with artificial labels.
# Raw acquisition, parent assignment and site inputs are mocked explicitly.
a <- commandArgs(trailingOnly = TRUE)
root <- if (length(a)) normalizePath(a[1]) else normalizePath(".")
if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("glmmTMB required")
suppressPackageStartupMessages({library(data.table); library(GenomicRanges)})
for (module in c("common", "config", "regulatory_common", "regulatory_fish",
    "regulatory_te_genes", "regulatory_context", "baseline_nb", "baseline_te", "baseline_handoff"))
  source(file.path(root, "R", paste0(module, ".R")))
checks <- 0L
check <- function(name, expr) {stopifnot(isTRUE(expr)); checks <<- checks + 1L; cat("PASS ", name, "\n", sep = "")}
fails <- function(expr) inherits(tryCatch(force(expr), error = function(e) e), "error")
home <- tempfile("baseline-workflow-test-"); dir.create(home)
data("Salamanders", package = "glmmTMB", envir = environment())
txmap <- data.table(transcript_id = c("fixture_txA", "fixture_txB"), gene_id = "fixture_gene",
  chrom = "chr1", strand = "+", start = c(10, 100), end = c(90, 180),
  tss = c(10, 100), first_donor = c(40, 140), first_start = c(10, 100),
  first_end = c(40, 140), n_exons = 2L)
paths <- setNames(file.path(home, paste0(c("libraries", "units", "tx_map", "event_map", "exons", "events"), ".tsv")),
  c("libraries", "units", "tx_map", "event_map", "exons", "events"))
metadata <- data.table(sample_id = sprintf("fixture%02d", 1:40),
  unit_id = sprintf("unit%02d", 1:40), family_id = sprintf("family%02d", 1:40),
  cluster_id = sprintf("family%02d", 1:40), independence_verified = TRUE,
  sex = rep(c("F", "M"), length.out = 40), ancestry = rep(c("A", "A", "B", "B"), length.out = 40),
  cross = NA_character_)
libraries <- metadata[, .(sample_id, unit_id, strand_mode = "none")]
units <- metadata[, .(unit_id, cluster_id, stage = "fixture_stage", tissue = "fixture_tissue")]
eventmap <- txmap[, .(transcript_id, gene_id, event_id = paste0("fixture_event", 1:2), window = 50L)]
for (nm in names(paths)) nb_write(switch(nm, libraries = libraries, units = units,
  tx_map = txmap, event_map = eventmap, exons = txmap, events = eventmap), paths[[nm]])
paths <- as.list(paths)
paths$model_gene <- file.path(home, "legacy", "model_gene.tsv")
paths$covariates <- file.path(home, "metadata.tsv")
nb_write(metadata[, .(sample_id, sex, ancestry)], paths$covariates)
paths$genome_fai <- file.path(home, "fixture.fa.fai")
writeLines("chr1\t1000\t0\t1000\t1001", paths$genome_fai)
paths$genome <- file.path(home, "fixture.fa")
writeLines(c(">chr1", paste(rep("ACGT", 250), collapse = "")), paths$genome)
paths$metrics <- file.path(home, "unused.tsv")
paths$scratch <- file.path(home, "scratch")
cfg <- list(profile = "human_trio", cohort_id = "SOFTWARE_FIXTURE_NOT_BIOLOGY", species = "human",
  paths = paths, reference = list(build = "software_fixture"),
  dosage = list(maternal = 1, paternal = 1), parents = list(focal = "M", alternate = "P"),
  regulatory = list(output_dir = file.path(home, "resources"), bootstrap = 200,
    te_genes = list(placebo_shifts = 1, permutations = 1)),
  environment_allowlist = character())
configfile <- file.path(home, "fixture.yaml")
yaml::write_yaml(list(software_fixture = TRUE), configfile); attr(cfg, "file") <- configfile
s <- list(output_dir = file.path(home, "baseline_nb_v1"), batches_gene = 2L, batches_transcript = 3L,
  min_snps = 3L, min_samples = 15L, min_total_counts = 50, min_per_cross = 2L,
  bootstrap = 2L, bootstrap_chunk = 1L, seed = 1147L, metadata = paths$covariates, max_dense_gb = 1)
input <- CJ(sample_id = metadata$sample_id, transcript_id = txmap$transcript_id, j = 1:4)
input[, i := seq_len(.N)]
input[, `:=`(gene_id = "fixture_gene",
  locus_id = ifelse(transcript_id == "fixture_txA", paste0("chr1:", j, ":A:G"),
    paste0("chr1:", j + 100, ":C:T")),
  maternal_count = Salamanders$count[((i - 1) %% nrow(Salamanders)) + 1L],
  paternal_count = Salamanders$count[((i + 159) %% nrow(Salamanders)) + 1L])]
rawfiles <- rbindlist(lapply(metadata$sample_id, function(sid) {
  dir <- file.path(home, "raw", sid); dir.create(dir, recursive = TRUE)
  for (nm in c("allele_counts", "thetas", "quant_receipt", "quant_log", "resource_log"))
    writeLines("software_fixture", file.path(dir, paste0(nm, ".txt")))
  writeLines(c("{", '"fixture_txA": 3,', '"fixture_txB": 1', "}"), file.path(dir, "tpms.txt"))
  data.table(sample_id = sid, input = c("allele_counts", "tpms", "quant_receipt"),
    path = file.path(dir, c("allele_counts.txt", "tpms.txt", "quant_receipt.txt")))
}))
rawfiles[, bytes := file.info(path)$size]
# Deliberate mocks isolate the downstream workflow from private upstream inputs.
nb_preflight <- function(cfg, s, root) list(manifests = list(libraries = libraries, units = units),
  metadata = copy(metadata), files = copy(rawfiles))
load_manifests <- function(cfg) list(libraries = libraries, units = units)
cfg_path <- function(c, key, sample, mode = "complete", strand = "none")
  file.path(home, "raw", sample, paste0(key, ".txt"))
stage_valid <- function(...) TRUE
assign_sample <- function(c, m, sample, return_transcripts = FALSE)
  copy(input[input$sample_id == sample, .(transcript_id, gene_id, locus_id, maternal_count, paternal_count)])
read_config_v2 <- function(path) cfg
nb_settings <- function(config) s
reg_promoter_gc <- function(genome, genes, contigs, flank) rep(.5, nrow(genes))
nb_prepare(cfg, s, root)
check("preparation receipt complete", !is.null(nb_require(s, "prepare")))
for (lv in c("gene", "transcript")) {
  manifest <- reg_read(nb_path(s, paste0(lv, "_batches.tsv")))
  for (b in manifest$batch) nb_fit_batch(cfg, s, root, lv, b)
}
nb_combine(cfg, s, root)
check("combined models have new version", all(reg_read(nb_path(s, "transcript_POE.tsv"))$version == NB_VERSION))
check("both baseline fixture transcripts fitted", nrow(nb_estimable(reg_read(nb_path(s, "transcript_POE.tsv")))) == 2)
dir.create(cfg$regulatory$output_dir, recursive = TRUE)
nb_write(data.table(chrom = "chr1", start = 8, end = 15, is_te = TRUE, is_unknown = FALSE,
  repeat_class = "DNA", repeat_family = "fixture_family", repeat_id = "fixture_repeat"),
  file.path(cfg$regulatory$output_dir, "repeat_intervals.tsv.gz"))
nb_write(data.table(resource_id = "repeats", scope = "full_reference"),
  file.path(cfg$regulatory$output_dir, "resource_manifest.tsv"))
nb_te_prepare(cfg, s, root)
check("disjoint-marker isoforms form a mixed-gene job", nrow(reg_read(nb_path(s, "within_gene_jobs.tsv"))) == 1)
nb_te_gene(cfg, s, root, 1)
nb_te_gene(cfg, s, root, 1)
nb_te_combine(cfg, s, root)
check("class and individual contrast exported", nrow(reg_read(nb_path(s, "within_gene_contrasts.tsv.gz"))) == 2)
nb_between_genes(cfg, s, root)
check("too few exposure genes reported explicitly",
  reg_read(nb_path(s, "between_genes", "te_genes_result.tsv"))$status == "insufficient_genes_per_group")
nb_context(cfg, s, root); nb_report_ready(cfg, s, root)
archive <- file.path(home, "software_fixture_results.tar.gz")
nb_export(configfile, archive, root)
listing <- utils::untar(archive, list = TRUE)
check("archive contains new gene model", any(endsWith(listing, "/gene_POE.tsv")))
check("private products excluded", !any(grepl("/private/|metadata.rds|legacy", listing)))
extract <- file.path(home, "unpacked"); dir.create(extract); utils::untar(archive, exdir = extract)
mf <- reg_read(list.files(extract, pattern = "^MANIFEST.tsv$", recursive = TRUE, full.names = TRUE)[1])
top <- file.path(extract, "software_fixture_results")
check("exported checksums match", all(vapply(file.path(top, mf$archive_path), reg_sha, "") == mf$sha256))
check("all stage statuses verify", all(nb_status(cfg, s, root)$status == "complete_verified"))
# A valid receipt cannot authorize a private identifier column for export.
report <- nb_path(s, "report_summary.tsv")
original <- reg_read(report); receipt <- nb_require(s, "report-ready")
bad <- copy(original); bad[, sample_id := "SOFTWARE_PRIVACY_FIXTURE"]
nb_write(bad, report); nb_finish(s, "report-ready", receipt$key, names(receipt$outputs))
check("privacy guard blocks identifier columns",
  fails(nb_export(configfile, file.path(home, "privacy_must_fail.tar.gz"), root)))
nb_write(original, report); nb_finish(s, "report-ready", receipt$key, names(receipt$outputs))
repeatfile <- file.path(cfg$regulatory$output_dir, "repeat_intervals.tsv.gz")
saved_repeat <- file.path(home, "saved_repeats.gz"); file.copy(repeatfile, saved_repeat)
altered <- reg_read(repeatfile); altered[, repeat_id := "changed_resource"]
nb_write(altered, repeatfile)
check("changed repeat resource invalidates TE results", fails(nb_te_require(cfg, s, root)))
check("changed repeat resource blocks final archive",
  fails(nb_export(configfile, file.path(home, "resource_must_fail.tar.gz"), root)))
file.copy(saved_repeat, repeatfile, overwrite = TRUE)
check("exact restored repeat resource verifies", !is.null(nb_te_require(cfg, s, root)))
old_contrasts <- getOption("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))
check("changed runtime contrasts invalidate prepared inputs", fails(nb_check_prepare(cfg, s, root)))
options(contrasts = old_contrasts)
batchfile <- reg_read(nb_path(s, "transcript_batches.tsv"))[n_features > 0]$input_file[1]
saved_batch <- file.path(home, "saved_batch.rds"); file.copy(batchfile, saved_batch)
changed_batch <- readRDS(batchfile); changed_batch$data[, maternal_count := maternal_count + 1]
saveRDS(changed_batch, batchfile)
check("mutated private batch blocks within-gene data retrieval",
  fails(nb_gene_data(s, txmap$transcript_id)))
file.copy(saved_batch, batchfile, overwrite = TRUE)
check("exact restored private batches verify", nrow(nb_gene_data(s, txmap$transcript_id)) > 0)
unlink(nb_path(s, "gene_POE.tsv"))
check("missing new model blocks archive instead of fallback",
  fails(nb_export(configfile, file.path(home, "must_fail.tar.gz"), root)))
# Second scenario: no estimable features, no fish repeats, no within-gene jobs.
cfg$profile <- "public_cross"; cfg$species <- "fish"
metadata[, cross := rep(c("focal_mother", "alternate_mother"), each = 20)]
s$output_dir <- file.path(home, "missing_repeats_run")
s$min_samples <- 100L  # Deliberate software failure-path gate, not production.
nb_prepare(cfg, s, root)
for (lv in c("gene", "transcript")) {
  manifest <- reg_read(nb_path(s, paste0(lv, "_batches.tsv")))
  for (b in manifest$batch) nb_fit_batch(cfg, s, root, lv, b)
}
nb_combine(cfg, s, root); nb_te_prepare(cfg, s, root); nb_te_combine(cfg, s, root)
nb_between_genes(cfg, s, root); nb_context(cfg, s, root); nb_report_ready(cfg, s, root)
check("missing fish repeats explicitly unavailable",
  reg_read(nb_path(s, "TE_availability.tsv"))$status == "TE_comparisons_unavailable_missing_repeats")
check("zero mixed-gene jobs combine without invented comparisons",
  nrow(reg_read(nb_path(s, "within_gene_contrasts.tsv.gz"))) == 0)
check("all filtered features remain in model tables",
  nrow(reg_read(nb_path(s, "transcript_POE.tsv"))) == 2 &&
    all(reg_read(nb_path(s, "transcript_POE.tsv"))$status == "insufficient_samples"))
nb_export(configfile, file.path(home, "unavailable_fixture_results.tar.gz"), root)
check("unavailable TE and no-eligible-feature archive is valid",
  file.info(file.path(home, "unavailable_fixture_results.tar.gz"))$size > 0)
cat("Completed ", checks, " end-to-end checks. Software fixture only.\n", sep = "")
unlink(home, recursive = TRUE)
