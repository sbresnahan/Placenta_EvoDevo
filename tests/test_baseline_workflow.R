#!/usr/bin/env Rscript
# End-to-end SOFTWARE harness, not biological POE data.
# Observed Salamanders example counts are reused with artificial labels.
# Raw acquisition, parent assignment and site inputs are mocked explicitly.
a <- commandArgs(trailingOnly = TRUE)
root <- if (length(a)) normalizePath(a[1]) else normalizePath(".")
source(file.path(root, "R/library_paths.R"))
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
paths$model_gene <- file.path(home, "promoter_model", "model_gene.tsv")
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
  regulatory = list(output_dir = file.path(home, "resources"), bootstrap = 200, te_tss_window_bp = 0L,
    te_genes = list(placebo_shifts = 1, permutations = 1)),
  environment_allowlist = character())
configfile <- file.path(home, "fixture.yaml")
yaml::write_yaml(list(software_fixture = TRUE), configfile); attr(cfg, "file") <- configfile
# Hypothesis #5 fixture: one lead cis-QTL variant at chr1:12, inside the fixture TE
# (chr1:8-15) that overlaps fixture_txA's TSS (10); fixture_txB's TSS (100) is non-TE.
qtl_dir <- file.path(home, "qtl"); dir.create(qtl_dir)
nb_write(data.table(modality = "isoforms", phenotype_id = "fixture_gene__fixture_txA",
  variant_id = "chr1:12:A:G", chrom = "chr1", pos = 12L, slope = 0.5, slope_se = 0.1),
  file.path(qtl_dir, "fixture_lead.tsv"))
nb_write(data.table(input_id = "fixture_lead", path = file.path(qtl_dir, "fixture_lead.tsv"),
  build = "hg38", kind = "lead", ancestry = "fixture", coordinate_system = "one_based"),
  file.path(qtl_dir, "qtl_inputs.tsv"))
nb_write(data.table(modality = "isoforms", phenotype_id = "fixture_gene__fixture_txA",
  gene_id = "fixture_gene", event_id = "fixture_event1", mapping_level = "event",
  mapping_status = "verified"), file.path(qtl_dir, "phenotype_map.tsv"))
cfg$regulatory$human_context <- list(require_qtl = TRUE,
  qtl_tables = file.path(qtl_dir, "qtl_inputs.tsv"),
  phenotype_map = file.path(qtl_dir, "phenotype_map.tsv"))
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
check("combined model schema", all(reg_read(nb_path(s, "transcript_POE.tsv"))$version == NB_VERSION))
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
noqtl <- cfg; noqtl$regulatory$human_context$qtl_tables <- NULL
gate <- tryCatch({nb_context(noqtl, s, root); "ran"}, error = function(e) conditionMessage(e))
check("require_qtl blocks human context without QTL tables", grepl("require_qtl", gate))
nb_context(cfg, s, root)
ctx <- reg_read(nb_path(s, "context", "human_transcript_context.tsv"))
check("QTL variant inside the TSS TE insertion is counted",
  ctx[transcript_id == "fixture_txA"]$n_qtl_variants_in_TSS_TE == 1)
check("untested transcript retains missing context",
  is.na(ctx[transcript_id == "fixture_txB"]$n_qtl_variants_in_TSS_TE))
check("verified event link counted", ctx[transcript_id == "fixture_txA"]$n_event_qtl_links == 1)
check("lead-only QTL leaves max PIP missing", is.na(ctx[transcript_id == "fixture_txA"]$max_event_link_pip))
qlinks <- reg_read(nb_path(s, "context", "human_qtl_links.tsv.gz"))
check("link flags TE-insertion overlap but not exact TSS",
  isTRUE(qlinks$variant_overlaps_TSS_TE[1]) && identical(qlinks$variant_at_transcript_TSS[1], FALSE))
check("own-isoform PANTRY phenotype flagged as exact transcript match",
  isTRUE(qlinks$exact_transcript_phenotype[1]) && identical(qlinks$pantry_transcript_id[1], "fixture_txA"))
check("exact isoform QTL counted per transcript", ctx[transcript_id == "fixture_txA"]$n_exact_isoform_qtl_links == 1)
check("PANTRY transcript parsing per modality", identical(nb_pantry_transcript(
  c("isoforms", "isoform_expression", "alt_TSS", "expression", "alt_TSS"),
  c("G__TX1.2", "G__TX2", "G__chr1_100_upstream_TX3.1", "G", "G__no_marker")),
  c("TX1", "TX2", "TX3", NA, NA)))
check("context provenance written", file.exists(nb_path(s, "context", "human_context_provenance.tsv")))
nb_report_ready(cfg, s, root)
archive <- file.path(home, "software_fixture_results.tar.gz")
nb_export(configfile, archive, root)
listing <- utils::untar(archive, list = TRUE)
check("archive contains gene POE model", any(endsWith(listing, "/gene_POE.tsv")))
check("private products excluded", !any(grepl("/private/|metadata.rds|promoter_model", listing)))
extract <- file.path(home, "unpacked"); dir.create(extract); utils::untar(archive, exdir = extract)
mf <- reg_read(list.files(extract, pattern = "^MANIFEST.tsv$", recursive = TRUE, full.names = TRUE)[1])
top <- file.path(extract, "software_fixture_results")
check("exported file sizes match", all(file.info(file.path(top, mf$archive_path))$size == mf$bytes))
check("all stage statuses verify", all(nb_status(cfg, s, root)$status == "complete_verified"))
# A valid receipt cannot authorize a private identifier column for export.
report <- nb_path(s, "report_summary.tsv")
original <- reg_read(report); receipt <- nb_require(s, "report-ready")
bad <- copy(original); bad[, sample_id := "SOFTWARE_PRIVACY_FIXTURE"]
nb_write(bad, report); nb_finish(s, "report-ready", receipt$key, names(receipt$outputs))
check("privacy guard blocks identifier columns",
  fails(nb_export(configfile, file.path(home, "privacy_must_fail.tar.gz"), root)))
nb_write(original, report); nb_finish(s, "report-ready", receipt$key, names(receipt$outputs))
# Completed receipts use required-output presence by design.
unlink(nb_path(s, "gene_POE.tsv"))
check("missing required model blocks archive",
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
nb_combine(cfg, s, root)
# Partial export: POE layer of a cohort whose TE stages have not run yet.
check("full export refuses a cohort with pending stages",
  fails(nb_export(configfile, file.path(home, "pending_must_fail.tar.gz"), root)))
check("partial export rejects unknown cohort IDs",
  fails(nb_export(configfile, file.path(home, "bad_partial.tar.gz"), root, partial = "NOT_A_COHORT")))
partial_archive <- file.path(home, "partial_fixture_results.tar.gz")
nb_export(configfile, partial_archive, root, partial = cfg$cohort_id)
listing <- utils::untar(partial_archive, list = TRUE)
check("partial export carries the POE layer", any(endsWith(listing, "/gene_POE.tsv")) &&
  any(endsWith(listing, "/transcript_POE.tsv")))
check("partial export omits pending TE outputs", !any(endsWith(listing, "/TE_availability.tsv")))
check("partial export records pending stages", any(endsWith(listing, "/STAGE_STATUS.tsv")))
nb_te_prepare(cfg, s, root); nb_te_combine(cfg, s, root)
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
