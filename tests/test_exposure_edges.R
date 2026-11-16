#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
suppressPackageStartupMessages(library(data.table))
for (m in c("exposure_inputs", "te_expression", "exposure_nb", "exposure_mediation",
            "exposure_followup", "exposure_workflow"))
  source(file.path(root, "R", paste0(m, ".R")))
checks <- character()
check <- function(label, expr) {
  if (!isTRUE(expr)) stop("FAIL ", label)
  checks <<- c(checks, label); cat("PASS ", label, "\n", sep = "")
}
fails <- function(expr) inherits(tryCatch({force(expr); NULL}, error = function(e) e), "error")
home <- tempfile("ex-edges-"); dir.create(home)
cfg <- list(profile = "public_cross", cohort_id = "fixture", species = "mouse",
  reference = list(build = "GRCm38"), dosage = list(maternal = 1, paternal = 1),
  exposure_nb = list(output_dir = file.path(home, "out"), baseline_dir = file.path(home, "base"),
    sample_crosswalk = file.path(home, "cw.tsv"), te_loci = file.path(home, "loci.tsv"),
    te_counts = file.path(home, "counts.tsv"), quant_qc = file.path(home, "qc.tsv"),
    libraries = file.path(home, "libraries.tsv")))
s <- ex_settings(cfg)
loci <- data.table(locus_id = "fixture-locus", species = "mouse", build = "GRCm38", chrom = "chr1",
  start = 100, end = 199, strand = "+", family_id = "Contains_n_and_t", repeat_name = "fixture",
  repeat_class = "fixture", source = "fixture")
check("native family names may contain n and t", !fails(ex_validate_loci(loci, cfg)))
ex_write(loci, s$te_loci)
lib <- data.table(sample_id = c("PE_reverse", "SE_forward", "unstranded"),
  bam = "/not/executed/fixture.bam", fastq1 = NA_character_, fastq2 = NA_character_,
  layout = c("PE", "SE", "PE"), strandedness = c("reverse", "forward", "unstranded"),
  star_index = NA_character_, telescope_report = "/not/executed/report.tsv",
  unique_telescope_report = "/not/executed/NH1_report.tsv", mapped_fragments = 1e6)
ex_write(lib, s$libraries)
ex_te_plan(cfg, s)
scripts <- file.path(s$output_dir, "quantification_plan", paste0(lib$sample_id, ".sh"))
check("generated quant scripts parse", all(vapply(scripts, function(p)
  system2("bash", c("-n", shQuote(p))) == 0L, logical(1))))
txt <- lapply(scripts, readLines)
check("PE reverse strand maps RF", any(grepl("--stranded_mode RF", txt[[1]], fixed = TRUE)))
check("SE forward strand maps F", any(grepl("--stranded_mode F ", txt[[2]], fixed = TRUE)))
check("unstranded omits strand option", !any(grepl("--stranded_mode", txt[[3]], fixed = TRUE)))
check("exact quantifier commit guarded", any(grepl(s$quantifier_commit, txt[[1]], fixed = TRUE)))
check("BAM tags probed before assignment", any(grepl("check_exposure_bam.R", txt[[1]], fixed = TRUE)))
check("strict NH1 filtering generated", any(grepl("NH:i:1", txt[[1]], fixed = TRUE)))
check("quantification not marked executed", all(!ex_read(ex_path(s, "quantification_plan", "jobs.tsv"))$executed))
bam <- system.file("extdata", "ex1.bam", package = "Rsamtools")
check("BAM missing NH HI AS nM rejected", nzchar(bam) && fails(ex_bam_probe(bam, 100)))

# Explicit empty-job object: test state handling without inventing cohort data.
jobs <- data.table(job_id = character(), batch = integer(), module = character(), level = character(),
  feature_id = character(), gene_id = character(), family_id = character(), exposure = character(),
  predictor = character())
ex_run_lock(s, cfg, s$te_loci)
ex_save(list(jobs = jobs), ex_path(s, "private", "prepared.rds"))
ex_write(jobs, ex_path(s, "private", "jobs.tsv"))
empty <- ex_combine(cfg, s)
check("empty model inventory combines safely", nrow(empty) == 0L)
check("empty status explicit", ex_status(s)$planned_batches == 0L)

# Multiplicity checks are arithmetic fixtures, not biological estimates.
d <- data.table(job_id = c("a", "b", "a", "b"), module = "h6",
  cohort_id = rep(c("cohort1", "cohort2"), each = 2), level = "transcript",
  predictor = "target_excluded", exposure = "matching_TE", status = "eligible",
  slope_p = c(.01, .03, .1, NA_real_))
d1 <- ex_finish(d[cohort_id == "cohort1"])
check("single cohort cannot claim pooled q", all(is.na(d1$slope_q_pooled_cohorts)))
paths <- file.path(home, c("cohort1.tsv", "cohort2.tsv"))
ex_write(d[cohort_id == "cohort1"], paths[1]); ex_write(d[cohort_id == "cohort2"], paths[2])
pool <- ex_pool_h6(paths, file.path(home, "pooled.tsv"))
check("cross-cohort BH uses all planned rows", isTRUE(all.equal(pool$slope_q_pooled_cohorts,
  p.adjust(d$slope_p, "BH", n = nrow(d)))))
check("pool rejects duplicate cohort jobs", fails(ex_pool_h6(c(paths, paths[1]), file.path(home, "bad.tsv"))))
g <- copy(d); g[, `:=`(gene_id = "fixture_gene", primary_p = slope_p,
  total_counts = 50, n_biological_units = 4, n_snps = 3)]
spans <- data.table(gene_id = "fixture_gene", start = 100, end = 200)
gs <- ex_gene_scores(g, spans)
check("gene omnibus retains planned family multiplicity",
  gs[cohort_id == "cohort2"]$gene_p == .2)
if (length(args) >= 2L) {
  dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(check = checks, status = "pass"), file.path(args[2], "exposure_edge_tests.tsv"), sep = "\t")
}
cat("Completed ", length(checks), " edge checks.\n", sep = "")
