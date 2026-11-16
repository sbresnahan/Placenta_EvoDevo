#!/usr/bin/env Rscript
# Observed glmmTMB data, relabelled exclusively as software fixtures.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages(library(data.table))
for (m in c("exposure_inputs", "te_expression", "exposure_nb", "exposure_mediation",
            "exposure_followup", "exposure_workflow", "exposure_plots"))
  source(file.path(root, "R", paste0(m, ".R")))
checks <- character()
check <- function(label, expr) {
  if (!isTRUE(expr)) stop("FAIL ", label)
  checks <<- c(checks, label); cat("PASS ", label, "\n", sep = "")
}
fails <- function(expr) inherits(tryCatch({force(expr); NULL}, error = function(e) e), "error")
data("Salamanders", package = "glmmTMB")
home <- tempfile("ex-workflow-"); dir.create(home)
cfg <- list(profile = "human_trio", cohort_id = "GUSTO", species = "human",
  reference = list(build = "GRCh38"), dosage = list(maternal = 1, paternal = 1),
  regulatory = list(te_tss_window_bp = 0L))
cfg$exposure_nb <- list(output_dir = file.path(home, "out"), baseline_dir = file.path(home, "base"),
  sample_crosswalk = file.path(home, "cw.tsv"), te_loci = file.path(home, "loci.tsv"),
  te_counts = file.path(home, "counts.tsv"), quant_qc = file.path(home, "qc.tsv"),
  pfas = file.path(home, "pf.tsv"), exposure_dictionary = file.path(home, "dict.tsv"),
  batch_jobs = 7L)
s <- ex_settings(cfg)
meta <- data.table(sample_id = sprintf("fixture%02d", 1:40), unit_id = sprintf("u%02d", 1:40),
  cluster_id = sprintf("family%02d", 1:40), SubjectID = sprintf("subject%02d", 1:40),
  independence_verified = TRUE, sex = rep(c("F", "M"), 20),
  ancestry = rep(c("A", "A", "B", "B"), 10))
ex_write(meta, s$sample_crosswalk)
pairs <- CJ(sample_id = meta$sample_id, SNPpos = paste0("chr1:", 1:4, ":A:G"))
pairs[, `:=`(feature_id = "txFixture", maternal_count = Salamanders$count[1:160],
  paternal_count = Salamanders$count[161:320], offset = 0)]
dir.create(s$baseline_dir)
for (lv in c("gene", "transcript")) {
  z <- copy(pairs); id <- if (lv == "gene") "geneFixture" else "txFixture"
  z[, feature_id := id]
  path <- file.path(s$baseline_dir, paste0(lv, ".rds"))
  saveRDS(list(data = z), path)
  ex_write(data.table(batch = 1L, level = lv, input_file = path,
    object_bytes = as.numeric(object.size(z))), file.path(s$baseline_dir, paste0(lv, "_batches.tsv")))
  ex_write(data.table(feature_id = id, batch = 1L),
    file.path(s$baseline_dir, paste0(lv, "_feature_batches.tsv")))
}
ann <- data.table(transcript_id = "txFixture", gene_id = "geneFixture", chrom = "chr1",
  start = 100, end = 500, tss = 100, tss_te_label = "TE",
  promoter_te_label = "TE", promoter_start = 100, promoter_end = 100, te_tss_window_bp = 0L)
ex_write(ann, file.path(s$baseline_dir, "transcript_TE_annotation.tsv"))
loci <- data.table(locus_id = c("local", "remote"), species = "human", build = "GRCh38",
  chrom = "chr1", start = c(100, 10000), end = c(199, 10100), strand = "+",
  family_id = "fixture_family", repeat_name = "fixture_subfamily", repeat_class = "fixture",
  source = "software fixture")
ex_write(loci, s$te_loci)
counts <- data.table(sample_id = meta$sample_id, locus_id = "remote",
  count = Salamanders$count[361:400], unique_count = Salamanders$count[361:400])
ex_write(counts, s$te_counts)
qc <- data.table(sample_id = meta$sample_id, mapped_fragments = 1e6, status = "complete",
  species = "human", build = "GRCh38", layout = "PE", strandedness = "reverse",
  quantifier = "observed package count fixture", quantifier_version = "fixture", multimap_limit = 100,
  annotation_id = "fixture_catalogue_v1", count_scope = "all_catalogue_loci_sparse_zero",
  unique_count_definition = "genomic_NH1")
ex_write(qc, s$quant_qc)
blocked <- ex_preflight(cfg, s)
check("absent dictionary blocks H7 not H6", !blocked$h7_ready && nrow(blocked$targets) == 2)
pf <- data.table(SubjectID = meta$SubjectID)
for (nm in EX_PFAS) pf[, (nm) := Salamanders$Wtemp[1:40]]
ex_write(pf, s$pfas)
dict <- data.table(exposure = EX_PFAS, biospecimen = "maternal fixture declaration",
  collection_time = "fixture only", original_units = "observed package variable",
  pre_z_transform = "none", z_reference = "software fixture only",
  z_center = 0, z_scale = 1, detection_limit_handling = "none", prior_imputation = "none")
ex_write(dict, s$exposure_dictionary)
pre <- ex_preflight(cfg, s)
check("dictionary enables H7", pre$h7_ready)
ex_prepare(cfg, s)
obj <- ex_load_prepared(s)
check("job IDs readable and unique", all(grepl("^job_[0-9]{8}$", obj$jobs$job_id)) &&
  !anyDuplicated(obj$jobs$job_id))
ex_prepare(cfg, s)
check("repeated preparation retains deterministic job IDs",
  identical(obj$jobs$job_id, ex_load_prepared(s)$jobs$job_id))
check("all feature PFAS and mediation jobs scheduled", nrow(obj$jobs) == 34L)
check("no significance gate on mediation", sum(obj$jobs$module == "h7_mediation") == 16L)
check("jobs bounded by requested batch size", max(obj$jobs[, .N, by = batch]$N) <= 7L)
check("H6 job exposure join", all(is.finite(ex_job_metadata(obj, obj$jobs[module == "h6"][1])$E)))
check("H7 paired exposure join", all(c("E", "M") %in%
  names(ex_job_metadata(obj, obj$jobs[module == "h7_mediation"][1]))))
check("combine blocks incomplete arrays", fails(ex_combine(cfg, s)))
for (batch_id in unique(obj$jobs$batch)) {
  mods <- obj$jobs[batch == batch_id]$module
  if (any(mods %in% c("h6", "h7_total"))) ex_fit_batch(cfg, s, batch_id)
  if (any(mods %in% c("h7_mediation", "h7_modification"))) ex_fit_batch(cfg, s, batch_id, TRUE)
}
results <- ex_combine(cfg, s)
check("every planned row survives combine", nrow(results) == 34L && setequal(results$job_id, obj$jobs$job_id))
check("workflow avoids hidden input model exceptions", !any(results$status == "input_or_model_error"))
check("status confirms all expected stage batches", all(ex_status(s)$planned_batches == ex_status(s)$completed_batches))
selection <- ex_read(ex_path(s, "bootstrap_selection.tsv"))
check("bootstrap selection bounded and deterministic", nrow(selection) <= 20L &&
  identical(selection$selection_index, seq_len(nrow(selection))))
inv <- ex_export(s)
exported <- ex_read(ex_path(s, "aggregate_export", "model_results.tsv"))
check("exports remove paths identities and error text",
  !any(grepl("sample_id|SubjectID|cluster_id|input_file|warning", names(exported))))
check("export inventory records nonempty file sizes", all(inv$bytes > 0) && identical(names(inv), c("file", "bytes")))
changed <- s; changed$seed <- s$seed + 1L
check("changed settings blocked", fails(ex_load_prepared(changed)))
ex_write(rbind(counts, counts[1]), s$te_counts)
check("changed prepared inputs blocked", fails(ex_load_prepared(s)))

# Source report shape from pinned Telescope, hand arithmetic counts.
report <- file.path(home, "telescope-report.tsv")
writeLines(c("## RunInfo\tversion:fixture", "transcript\tfinal_count\tunique_count",
  "__no_feature\t0\t0", "remote\t2.5\t2"), report)
imported <- ex_import_telescope(data.table(sample_id = "fixture01", telescope_report = report,
  mapped_fragments = 1e6), loci)
check("Telescope commented header imported by name", nrow(imported) == 1L && imported$count == 2.5)
check("annotation unique count not mislabelled unique mapper", !"unique_count" %in% names(imported) &&
  imported$annotation_unique_count == 2)
u_report <- file.path(home, "NH1-report.tsv")
writeLines(c("## RunInfo\tversion:fixture", "transcript\tfinal_count\tunique_count",
  "remote\t1\t1"), u_report)
nh1 <- ex_import_telescope(data.table(sample_id = "fixture01", telescope_report = report,
  unique_telescope_report = u_report, mapped_fragments = 1e6), loci)
check("strict NH1 report supplies unique sensitivity", nh1$unique_count == 1 && nh1$count == 2.5)

# Direct fit parity uses observed counts/environmental values, not new biology.
m <- copy(meta); m[, `:=`(E = Salamanders$Wtemp[1:40], M = Salamanders$cover[1:40])]
fit <- ex_nb(pairs, m, cfg, s)
direct <- glmmTMB::glmmTMB(count ~ allele * E + sex + ancestry + (1|SNPpos) + (1|sample_id),
  data = fit$data, family = glmmTMB::nbinom2(), offset = offset,
  control = glmmTMB::glmmTMBControl(optimizer = nlminb, parallel = 1L,
    optCtrl = list(rel.tol = 1e-5, eval.max = 300L, iter.max = 300L)))
check("H6 coefficients match direct glmmTMB", max(abs(fit$coefficients - glmmTMB::fixef(direct)$cond)) < 1e-10)
check("residual concentration diagnostics finite", all(is.finite(unlist(fit$row[,
  .(pearson_RMS, pearson_abs_max, largest_family_count_fraction, largest_family_residual_fraction)]))))
pub <- modifyList(cfg, list(profile = "public_cross", cohort_id = "fixture", species = "mouse"))
top_units <- fit$data[, .(depth = sum(count)), by = sample_id][order(-depth)]$sample_id[1:4]
pm <- m[sample_id %in% as.character(top_units), !"ancestry"]; pm[, sex := NA_character_]
pm[, `:=`(cross = rep(c("focal_mother", "alternate_mother"), each = 2), E = c(-1, 1, -2, 2))]
small <- ex_nb(pairs[sample_id %in% pm$sample_id], pm, pub, s)
check("identifiable four unit cohort attempted", small$row$status %in%
  c("eligible", "fit_error", "failed_convergence_or_Hessian", "nonfinite_uncertainty"))
saturated <- copy(pm); saturated[, sex := c("F", "F", "F", "M")]
check("saturated cohort unavailable", ex_nb(pairs[sample_id %in% pm$sample_id],
  saturated, pub, s)$row$status == "no_residual_biological_unit_support")

# Deterministic arithmetic tables verify enrichment contracts and sparse H5.
scores <- data.table(gene_id = paste0("fixture_gene_", 1:20), gene_p = (1:20) / 21,
  total_counts = 100, n_biological_units = 20, n_snps = 4, gene_length = 1000, n_tested = 1)
man <- data.table(set_id = "fixture_set", category = "fixture", name = "arithmetic",
  species = "human", source = "fixture", source_version = "1", retrieved_date = "2026-11-15")
members <- data.table(set_id = "fixture_set", gene_id = scores$gene_id[1:5], species = "human")
es <- s; es$permutations <- 19L
en <- ex_enrich_group(scores, members, man, es)
check("matched gene permutations retain planned count", en$status == "available" && en$permutations == 19L)
context <- data.table(transcript_id = "txFixture", gene_id = "geneFixture",
  n_exact_isoform_qtl_links = 1, n_event_qtl_links = 0, n_qtl_variants_in_TSS_TE = 0,
  n_peak_overlaps = 1, qtl_assessed = FALSE, chromatin_assessed = TRUE, TE_insertion_assessed = FALSE)
h5 <- ex_h5_adapter(context)
check("unassessed H5 support stays missing", is.na(h5$exact_QTL) && isTRUE(h5$cCRE))
h <- ex_h5_followup(results, h5, cfg)
check("sparse H5 enrichment flagged", all(h$status == "sparse_or_constant_annotation"))

# Render to scratch only: labelled numerical fixtures, not delivered biology.
fd <- copy(results)
ex_plot_forest(fd[module == "h6"], file.path(home, "forest.png"), mediation = FALSE)
ex_plot_beta_curve(fit, file.path(home, "curve.png"))
ex_plot_signed_magnitude(fd, file.path(home, "magnitude.png"))
ex_plot_joint_support(m, file.path(home, "support.png"))
ex_plot_pfas_correlations(pf, file.path(home, "correlation.png"))
check("five PNG helpers render", all(file.info(file.path(home,
  paste0(c("forest", "curve", "magnitude", "support", "correlation"), ".png")))$size > 0))
check("empty forest returns without fabricating plot", is.null(ex_plot_forest(fd[0], file.path(home, "empty.png"))))

# Full existing YAML remains authoritative; only new settings are merged.
full <- yaml::read_yaml(file.path(root, "config/GUSTO.example.yaml"))
full$exposure_nb <- cfg$exposure_nb
yml <- file.path(home, "fixture.yaml"); yaml::write_yaml(full, yml)
status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(file.path(root, "poe_exposure.R")),
  shQuote(yml), "status"), stdout = TRUE, stderr = TRUE)
check("standalone CLI status runs with full schema2", is.null(attr(status, "status")))
if (length(args) >= 2L) {
  dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(check = checks, status = "pass"), file.path(args[2], "exposure_workflow_tests.tsv"), sep = "\t")
  fwrite(results[, .N, by = .(module, status)],
    file.path(args[2], "workflow_fixture_coverage.tsv"), sep = "\t")
}
cat("Completed ", length(checks), " workflow checks. All identifiers are fixtures.\n", sep = "")
