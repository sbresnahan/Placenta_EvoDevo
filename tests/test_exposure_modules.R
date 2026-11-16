#!/usr/bin/env Rscript
# Hand-specified arithmetic and real package-data fixtures only.
# Artificial identifiers/labels below are software fixtures, never POE evidence.
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
cfg <- list(profile = "human_trio", cohort_id = "GUSTO", species = "human",
  reference = list(build = "GRCh38"), dosage = list(maternal = 1, paternal = 1))
home <- tempfile("ex-tests-"); dir.create(home)
cfg$exposure_nb <- list(output_dir = file.path(home, "out"), baseline_dir = file.path(home, "base"),
  sample_crosswalk = file.path(home, "crosswalk.tsv"), te_loci = file.path(home, "loci.tsv"),
  te_counts = file.path(home, "counts.tsv"), quant_qc = file.path(home, "qc.tsv"))
s <- ex_settings(cfg)
check("new output namespace", s$version == EX_VERSION)
bad <- cfg; bad$exposure_nb$output_dir <- file.path(home, "base", "new")
check("reject writing inside baseline", fails(ex_settings(bad)))
check("PFAS rejects non-GUSTO", fails(ex_require_human(modifyList(cfg, list(cohort_id = "mouse")))))

# Explicit arithmetic inputs test the contracts and exclusion formulas.
cw <- data.table(sample_id = paste0("fixture", 1:4), unit_id = paste0("u", 1:4),
  cluster_id = paste0("c", 1:4), SubjectID = paste0("subject", 1:4),
  sex = c("F", "M", "F", "M"), ancestry = c("A", "A", "B", "B"),
  independence_verified = TRUE)
check("valid explicit crosswalk", nrow(ex_validate_crosswalk(cw, cfg)) == 4)
check("duplicate RNA IDs rejected", fails(ex_validate_crosswalk(rbind(cw, cw[1]), cfg)))
bad <- copy(cw); bad$unit_id[2] <- bad$unit_id[1]
check("technical replicates rejected", fails(ex_validate_crosswalk(bad, cfg)))
dict <- data.table(exposure = EX_PFAS, biospecimen = "maternal blood",
  collection_time = "software fixture declaration", original_units = "fixture units",
  pre_z_transform = "none", z_reference = "fixture", z_center = 0, z_scale = 1,
  detection_limit_handling = "none", prior_imputation = "none")
check("complete dictionary accepted", !fails(ex_validate_dictionary(dict)))
bad <- copy(dict); bad$pre_z_transform[1] <- "unknown"
check("undeclared preprocessing blocks PFAS", fails(ex_validate_dictionary(bad)))
bad <- copy(dict); bad$z_scale[1] <- 0
check("invalid standardization blocks PFAS", fails(ex_validate_dictionary(bad)))
pf <- data.table(SubjectID = cw$SubjectID)
for (nm in EX_PFAS) pf[, (nm) := c(-1, 0, .5, 1.5)]
joined <- ex_join_pfas(cw, pf, dict, cfg)
check("PFAS z-scores preserved", identical(joined$PFBA_mat, pf$PFBA_mat))
bad <- copy(pf); bad$PFBA_mat[1] <- Inf
check("nonfinite PFAS rejected", fails(ex_join_pfas(cw, bad, dict, cfg)))
bad <- copy(pf); bad$PFBA_mat[1] <- NA_real_
common <- ex_join_pfas(cw, bad, dict, cfg, TRUE)
check("common-panel mode excludes incomplete panel", all(is.na(as.matrix(common[1, ..EX_PFAS]))))
loci <- data.table(locus_id = c("local", "remote", "other"), species = "human", build = "GRCh38",
  chrom = "chr1", start = c(100, 10000, 20000), end = c(199, 10100, 20100), strand = "+",
  family_id = c("F1", "F1", "F2"), repeat_name = c("sub1", "sub2", "sub3"),
  repeat_class = "LTR", source = "arithmetic fixture")
check("locus schema accepted", !fails(ex_validate_loci(loci, cfg)))
bad <- copy(loci); bad$build[1] <- "hg19"
check("build mismatch rejected", fails(ex_validate_loci(bad, cfg)))
annotation <- data.table(transcript_id = c("txTE", "txNon"), gene_id = "g1", chrom = "chr1",
  start = 100, end = 500, tss = c(100, 400), tss_te_label = c("TE", "nonTE"))
targets <- ex_match_targets(annotation, loci); spans <- ex_gene_spans(annotation)
check("only exact TE TSS creates matching pair", nrow(targets) == 2 && all(targets$family_id == "F1"))
qc <- data.table(sample_id = cw$sample_id, mapped_fragments = 1e6, status = "complete",
  species = "human", build = "GRCh38", layout = "PE", strandedness = "reverse",
  quantifier = "arithmetic fixture", quantifier_version = "fixture", multimap_limit = 100,
  annotation_id = "fixture_catalogue_v1", count_scope = "all_catalogue_loci_sparse_zero",
  unique_count_definition = "genomic_NH1")
counts <- CJ(sample_id = cw$sample_id, locus_id = loci$locus_id)
counts[, count := c(10, 5, 2, 20, 6, 4, 30, 7, 6, 40, 8, 8)]
counts[, unique_count := count / 2]
check("valid sparse TE inputs", !fails(ex_validate_quant(counts, qc, loci, cw, cfg, s)))
bad_qc <- copy(qc); bad_qc[, unique_count_definition := NULL]
check("annotation uniqueness cannot stand in for NH1", fails(ex_validate_quant(counts, bad_qc, loci, cw, cfg, s)))
bad <- copy(counts); bad$count[1] <- -1
check("negative TE counts rejected", fails(ex_validate_quant(bad, qc, loci, cw, cfg, s)))
pred <- ex_predictors(counts, qc, loci, targets, spans, s)
remote <- counts[locus_id == "remote"][order(sample_id)]$count
check("target gene locus excluded", identical(pred[predictor == "target_excluded"]$family_count, remote))
check("z-score centered across matched samples", abs(mean(pred[predictor == "target_excluded"]$TE_abundance)) < 1e-12)
check("z-score unit scale", abs(sd(pred[predictor == "target_excluded"]$TE_abundance) - 1) < 1e-12)
check("inclusive sensitivity retained", all(c("inclusive", "intergenic", "unique_target_excluded") %in% pred$predictor))
gtf <- file.path(home, "te.gtf"); aqc <- ex_write_te_gtf(loci, gtf)
check("GTF segments have unique ownership", aqc$n_segments == 3 && length(readLines(gtf)) == 3)
overlap <- rbind(loci, transform(loci[1], locus_id = "nested", start = 120, end = 140, family_id = "conflict"))
aqc <- ex_write_te_gtf(overlap, gtf)
check("ambiguous nested family segment excluded", aqc$ambiguous_segments_excluded == 1)
ci <- ex_product_interval(-2, 1, -3, 4)
check("product interval uses all corners", identical(unname(ci), c(-8, 6)))
check("reversed intervals rejected", fails(ex_product_interval(2, 1, 0, 1)))
check("absolute interval includes zero across sign", identical(unname(ex_abs_range(-2, 3)), c(0, 3)))

# Observed Salamanders counts/environmental values are used only to verify
# numerical fitting. Assigned allele/sex/ancestry labels are not biology.
stopifnot(requireNamespace("glmmTMB", quietly = TRUE))
data("Salamanders", package = "glmmTMB")
pairs <- CJ(sample_id = sprintf("real_fixture_%02d", 1:40), SNPpos = paste0("chr1:", 1:4, ":A:G"))
pairs[, `:=`(feature_id = "numerical_fixture", maternal_count = Salamanders$count[seq_len(.N)],
  paternal_count = Salamanders$count[seq_len(.N) + 160], offset = 0)]
meta <- data.table(sample_id = sprintf("real_fixture_%02d", 1:40), cluster_id = sprintf("f%02d", 1:40),
  independence_verified = TRUE, sex = rep(c("F", "M"), 20), ancestry = rep(c("A", "A", "B", "B"), 10),
  E = Salamanders$Wtemp[1:40], M = Salamanders$cover[1:40])
check("NB long table retains exposures", all(c("E", "M") %in% names(ex_long(pairs, meta))))
test <- ex_nb(pairs, meta, cfg, s, "h6")
check("real-data numerical H6 fixture converges", test$row$status == "eligible")
check("H6 exposure interaction coefficient extracted", is.finite(test$row$slope))
check("magnitude follows absolute endpoint definition", abs(test$row$delta_magnitude -
  (abs(test$row$beta_p90) - abs(test$row$beta_p10))) < 1e-10)
check("dosage adjustment changes prediction only", {
  cfg_maize_dosage <- cfg; cfg_maize_dosage$dosage <- list(maternal = 2, paternal = 1)
  alt <- ex_nb(pairs, meta, cfg_maize_dosage, s, "h6")
  abs(alt$row$slope - test$row$slope) < 1e-10 &&
    abs(alt$row$beta_p50 - test$row$beta_p50 - log(2)) < 1e-8
})
bad <- copy(meta); bad$E <- 1
check("constant exposure unavailable", ex_nb(pairs, bad, cfg, s, "h6")$row$status == "constant_or_insufficient_exposure_values")
bad <- copy(meta); bad$M <- bad$E
check("collinear mediation unavailable", ex_nb(pairs, bad, cfg, s, "mediation")$row$status == "rank_deficient_sample_design")
med <- ex_mediation(pairs, meta, cfg, s)
check("real-data numerical mediation fixture converges", med$row$status == "eligible")
check("indirect product arithmetic", abs(med$row$indirect_beta - med$row$a * med$row$b) < 1e-12)
check("joint-significance p-value", med$row$indirect_path_p == max(med$row$a_p, med$row$b_p))
check("direct and indirect signed scenarios add", abs(med$row$natural_direct_beta +
  med$row$natural_indirect_beta - med$row$natural_total_beta) < 1e-10)
check("direct and indirect magnitude scenarios add", abs(med$row$natural_direct_magnitude +
  med$row$natural_indirect_magnitude - med$row$natural_total_magnitude) < 1e-10)
check("absolute indirect is not absolute signed indirect", abs(med$row$natural_indirect_magnitude -
  abs(med$row$natural_indirect_beta)) > 1e-10)
check("nonhuman PFAS model blocked", fails(ex_mediation(pairs, meta, modifyList(cfg, list(cohort_id = "mouse")), s)))
mod <- ex_mediation(pairs, meta, cfg, s, modification = TRUE)
check("optional modification numerical fixture converges", mod$row$status == "eligible")
check("modification retains all lower-order terms", all(c("E:M", "allelepaternal:E",
  "allelepaternal:M", "allelepaternal:E:M") %in% names(mod$coefficients)))
check("simple product not misreported under modification", is.na(mod$row$indirect_beta) &&
  is.na(mod$row$indirect_path_p) && is.finite(mod$row$natural_indirect_beta))
check("modification standardization decomposes signed total", abs(mod$row$natural_direct_beta +
  mod$row$natural_indirect_beta - mod$row$natural_total_beta) < 1e-10)
bad <- copy(meta); bad[, independence_verified := FALSE]
check("unverified independence blocks fitting", ex_nb(pairs, bad, cfg, s)$row$status == "unverified_biological_independence")
bad <- copy(meta); bad[, sex := NA_character_]
check("missing human covariate not silently dropped", ex_nb(pairs, bad, cfg, s)$row$status == "missing_sex")
public_cfg <- modifyList(cfg, list(profile = "public_cross", cohort_id = "fixture_cross", species = "mouse"))
public_meta <- copy(meta); public_meta[, cross := rep(c("focal_mother", "alternate_mother"), each = 20)]
bad <- copy(public_meta); bad[, E := as.numeric(factor(cross))]
check("between-cross variation cannot identify H6", ex_nb(pairs, bad, public_cfg, s, "h6")$row$status ==
  "no_within_cross_exposure_variation")
check("small-cohort unit gate respected", ex_nb(pairs[sample_id %in% meta$sample_id[1:3]],
  meta[1:3], cfg, s, "h6")$row$status == "insufficient_biological_units")

d <- data.table(module = rep("h7_mediation", 4), cohort_id = "GUSTO", level = "transcript",
  predictor = "target_excluded", exposure = rep(EX_PFAS[1:2], each = 2),
  status = "eligible", indirect_path_p = c(.01, .04, .02, .9), slope_p = c(.5, .5, .5, .5))
finished <- ex_finish(d)
check("H7 FDR pools compounds", identical(finished$indirect_q, p.adjust(d$indirect_path_p, "BH")))
check("within-compound FDR separate", all(finished$indirect_q_within_exposure[1:2] == p.adjust(d$indirect_path_p[1:2], "BH")))

# Minimal paired-refit check only; production default remains 2,000.
s_fixture <- s; s_fixture$bootstrap_min_valid <- .01
boot <- ex_bootstrap_pair(pairs, meta, cfg, s_fixture, B = 2L)
check("paired bootstrap returns every draw", nrow(boot$draws) == 2L)
check("paired bootstrap reports refit counts", boot$summary$B == 2L && boot$summary$n_valid >= 1L)

if (length(args) >= 2L) {
  dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(check = checks, status = "pass"), file.path(args[2], "exposure_unit_tests.tsv"), sep = "\t")
  saveRDS(list(pairs = pairs, metadata = meta, cfg = cfg, s = s, mediation = med$row),
    file.path(args[2], "numerical_fixture_private.rds"))
}
cat("Completed ", length(checks), " exposure software checks. No biological results were generated.\n", sep = "")
