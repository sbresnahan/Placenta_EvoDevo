#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1]) else normalizePath(".")
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages(library(data.table))
for (m in c("common", "config", "regulatory_common", "baseline_nb", "baseline_te", "baseline_handoff",
    "window_prepare", "exposure_inputs", "te_expression"))
  source(file.path(root, "R", paste0(m, ".R")))
checks <- 0L
check <- function(label, expression) {
  stopifnot(isTRUE(expression)); checks <<- checks + 1L; cat("PASS ", label, "\n", sep = "")
}
fails <- function(expression) inherits(tryCatch({force(expression); NULL}, error = identity), "error")
tx <- data.table(transcript_id = c("t1", "t2", "t3", "t4"), gene_id = c("g", "g", "h", "h"),
  chrom = "c", tss = c(1, 5000, 9000, 15000))
ct <- data.table(chrom = "c", length = 16000)
rr <- data.table(chrom = "c", start = c(3001, 3002, 8500), end = c(3001, 3010, 8600),
  is_te = TRUE, is_unknown = FALSE, repeat_id = c("a", "b", "d"),
  repeat_family = c("F", "F", "G"), repeat_class = "LINE")
a <- nb_te_annotation(tx, rr, TRUE, 3000, ct)
check("inclusive clipped windows", a$promoter_start[1] == 1 && a$promoter_end[1] == 3001 &&
  a$promoter_end[4] == 16000)
check("exact annotation retained", all(a$tss_te_label == "nonTE"))
check("active annotation and distance", all(a$promoter_te_label == c("TE", "TE", "TE", "nonTE")) &&
  a$promoter_min_TE_distance[1] == 3000)
check("radius zero sensitivity", all(nb_te_annotation(tx, rr, TRUE, 0, ct)$promoter_te_label == "nonTE"))
check("one kb sensitivity", identical(nb_te_annotation(tx, rr, TRUE, 1000, ct)$promoter_te_label,
  c("nonTE", "nonTE", "TE", "nonTE")))
check("native lengths required", fails(reg_tss_windows(tx, 3000)))
unknown <- copy(rr); unknown[1, is_unknown := TRUE]
check("unknown repeat precedence", nb_te_annotation(tx, unknown, TRUE, 3000, ct)$promoter_te_label[1] == "unknown_repeat")
loci <- copy(rr)[, `:=`(family_id = repeat_family, locus_id = repeat_id)]
p <- ex_match_targets(a, loci)
check("family deduplication", p[level == "transcript" & feature_id == "t2", .N] == 1 &&
  p[level == "transcript" & feature_id == "t2"]$n_insertions == 2)
check("insertion distance provenance", all(p$proximal_only) && all(!p$exact_TSS))
check("all insertion rows retain coordinates and distance", nrow(attr(p, "insertion_links")) == 4L &&
  all(c("insertion_start", "insertion_end", "distance", "transcript_tss") %in%
    names(attr(p, "insertion_links"))))
check("no downstream keys in reuse contract", !any(c("regulatory", "exposure_nb") %in%
  names(nb_poe_contract(list(), list()))))
home <- tempfile("poe-pointer-contract-"); dir.create(home)
cfg <- list(cohort_id = "SOFTWARE_CONTRACT", species = "mouse", profile = "public_cross",
  paths = list(model_gene = file.path(home, "model.tsv"), covariates = file.path(home, "meta.tsv"),
    scratch = file.path(home, "scratch")),
  reference = list(build = "software"), dosage = list(maternal = 1, paternal = 1),
  regulatory = list(te_tss_window_bp = 3000), environment_allowlist = character(),
  baseline_nb = list(output_dir = file.path(home, "source")))
attr(cfg, "file") <- file.path(home, "source.yaml")
s <- nb_settings(cfg)
for (f in NB_POE_LAYER_FILES) nb_write(data.table(software_contract = TRUE), nb_path(s, f))
nb_save(data.table(software_contract = TRUE), nb_path(s, "private", "metadata.rds"))
pointer <- file.path(home, "immutable.rds")
saveRDS(list(data = data.table(software_contract = TRUE)), pointer)
for (lv in c("gene", "transcript")) {
  nb_write(data.table(batch = 1L, input_file = pointer), nb_path(s, paste0(lv, "_batches.tsv")))
  nb_write(data.table(feature_id = "fixture", batch = 1L), nb_path(s, paste0(lv, "_feature_batches.tsv")))
}
nb_write(data.table(software_contract = TRUE), nb_path(s, "within_gene_contrasts.tsv.gz"))
nb_finish(s, "prepare", list(software = TRUE), nb_path(s, NB_POE_LAYER_FILES))
nb_finish(s, "combine", list(software = TRUE), nb_path(s, c("gene_POE.tsv", "transcript_POE.tsv")))
dest <- cfg; dest$baseline_nb$output_dir <- file.path(home, "dest")
attr(dest, "file") <- file.path(home, "dest.yaml")
ds <- nb_settings(dest)
nb_prepare_window(dest, ds, cfg, root)
check("count pointers remain immutable", identical(reg_read(nb_path(ds, "transcript_batches.tsv"))$input_file, pointer))
check("downstream outputs excluded", !file.exists(nb_path(ds, "within_gene_contrasts.tsv.gz")) &&
  !dir.exists(nb_path(ds, "private", "batches")))
check("destination receipts point to destination", all(startsWith(names(nb_require(ds, "prepare")$outputs), ds$output_dir)))
check("nonempty destination refused", fails(nb_prepare_window(dest, ds, cfg, root)))
bad <- dest; bad$dosage$paternal <- 2
check("incompatible dosage refused", fails(nb_prepare_window(bad, ds, cfg, root)))
unlink(home, recursive = TRUE)
cat("Promoter contract checks: ", checks, "\n", sep = "")
