#!/usr/bin/env Rscript
args <- commandArgs(TRUE); root <- normalizePath(if (length(args)) args[1] else ".")
source(file.path(root, "R/library_paths.R"))
suppressPackageStartupMessages(library(data.table))
for (m in c("common", "config", "regulatory_common", "regulatory_context", "baseline_nb",
    "baseline_handoff", "gxe_ppbmi", "exposure_inputs", "exposure_followup"))
  source(file.path(root, "R", paste0(m, ".R")))
checks <- 0L
check <- function(label, expr) {stopifnot(isTRUE(expr)); checks <<- checks + 1L; cat("PASS ", label, "\n", sep = "")}
fails <- function(expr) inherits(tryCatch({force(expr); NULL}, error = identity), "error")
home <- tempfile("gxe-format-"); dir.create(home)
cfg <- list(environment_allowlist = character()); attr(cfg, "file") <- file.path(home, "cfg.yaml")
d <- data.table(modality = "isoforms", phenotype_id = c("G__T1", "G__T2"),
  exposure = "ppBMI", ancestry = "EAS", variant_id = c("chr1:30:A:G", "chr1:500:A:G"),
  b_gi = c(.2, -.1), b_gi_se = .1, b_g = 9, b_e = 8,
  pval_nominal = c(.001, .5), pval_beta = c(.01, .8), qval = c(.04, .8))
path <- file.path(home, "gxe_tier1_results.tsv.gz"); nb_write(d, path)
spec <- data.table(input_id = "one", path = path, build = "GRCh38", tier = "tier1",
  source_format = "CHROM:POS:REF:ALT", coordinate_system = "one_based")
q <- gxe_read(spec, NULL, "chr1")
check("interaction effects use b_gi", identical(q$interaction_beta, d$b_gi) && all(q$source_b_g == 9))
check("all tested features retained", nrow(q) == 2 && identical(q$primary_source_support, c(TRUE, FALSE)))
shards <- copy(d); shards[, `:=`(exposure = "ppBMI.chr1", variant_id = "SHARD_NOT_FOR_CONTEXT")]
nb_write(rbind(d, shards), path)
filtered <- gxe_read(spec, NULL, "chr1")
check("merged ppBMI rows exclude chromosome-labelled duplicates before validation",
  nrow(filtered) == 2L && identical(filtered$interaction_beta, d$b_gi))
nb_write(shards, path)
check("chromosome-labelled rows cannot substitute for merged ppBMI",
  fails(gxe_read(spec, NULL, "chr1")))
bad <- copy(d); bad[, b_gi := NULL]; nb_write(bad, path)
check("marginal effect cannot substitute interaction", fails(gxe_read(spec, NULL, "chr1")))
bad <- copy(d); bad[, pos := c(31, 500)]; nb_write(bad, path)
check("supplied coordinates validated", fails(gxe_read(spec, NULL, "chr1")))
nb_write(d, path)
rr <- data.table(chrom = "chr1", start = c(10, 30, 1000), end = c(12, 35, 1001), is_te = TRUE)
tx <- data.table(transcript_id = c("T1", "T2", "T3"), gene_id = "G", chrom = "chr1",
  tss = c(10, 100, 300), promoter_start = c(1, 50, 250), promoter_end = c(50, 150, 350),
  promoter_te_label = c("TE", "nonTE", "nonTE"), te_tss_window_bp = 3000, tss_te_label = c("TE", "nonTE", "nonTE"))
links <- data.table(event_id = "T1", kind = "lead", chrom = "chr1", pos = 30)
unit <- copy(tx); unit[, event_id := transcript_id]
ll <- reg_context_variant_links(links, unit, rr, "chr1")
check("same insertion requires physical variant overlap", !ll$variant_overlaps_TSS_TE &&
  ll$variant_overlaps_promoter_TE && ll$variant_overlaps_proximal_only_TE)
links[, pos := 40]
ll <- reg_context_variant_links(links, unit, rr, "chr1")
check("window position alone is insufficient", !ll$variant_overlaps_promoter_TE)
meta <- data.table(modality = "isoforms", phenotype_id = c("G__T1", "G__T2"),
  exposure = "ppBMI", meta_variant_id = NA_character_, beta_meta = NA_real_, se_meta = NA_real_,
  scope = c("shared_discordant_lead", "ancestry_only"), lead_concordant = FALSE,
  n_ancestries = c(2, 1), p_feature_acat = c(.01, .5), q_feature_acat = c(.04, NA_real_))
mp <- file.path(home, "gxe_tier1_meta.tsv.gz"); nb_write(meta, mp)
ms <- copy(spec); ms[, `:=`(tier = "tier1_meta", path = mp)]
mq <- gxe_read(ms, NULL, "chr1")
check("ancestry-only synthesis retained", nrow(mq) == 2 && all(is.na(mq$variant_id)))
meta[1, beta_meta := .5]; nb_write(meta, mp)
check("discordant leads cannot have common variant effect", fails(gxe_read(ms, NULL, "chr1")))
meta[1, beta_meta := NA_real_]; nb_write(meta, mp)
tt <- data.table(modality = "isoforms", phenotype_id = "G__T1", exposure = "ppBMI",
  variant_id = "chr1:30:A:G", beta_meta = .2, se_meta = .1, p_meta = .01, q_meta = .04)
tp <- file.path(home, "gxe_tier2_results.tsv.gz"); nb_write(tt, tp)
ts <- copy(spec); ts[, `:=`(tier = "tier2", path = tp)]
tq <- gxe_read(ts, NULL, "chr1")
check("tier2 remains exploratory", is.na(tq$primary_source_support) && tq$exploratory_source_flag &&
  grepl("uncalibrated", tq$statistical_scope))
sp <- file.path(home, "inputs.tsv"); nb_write(spec, sp)
writeLines("producing ppBMI configuration declaration", file.path(home, "producer.tsv"))
writeLines("producing runtime provenance declaration", file.path(home, "producer.txt"))
ep <- file.path(home, "exposures.tsv")
e <- data.table(exposure = "ppBMI", exposure_type = "continuous", transformation = "none",
  scale_definition = "pooled across ancestries", range_definition = "declared producing range",
  producing_config = file.path(home, "producer.tsv"), producing_provenance = file.path(home, "producer.txt"),
  enabled = TRUE)
nb_write(e, ep)
map <- data.table(transcript_id = tx$transcript_id, event_id = paste0("E", 1:3), gene_id = "G")
pm <- file.path(home, "map.tsv")
nb_write(data.table(modality = "isoforms", phenotype_id = d$phenotype_id, gene_id = "G",
  event_id = c("E1", "E2"), mapping_level = "event", mapping_status = "verified"), pm)
z <- list(human_context = list(gxe_tables = sp, exposure_manifest = ep, phenotype_map = pm))
check("manifest and producing declarations validated", nrow(gxe_specs(cfg, z)$tables) == 1)
e[, enabled := FALSE]; nb_write(e, ep)
check("disabled ppBMI blocked", fails(gxe_specs(cfg, z)))
e[, enabled := TRUE]; nb_write(e, ep)
s <- list(output_dir = file.path(home, "out"))
files <- gxe_context(cfg, z, tx, map, rr, "chr1", NULL, s)
ctx <- reg_read(files[1])
check("tested nonsupport differs from untested", identical(ctx$exact_QTL, c(TRUE, FALSE, NA)))
check("proximal same-insertion support separated", ctx$proximal_same_TE_lead[1] && !ctx$same_TE_lead[1])
check("context retains separate source strata", all(ctx$ancestry == "EAS") && all(ctx$tier == "tier1"))
cfg$cohort_id <- "GUSTO"; cfg$profile <- "human_trio"
results <- data.table(feature_id = c("T1", "T2"), gene_id = "G", level = "transcript",
  status = "eligible", module = "h6", cohort_id = "GUSTO", predictor = "target_excluded",
  exposure = "matching_TE", discovery = FALSE, primary_p = c(.2, .8))
follow <- ex_h5_gxe_followup(results, ctx, cfg)
check("context exposure stays separate from modeled exposure", !anyDuplicated(names(follow)) &&
  all(follow$context_exposure == "ppBMI") && all(follow$exposure == "matching_TE"))
check("sparse follow-up remains explicit", all(follow$status == "sparse_or_constant_annotation"))
unlink(home, recursive = TRUE)
cat("GxE contract checks: ", checks, "\n", sep = "")
