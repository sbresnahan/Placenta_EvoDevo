#!/usr/bin/env Rscript
# Results archive for downstream reporting: per-cohort primary results, QC,
# additional_v1 analyses, hypothesis #1 (TE-TSS POE, te-tss) and hypothesis #5
# (human regulatory context), plus the method docs needed to interpret them.
#
# Usage (from the scripts directory, inside the R container):
#   bin/Rscript_sif bin/export_results.R --out /path/results.tar.gz [--force] CONFIG [CONFIG ...]
#
# Without --force every applicable part must be complete for every cohort.
# With --force, parts that are not complete yet (e.g. fish te-tss while repeat
# discovery runs) are skipped and recorded as unavailable in the archive.
# Human (GUSTO) tables are aggregate-only: any table carrying sample, family,
# child or parent identifiers is withheld, and no human manifest is exported.
f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]), fixed = TRUE)
root <- normalizePath(file.path(dirname(f), ".."), mustWork = TRUE)
source(file.path(root, "R/library_paths.R"))

usage <- "Usage: bin/Rscript_sif bin/export_results.R --out FILE.tar.gz [--force] CONFIG [CONFIG ...]"
args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args
args <- args[args != "--force"]
k <- match("--out", args)
if (is.na(k) || k == length(args)) stop(usage, call. = FALSE)
out <- args[k + 1L]
configs <- args[-c(k, k + 1L)]
if (!length(configs) || any(startsWith(configs, "--"))) stop(usage, call. = FALSE)
if (!grepl("\\.tar\\.gz$", out)) stop("--out must end in .tar.gz", call. = FALSE)
out <- file.path(normalizePath(dirname(out), mustWork = TRUE), basename(out))
if (file.exists(out)) stop("Refusing to overwrite existing archive: ", out, call. = FALSE)

for (m in c("config", "common", "workflow", "additional", "regulatory_common", "regulatory_fish"))
  source(file.path(root, "R", paste0(m, ".R")))

# Columns that identify individuals; tables containing them are not exported for humans.
ID_COLUMNS <- c("sample_id", "family_id", "child_id", "mother_id", "father_id",
  "run_id", "infant_id", "unit_id", "cluster_id", "library_name", "sample_accession")
TE_TSS <- c("transcript_TE_annotation.tsv", "event_TE_annotation.tsv",
  "TE_annotation_coverage.tsv", "te_tss_event_pairs.tsv", "te_tss_gene_summary.tsv",
  "te_tss_cohort_summary.tsv", "te_tss_provenance.tsv")
HUMAN_CONTEXT <- c("human_regulatory_context.tsv", "human_qtl_links.tsv.gz",
  "human_peak_overlaps.tsv.gz", "human_mapping_coverage.tsv",
  "human_resource_coverage.tsv", "human_context_provenance.tsv")
DOCS <- c(README.md = "README.md", REGULATORY_MODULES.md = "docs/REGULATORY_MODULES.md",
  REVISED_RUNBOOK.md = "docs/REVISED_RUNBOOK.md", SCHEMA.md = "config/SCHEMA.md",
  evidence_review.md = "evidence/review.md", dataset_evidence.csv = "evidence/dataset_evidence.csv",
  ATTRIBUTION.md = "ATTRIBUTION.md")

named <- function(paths) stats::setNames(paths, basename(paths))
present <- function(paths) length(paths) > 0L && all(file.exists(paths)) && all(file.info(paths)$size > 0)
header_of <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  h <- readLines(con, n = 1L, warn = FALSE)
  if (!length(h)) return(character())
  gsub('"', "", strsplit(h, "\t", fixed = TRUE)[[1]])
}

# Each part: whether it applies, whether it is complete, required files and optional extras.
cohort_parts <- function(cfg) {
  p <- cfg$paths
  human <- cfg$profile == "human_trio"
  z <- reg_settings(cfg)
  models <- unlist(p[c("model_gene", "model_promoter", "model_pairs")])
  te <- file.path(z$output_dir, TE_TSS)
  hc <- file.path(z$output_dir, HUMAN_CONTEXT)
  h1_extra <- file.path(z$output_dir, "resource_manifest.tsv")
  if (identical(cfg$species, "fish")) h1_extra <- c(h1_extra, reg_fish_paths(z)$provenance)
  h5_extra <- character()
  hcfg <- z$human_context
  if (!is.null(hcfg$qtl_tables)) {
    d <- dirname(reg_path(hcfg$qtl_tables, cfg))
    h5_extra <- c(h5_extra, qtl_input_summary.tsv = file.path(d, "summary.tsv"),
      qtl_phenotype_map.tsv = file.path(d, "phenotype_map.tsv"))
  }
  if (!is.null(hcfg$beds)) {
    d <- dirname(reg_path(hcfg$beds, cfg))
    h5_extra <- c(h5_extra, peak_input_summary.tsv = file.path(d, "summary.tsv"))
  }
  list(
    primary = list(applies = TRUE, ok = stage_valid(p$model_receipt, models),
      files = named(models), extra = character(), why = "primary models receipt incomplete"),
    qc = list(applies = TRUE, ok = present(c(p$pilot_report, p$metrics)),
      files = c(pilot_report.tsv = p$pilot_report, quantification_metrics.tsv = p$metrics),
      extra = character(), why = "pilot report or quantification metrics missing"),
    manifests = list(applies = !human,
      ok = present(unlist(p[c("sources", "libraries", "units", "marker_sources")])),
      files = named(unlist(p[c("sources", "libraries", "units", "marker_sources")])),
      extra = character(), why = "cohort manifests missing"),
    additional_v1 = list(applies = TRUE,
      ok = stage_valid(additional_receipt(cfg), additional_paths(cfg)),
      files = named(additional_paths(cfg)), extra = character(),
      why = "additional_v1 receipt incomplete"),
    h1_te_tss = list(applies = TRUE, ok = reg_valid(z, "te-tss"), files = named(te),
      extra = named(h1_extra), why = if (identical(cfg$species, "fish"))
        "te-tss incomplete (fish needs fish-repeats discovery and masking first)" else
        "te-tss incomplete"),
    h5_human_context = list(applies = human, ok = human && reg_valid(z, "human-context"),
      files = named(hc), extra = h5_extra, why = "human-context incomplete"))
}

cohort_summary <- function(cfg) {
  kv <- list(cohort_id = cfg$cohort_id, species = cfg$species, profile = cfg$profile,
    reference_build = cfg$reference$build,
    focal_parent = cfg$parents$focal, alternate_parent = cfg$parents$alternate,
    dosage_maternal = cfg$dosage$maternal, dosage_paternal = cfg$dosage$paternal,
    marker_mode = cfg$markers$mode, marker_proxy = isTRUE(cfg$markers$proxy),
    marker_proxy_reason = cfg$markers$proxy_reason, fish_stage_policy = cfg$fish_stage_policy,
    models_min_depth = cfg$models$min_depth, models_min_callable_units = cfg$models$min_callable_units,
    models_min_families = cfg$models$min_families, models_bootstrap = cfg$models$bootstrap,
    models_seed = cfg$models$seed, regulatory_bootstrap = cfg$regulatory$bootstrap,
    regulatory_min_depth = cfg$regulatory$min_depth)
  if (cfg$profile == "public_cross" && file.exists(cfg$paths$units)) {
    u <- fread(cfg$paths$units)
    kv <- c(kv, list(n_units = nrow(u), tissue = paste(unique(u$tissue), collapse = ";"),
      stage = paste(unique(u$stage), collapse = ";"),
      units_by_maternal_parent = paste(names(table(u$maternal_parent)),
        as.integer(table(u$maternal_parent)), sep = "=", collapse = ";")))
  }
  kv <- kv[!vapply(kv, is.null, logical(1))]
  # Not "key": data.table() treats a `key` argument as the table's sort key.
  data.table(setting = names(kv), value = vapply(kv, function(v) paste(v, collapse = ";"), ""))
}

# Pass 1: read every config and decide what can be exported.
cohorts <- lapply(configs, function(path) {
  cfg <- read_config_v2(path)
  list(cfg = cfg, parts = cohort_parts(cfg))
})
ids <- vapply(cohorts, function(x) x$cfg$cohort_id, "")
if (anyDuplicated(ids)) {
  dup <- ids %in% ids[duplicated(ids)]
  stop("Duplicate cohort_id across configs:\n  ", paste(sprintf("%s  (%s)", ids[dup], configs[dup]),
    collapse = "\n  "), "\nPass one config per cohort.", call. = FALSE)
}
blocked <- unlist(lapply(cohorts, function(x) {
  bad <- Filter(function(part) part$applies && !part$ok, x$parts)
  if (length(bad)) paste0(x$cfg$cohort_id, " / ", names(bad), ": ", vapply(bad, `[[`, "", "why"))
}))
if (length(blocked) && !force)
  stop("Incomplete parts (rerun with --force to skip them):\n  ",
    paste(blocked, collapse = "\n  "), call. = FALSE)

# Pass 2: stage files, then tar.
name <- sub("\\.tar\\.gz$", "", basename(out))
stage <- file.path(dirname(out), paste0(".", name, ".partial"))
if (file.exists(stage)) stop("Remove the incomplete previous attempt first: ", stage, call. = FALSE)
top <- file.path(stage, name)
dir.create(file.path(top, "docs"), recursive = TRUE)
build_archive <- function() {
  manifest <- list(); availability <- list()

  copy_in <- function(cohort, part, src, rel, human) {
    dest <- file.path(top, rel)
    row <- function(state) data.table(cohort_id = cohort, part = part, archive_path = rel,
      source_path = src, state = state,
      bytes = if (state == "included") file.info(dest)$size else NA_real_,
      sha256 = if (state == "included") digest::digest(file = dest, algo = "sha256") else NA_character_)
    if (!file.exists(src) || file.info(src)$size <= 0) return(row("missing"))
    if (human && length(intersect(header_of(src), ID_COLUMNS)))
      return(row("withheld_individual_level"))
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(src, dest, copy.date = TRUE)) stop("Copy failed: ", src)
    row("included")
  }

  for (x in cohorts) {
    cfg <- x$cfg; human <- cfg$profile == "human_trio"; base <- file.path("cohorts", cfg$cohort_id)
    dir.create(file.path(top, base), recursive = TRUE, showWarnings = FALSE)
    write_tsv(cohort_summary(cfg), file.path(top, base, "cohort.tsv"))
    for (part in names(x$parts)) {
      s <- x$parts[[part]]
      if (!s$applies) {
        availability[[length(availability) + 1L]] <- data.table(cohort_id = cfg$cohort_id,
          part = part, state = "not_applicable", detail = "")
        next
      }
      if (!s$ok) {
        availability[[length(availability) + 1L]] <- data.table(cohort_id = cfg$cohort_id,
          part = part, state = "unavailable_skipped", detail = s$why)
        next
      }
      rows <- c(
        lapply(names(s$files), function(n) copy_in(cfg$cohort_id, part, s$files[[n]],
          file.path(base, part, n), human)),
        lapply(names(s$extra), function(n) copy_in(cfg$cohort_id, part, s$extra[[n]],
          file.path(base, part, n), human)))
      rows <- rbindlist(rows)
      required <- rows[seq_along(s$files)]
      if (any(required$state == "missing") && !force)
        stop(cfg$cohort_id, " / ", part, ": receipt complete but file missing: ",
          paste(required[state == "missing"]$source_path, collapse = ", "), call. = FALSE)
      manifest[[length(manifest) + 1L]] <- rows
      states <- unique(rows$state)
      availability[[length(availability) + 1L]] <- data.table(cohort_id = cfg$cohort_id,
        part = part, state = if (identical(states, "included")) "included" else "included_partial",
        detail = paste(sprintf("%s=%d", names(table(rows$state)), as.integer(table(rows$state))),
          collapse = ";"))
    }
  }
  for (n in names(DOCS)) {
    src <- file.path(root, DOCS[[n]])
    if (file.exists(src)) {
      file.copy(src, file.path(top, "docs", n))
      manifest[[length(manifest) + 1L]] <- data.table(cohort_id = "all", part = "docs",
        archive_path = file.path("docs", n), source_path = DOCS[[n]], state = "included",
        bytes = file.info(src)$size, sha256 = digest::digest(file = src, algo = "sha256"))
    }
  }
  manifest <- rbindlist(manifest, fill = TRUE)
  availability <- rbindlist(availability)
  write_tsv(manifest, file.path(top, "MANIFEST.tsv"))
  write_tsv(availability, file.path(top, "AVAILABILITY.tsv"))

  commit <- tryCatch(system2("git", c("-C", shQuote(root), "rev-parse", "--short", "HEAD"),
    stdout = TRUE, stderr = FALSE), error = function(e) "unknown", warning = function(w) "unknown")
  if (!length(commit) || !nzchar(commit[1])) commit <- "unknown"
  avail_md <- c("| Cohort | Part | State | Detail |", "|---|---|---|---|",
    sprintf("| %s | %s | %s | %s |", availability$cohort_id, availability$part,
      availability$state, availability$detail))
  writeLines(c(
    "# Placenta_EvoDevo results archive", "",
    sprintf("Created %s from scripts commit `%s`%s.", format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
      paste(commit, collapse = ""), if (force) " with `--force` (incomplete parts skipped)" else ""),
    "", "## What is being tested", "",
    "- **Primary parent-of-origin (POE) models** (`primary/`): dosage-adjusted maternal vs paternal allelic expression per gene (`model_gene.tsv`) and per promoter/TSS event (`model_promoter.tsv`); `model_pairs.tsv` holds the paired event inputs. Public crosses (maize, mouse, fish) use reciprocal-cross intercepts conditional on parental markers; GUSTO (human) uses trio parental-origin assignments, cis orientation and covariates. GUSTO `model_gene.tsv` is a placeholder (`deferred_unimplemented`).",
    "- **additional_v1** (`additional_v1/`): exploratory reanalysis of the frozen counts: gene-level POE (`gene_POE.tsv`), TSS-vs-first-donor separation (`separated_TSS_donor_POE.tsv`), TE annotation and TE comparison. See `docs/REVISED_RUNBOOK.md` section 0 (Q1-Q4).",
    "- **Hypothesis #1** (`h1_te_tss/`): within a gene, conditional on the same first splice donor, do TSSs inside transposable elements show a different POE (signed difference beta_TE - beta_nonTE) or a larger POE magnitude (|beta_TE| - |beta_nonTE|) than non-TE TSSs? Start with `te_tss_cohort_summary.tsv`, then `te_tss_gene_summary.tsv`; `TE_annotation_coverage.tsv` shows how many TSSs could be classified. Magnitude inference is exploratory. See `docs/REGULATORY_MODULES.md`, \"#1\".",
    "- **Hypothesis #5** (`h5_human_context/`, GUSTO only): do POE TSS events, particularly TE-derived ones, carry placental cis-QTL (PANTRY, EAS and EUR, nine modalities) and placental chromatin evidence (ENCODE ENCSR057ORE cCREs)? `human_regulatory_context.tsv` is one row per event; `human_mapping_coverage.tsv` and `human_resource_coverage.tsv` separate 'not measured' from 'measured, no link'. Gene-level QTL links are context, not validation of a TSS phenotype. This is evidence integration, not causal inference. See `docs/REGULATORY_MODULES.md`, \"#5\".",
    "", "## Interpretation rules", "",
    "- Analyse each cohort/species separately. Cohorts differ in tissue, developmental stage, ploidy and design; do not pool across species or treat several analyses of the same allocations as independent replication.",
    "- Maize endosperm is triploid: expected maternal:paternal dosage is 2:1 and is subtracted before POE testing (see `cohort.tsv`).",
    "- `cohort.tsv` records each cohort's design (parents, dosage, marker proxies, stage policy, model settings). Marker proxies and the fish mixed-stage, small-N design are flagged there and in the tables' uncertainty columns.",
    "- Human tables are aggregate-only. Any table with sample, family, child or parent identifiers was withheld (`withheld_individual_level` in `MANIFEST.tsv`); GUSTO `model_pairs.tsv` (per-sample inputs) is withheld by design.",
    "", "## Availability", "", avail_md,
    "", "## Layout", "",
    "- `cohorts/<cohort_id>/cohort.tsv`: design and settings",
    "- `cohorts/<cohort_id>/{primary,qc,manifests,additional_v1,h1_te_tss,h5_human_context}/`",
    "- `docs/`: method documentation (hypothesis definitions in `REGULATORY_MODULES.md`, additional_v1 questions in `REVISED_RUNBOOK.md`, settings in `SCHEMA.md`, dataset provenance in `evidence_*`)",
    "- `MANIFEST.tsv`: every file with source path, size and SHA-256; `AVAILABILITY.tsv`: the table above"),
    file.path(top, "README.md"))

  tmp <- paste0(out, ".partial")
  owd <- setwd(stage)
  status <- tryCatch(utils::tar(tmp, files = name, compression = "gzip",
    tar = if (nzchar(Sys.which("tar"))) "tar" else "internal"), finally = setwd(owd))
  if ((length(status) && !identical(as.integer(status), 0L)) || !file.exists(tmp))
    stop("tar failed")
  if (!file.rename(tmp, out)) stop("Cannot finalize ", out)
  unlink(stage, recursive = TRUE)
  print(availability)
  cat(sprintf("\nArchive: %s (%d files, %.1f MB)\n", out, nrow(manifest[state == "included"]),
    file.info(out)$size / 1e6))
}
tryCatch(build_archive(), error = function(e) {
  unlink(c(stage, paste0(out, ".partial")), recursive = TRUE)
  stop(conditionMessage(e), call. = FALSE)
})
