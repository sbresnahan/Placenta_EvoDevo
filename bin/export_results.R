#!/usr/bin/env Rscript
# Results archive for downstream reporting: per-cohort primary results, QC,
# additional_v1 analyses, hypothesis #1 (TE-TSS POE: within-gene te-tss and
# between-gene te-genes) and hypothesis #5
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
TE_GENES <- c("te_genes_result.tsv", "te_genes_balance.tsv", "te_genes_gene_table.tsv",
  "te_genes_placebo.tsv", "te_genes_provenance.tsv")
HUMAN_CONTEXT <- c("human_regulatory_context.tsv", "human_qtl_links.tsv.gz",
  "human_peak_overlaps.tsv.gz", "human_mapping_coverage.tsv",
  "human_resource_coverage.tsv", "human_context_provenance.tsv")
DOCS <- c(README.md = "README.md", REGULATORY_MODULES.md = "docs/REGULATORY_MODULES.md",
  REVISED_RUNBOOK.md = "docs/REVISED_RUNBOOK.md", SCHEMA.md = "config/SCHEMA.md",
  evidence_review.md = "evidence/review.md", dataset_evidence.csv = "evidence/dataset_evidence.csv",
  ATTRIBUTION.md = "ATTRIBUTION.md")

named <- function(paths) stats::setNames(paths, basename(paths))
read_or_null <- function(path) {
  if (is.null(path) || !file.exists(path) || file.info(path)$size <= 0) return(NULL)
  tryCatch(fread(path), error = function(e) NULL)
}
# Analysis settings without paths or private locations (values or keys that point at files).
SETTINGS_SECTIONS <- c("reference", "parents", "dosage", "markers", "ornaments", "pilot",
  "fragments", "models", "additional_analysis", "phase_origin", "gusto_precomputed",
  "regulatory", "fish_stage_policy")
# Whole-word key segments only, so e.g. min_per_direction is kept but output_dir is not.
PATHLIKE_KEY <- "(^|[._])(paths?|files?|dir|vcf|matrix|map|tables|beds|chain|source|image|executable|resources|aliases|fastq|sif|rmsk)($|[._])"
flatten_settings <- function(x, prefix) {
  out <- character()
  for (k in names(x)) {
    v <- x[[k]]; key <- paste0(prefix, k)
    if (is.list(v)) out <- c(out, flatten_settings(v, paste0(key, ".")))
    else out[key] <- paste(as.character(v), collapse = ";")
  }
  out
}
cohort_settings <- function(cfg) {
  v <- unlist(lapply(SETTINGS_SECTIONS, function(k) {
    x <- cfg[[k]]
    if (is.null(x)) character() else if (is.list(x)) flatten_settings(x, paste0(k, ".")) else
      stats::setNames(paste(as.character(x), collapse = ";"), k)
  }))
  if (!length(v)) return(data.table(setting = character(), value = character()))
  keep <- !grepl("^/|/rsrch|\\$\\{", v) & !grepl(PATHLIKE_KEY, names(v), ignore.case = TRUE)
  data.table(setting = names(v)[keep], value = unname(v[keep]))
}
# Per-metric n/median/IQR, by pilot mode/strand where present. For human cohorts,
# cells with fewer than min_n samples report n only.
qc_summary <- function(path, human, min_n = 5L) {
  d <- read_or_null(path)
  if (is.null(d) || !nrow(d)) return(NULL)
  grp <- intersect(c("mode", "strand"), names(d))
  num <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], ID_COLUMNS)
  if (!length(num)) return(NULL)
  long <- rbindlist(lapply(num, function(v) {
    x <- if (length(grp)) d[, grp, with = FALSE] else data.table(group = rep("all", nrow(d)))
    x[, `:=`(metric = v, value = as.numeric(d[[v]]))]
    x
  }))
  by <- c(if (length(grp)) grp else "group", "metric")
  out <- long[is.finite(value), .(n = .N, median = stats::median(value),
    q25 = stats::quantile(value, .25, names = FALSE),
    q75 = stats::quantile(value, .75, names = FALSE)), by = by]
  if (human) out[n < min_n, `:=`(median = NA_real_, q25 = NA_real_, q75 = NA_real_)]
  out
}
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
    manifests = list(applies = TRUE, withheld = human,
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
    h1b_te_genes = list(applies = TRUE, ok = reg_valid(z, "te-genes"),
      files = named(file.path(z$output_dir, TE_GENES)), extra = character(),
      why = "te-genes not run or incomplete (needs te-tss; GUSTO also needs additional_v1)"),
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
  kv$models_min_orientation <- cfg$models$min_orientation
  covs <- unlist(cfg$models$covariates)
  kv$models_covariates <- if (length(covs)) paste(covs, collapse = ";") else "none"
  # Design counts are aggregate for every cohort; unit-level detail stays public-only.
  u <- read_or_null(cfg$paths$units); l <- read_or_null(cfg$paths$libraries)
  r <- read_or_null(cfg$paths$sources)
  if (!is.null(u)) {
    kv$n_units <- nrow(u)
    for (k in c("tissue", "stage")) if (k %in% names(u))
      kv[[k]] <- {
        v <- unique(stats::na.omit(as.character(u[[k]]))); v <- v[nzchar(v)]
        if (length(v)) paste(v, collapse = ";") else "not_recorded_in_units_manifest"
      }
    if (cfg$profile == "public_cross" && "maternal_parent" %in% names(u))
      kv$units_by_maternal_parent <- paste(names(table(u$maternal_parent)),
        as.integer(table(u$maternal_parent)), sep = "=", collapse = ";")
  }
  if (!is.null(l)) kv$n_libraries <- nrow(l)
  if (!is.null(r)) kv$n_runs <- nrow(r)
  if (cfg$profile == "human_trio" && !is.null(l) && "child_id" %in% names(l)) {
    t <- read_or_null(cfg$paths$trios)
    if (!is.null(t) && all(c("child_id", "family_id") %in% names(t)))
      kv$n_families <- uniqueN(t$family_id[t$child_id %in% l$child_id])
  }
  kv <- kv[!vapply(kv, function(v) is.null(v) || !length(v), logical(1))]
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
  bad <- Filter(function(part) part$applies && !isTRUE(part$withheld) && !part$ok, x$parts)
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
    write_tsv(cohort_settings(cfg), file.path(top, base, "settings.tsv"))
    for (f in c("cohort.tsv", "settings.tsv")) {
      dest <- file.path(top, base, f)
      manifest[[length(manifest) + 1L]] <- data.table(cohort_id = cfg$cohort_id, part = "design",
        archive_path = file.path(base, f), source_path = paste0("derived:", attr(cfg, "file")),
        state = "included", bytes = file.info(dest)$size,
        sha256 = digest::digest(file = dest, algo = "sha256"))
    }
    for (part in names(x$parts)) {
      s <- x$parts[[part]]
      if (!s$applies) {
        availability[[length(availability) + 1L]] <- data.table(cohort_id = cfg$cohort_id,
          part = part, state = "not_applicable", detail = "")
        next
      }
      if (isTRUE(s$withheld)) {
        availability[[length(availability) + 1L]] <- data.table(cohort_id = cfg$cohort_id,
          part = part, state = "withheld_individual_level",
          detail = "human manifests list participant and sample identifiers; aggregate counts are in cohort.tsv")
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
      if (part == "qc") {
        for (k in c("pilot_report", "metrics")) {
          sm <- qc_summary(cfg$paths[[k]], human)
          if (is.null(sm)) next
          rel <- file.path(base, part, paste0(if (k == "metrics") "quantification_metrics" else k,
            "_summary.tsv"))
          dir.create(dirname(file.path(top, rel)), recursive = TRUE, showWarnings = FALSE)
          write_tsv(sm, file.path(top, rel))
          rows <- rbind(rows, data.table(cohort_id = cfg$cohort_id, part = part, archive_path = rel,
            source_path = paste0("derived_summary:", cfg$paths[[k]]), state = "included",
            bytes = file.info(file.path(top, rel))$size,
            sha256 = digest::digest(file = file.path(top, rel), algo = "sha256")), fill = TRUE)
        }
      }
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

  # Seadragon's scripts folder is not a git checkout: fingerprint the code itself.
code <- c(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE),
  list.files(file.path(root, "bin"), full.names = TRUE, recursive = TRUE),
  file.path(root, c("poe.R", "run_pipeline.R", "runtime.R", "poe_regulatory.R")))
code <- sort(code[file.exists(code) & !dir.exists(code) & !grepl("__pycache__", code)])
scripts <- data.table(path = substring(code, nchar(root) + 2L),
  sha256 = vapply(code, function(f) digest::digest(file = f, algo = "sha256"), ""))
write_tsv(scripts, file.path(top, "SCRIPTS.tsv"))
fingerprint <- digest::digest(paste(scripts$path, scripts$sha256, collapse = "\n"),
  algo = "sha256", serialize = FALSE)
commit <- tryCatch(system2("git", c("-C", shQuote(root), "rev-parse", "--short", "HEAD"),
    stdout = TRUE, stderr = FALSE), error = function(e) "unknown", warning = function(w) "unknown")
  if (!length(commit) || !nzchar(commit[1])) commit <- "unknown"
  avail_md <- c("| Cohort | Part | State | Detail |", "|---|---|---|---|",
    sprintf("| %s | %s | %s | %s |", availability$cohort_id, availability$part,
      availability$state, availability$detail))
  writeLines(c(
    "# Placenta_EvoDevo results archive", "",
    sprintf("Created %s from scripts commit `%s`, code fingerprint `%s` (SHA-256 over every pipeline file; per-file hashes in `SCRIPTS.tsv`)%s.",
      format(Sys.time(), "%Y-%m-%d %H:%M %Z"), paste(commit, collapse = ""), substr(fingerprint, 1, 16),
      if (force) " with `--force` (incomplete parts skipped)" else ""),
    "", "## What is being tested", "",
    "- **Primary parent-of-origin (POE) models** (`primary/`): dosage-adjusted maternal vs paternal allelic expression per gene (`model_gene.tsv`) and per promoter/TSS event (`model_promoter.tsv`); `model_pairs.tsv` holds the paired event inputs. Public crosses (maize, mouse, fish) use reciprocal-cross intercepts conditional on parental markers; GUSTO (human) uses trio parental-origin assignments, cis orientation and covariates. GUSTO `model_gene.tsv` is a placeholder (`deferred_unimplemented`).",
    "- **additional_v1** (`additional_v1/`): exploratory reanalysis of the frozen counts: gene-level POE (`gene_POE.tsv`), TSS-vs-first-donor separation (`separated_TSS_donor_POE.tsv`), and TE annotation. See `docs/REVISED_RUNBOOK.md` section 0 (Q1-Q4).",
    "- **Hypothesis #1** (`h1_te_tss/`): within a gene, conditional on the same first splice donor, do TSSs inside transposable elements show a different POE (signed difference beta_TE - beta_nonTE) or a larger POE magnitude (|beta_TE| - |beta_nonTE|) than non-TE TSSs? Start with `te_tss_cohort_summary.tsv`, then `te_tss_gene_summary.tsv`; `TE_annotation_coverage.tsv` shows how many TSSs could be classified. Magnitude inference is exploratory. See `docs/REGULATORY_MODULES.md`, \"#1\".",
  "- **Hypothesis #1, between genes** (`h1b_te_genes/`): do genes with at least one TE-derived TSS show a different gene-level POE than comparable genes whose TSSs are all non-TE? Start with `te_genes_result.tsv` (magnitude = difference in POE^2 - SE^2, which removes noise inflation; signed = difference in POE), then `te_genes_balance.tsv` (standardized mean differences before/after propensity weighting; |SMD| > 0.1 means residual imbalance). Confounders adjusted: depth, informative loci, units, SE, number of TSS events and transcripts, gene length, local TE density, promoter GC. p-values: genomic block bootstrap, within-stratum permutation, and a placebo that shifts the TE track (`te_genes_placebo.tsv`). Observational and complementary to the within-gene pairs, which cancel these confounders by design. See `docs/REGULATORY_MODULES.md`, \"#1b\".",
    "- **Hypothesis #5** (`h5_human_context/`, GUSTO only): do POE TSS events, particularly TE-derived ones, carry placental cis-QTL (PANTRY, EAS and EUR, nine modalities) and placental chromatin evidence (ENCODE ENCSR057ORE cCREs)? `human_regulatory_context.tsv` is one row per event; `human_mapping_coverage.tsv` and `human_resource_coverage.tsv` separate 'not measured' from 'measured, no link'. Gene-level QTL links are context, not validation of a TSS phenotype. This is evidence integration, not causal inference. See `docs/REGULATORY_MODULES.md`, \"#5\".",
    "", "## Interpretation rules", "",
    "- Analyse each cohort/species separately. Cohorts differ in tissue, developmental stage, ploidy and design; do not pool across species or treat several analyses of the same allocations as independent replication.",
    "- Maize endosperm is triploid: expected maternal:paternal dosage is 2:1 and is subtracted before POE testing (see `cohort.tsv`).",
    "- `cohort.tsv` records each cohort's design (parents, dosage, marker proxies, stage policy, model settings). Marker proxies and the fish mixed-stage, small-N design are flagged there and in the tables' uncertainty columns.",
    "- Human tables are aggregate-only. Any table with sample, family, child or parent identifiers was withheld (`withheld_individual_level` in `MANIFEST.tsv`); GUSTO `model_pairs.tsv` (per-sample inputs), per-library QC and manifests are withheld by design; their aggregate equivalents are `cohort.tsv` (design counts) and `qc/*_summary.tsv`.",
    "", "## Availability", "", avail_md,
    "", "## Layout", "",
    "- `cohorts/<cohort_id>/cohort.tsv`: design summary (species, tissue, stage, units, libraries, runs, families, dosage, markers, model thresholds, covariate names)",
"- `cohorts/<cohort_id>/settings.tsv`: all analysis settings from the cohort configuration, without file paths",
"- `cohorts/<cohort_id>/qc/*_summary.tsv`: per-metric n, median and IQR (by pilot mode/strand); for human cohorts cells with fewer than 5 samples show n only",
    "- `cohorts/<cohort_id>/{primary,qc,manifests,additional_v1,h1_te_tss,h1b_te_genes,h5_human_context}/`",
    "- `docs/`: method documentation (hypothesis definitions in `REGULATORY_MODULES.md`, additional_v1 questions in `REVISED_RUNBOOK.md`, settings in `SCHEMA.md`, dataset provenance in `evidence_*`)",
    "- `MANIFEST.tsv`: every file with source path, size and SHA-256; `AVAILABILITY.tsv`: the table above; `SCRIPTS.tsv`: SHA-256 of each pipeline code file used for this export"),
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
