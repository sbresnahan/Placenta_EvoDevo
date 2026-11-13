reg_empty_repeats <- function() data.table(chrom = character(), start = numeric(),
  end = numeric(), repeat_strand = character(), repeat_id = character(),
  repeat_name = character(), repeat_class = character(), repeat_family = character(),
  is_te = logical(), is_unknown = logical(), source = character())
reg_te_class <- function(cl) {
  base <- sub("/.*", "", cl)
  list(te = base %in% c("DNA", "LINE", "SINE", "LTR", "RC", "Retroposon",
    "PLE", "Penelope", "DIRS", "Crypton", "Helitron", "Maverick"),
    unknown = is.na(cl) | base %in% c("", "Unknown", "Unclassified", "Other"))
}
reg_rmsk <- function(path) {
  x <- fread(cmd = paste("gzip -dc", shQuote(path)), header = FALSE,
    select = c(6L, 7L, 8L, 10L, 11L, 12L, 13L, 17L),
    col.names = c("chrom", "start0", "end", "repeat_strand", "repeat_name",
      "repeat_class", "repeat_family", "repeat_id"), showProgress = FALSE)
  x[, start := as.numeric(start0) + 1]
  x[, c("start0") := NULL]
  cl <- reg_te_class(x$repeat_class)
  x[, `:=`(is_te = cl$te, is_unknown = cl$unknown, source = "UCSC_rmsk")]
  reg_intervals(x, "UCSC rmsk")
}
reg_gff_attr <- function(x, key) {
  pat <- paste0("(^|;)", key, "=([^;]*)")
  hit <- regexec(pat, x); v <- regmatches(x, hit)
  vapply(v, function(a) if (length(a) >= 3L) utils::URLdecode(a[3]) else NA_character_, "")
}
reg_maize_gff <- function(path) {
  con <- if (endsWith(path, ".gz")) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con)); blocks <- list()
  repeat {
    lines <- readLines(con, n = 50000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    if (!length(lines)) next
    x <- fread(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE,
      quote = "", fill = FALSE, showProgress = FALSE)
    if (ncol(x) != 9L) stop("TE GFF3 must have exactly nine columns")
    setnames(x, c("chrom", "method", "type", "start", "end", "score",
      "repeat_strand", "phase", "attributes"))
    x[, `:=`(repeat_id = reg_gff_attr(attributes, "ID"),
      parent = reg_gff_attr(attributes, "Parent"),
      repeat_name = reg_gff_attr(attributes, "Name"),
      repeat_class = reg_gff_attr(attributes, "Classification"))]
    # Native v4: top-level LTR_retrotransposon/helitron/etc; children are
    # parts of the same insertion, not additional independent insertions.
    x <- x[is.na(parent)]
    x[is.na(repeat_class) & grepl("LTR", type, ignore.case = TRUE), repeat_class := "LTR"]
    x[is.na(repeat_class) & grepl("helitron", type, ignore.case = TRUE), repeat_class := "RC/Helitron"]
    x[is.na(repeat_class) & grepl("DNA|terminal_inverted|MITE|transposon", type,
      ignore.case = TRUE), repeat_class := "DNA"]
    x[is.na(repeat_class) & grepl("LINE", type, ignore.case = TRUE), repeat_class := "LINE"]
    x[is.na(repeat_class) & grepl("SINE", type, ignore.case = TRUE), repeat_class := "SINE"]
    x[is.na(repeat_class), repeat_class := "Unknown"]
    x[, repeat_family := sub("_NA_.*|_LTRsimilarity.*", "", repeat_name)]
    x[is.na(repeat_name), repeat_name := repeat_id]
    cl <- reg_te_class(x$repeat_class)
    x[, `:=`(is_te = cl$te, is_unknown = cl$unknown, source = "MaizeGDB_v4_GFF3")]
    blocks[[length(blocks) + 1L]] <- x[, names(reg_empty_repeats()), with = FALSE]
  }
  x <- rbindlist(blocks)
  if (!nrow(x)) stop("No top-level TE insertions in GFF3")
  assert_unique(x, c("chrom", "repeat_id"), "Top-level maize TE IDs")
  reg_intervals(x)
}
reg_repeatmasker_out <- function(path) {
  con <- file(path, "rt"); on.exit(close(con)); blocks <- list()
  repeat {
    lines <- readLines(con, n = 50000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[grepl("^\\s*[0-9]+\\s+", lines)]
    if (!length(lines)) next
    words <- strsplit(trimws(lines), "\\s+")
    if (any(lengths(words) < 15L)) stop("Malformed RepeatMasker .out record")
    x <- rbindlist(lapply(words, function(a) data.table(chrom = a[5],
      start = as.numeric(a[6]), end = as.numeric(a[7]),
      repeat_strand = if (a[9] == "C") "-" else a[9],
      repeat_id = a[15], repeat_name = a[10], repeat_class = a[11],
      repeat_family = sub("^[^/]+/?", "", a[11]))))
    if (any(!x$repeat_strand %in% c("+", "-"))) stop("Invalid RepeatMasker strand")
    cl <- reg_te_class(x$repeat_class)
    x[, `:=`(is_te = cl$te, is_unknown = cl$unknown, source = "fish_de_novo_RepeatMasker")]
    blocks[[length(blocks) + 1L]] <- x
  }
  x <- rbindlist(blocks)
  if (!nrow(x)) stop("RepeatMasker output has no repeat hits; inspect tools/library before inference")
  reg_intervals(x)
}
reg_bed <- function(path) {
  con <- if (endsWith(path, ".gz")) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con)); blocks <- list()
  repeat {
    a <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(a)) break
    a <- a[nzchar(a) & !grepl("^(#|track |browser )", a)]
    if (!length(a)) next
    d <- fread(text = paste(a, collapse = "\n"), header = FALSE, sep = "\t",
      fill = TRUE, quote = "", select = 1:3, showProgress = FALSE)
    if (ncol(d) != 3L) stop("BED needs >=3 tab-separated columns")
    setnames(d, c("chrom", "start", "end"))
    if (!is.numeric(d$start) || !is.numeric(d$end)) stop("Nonnumeric BED coordinates")
    d[, start := start + 1]
    blocks[[length(blocks) + 1L]] <- d
  }
  if (!length(blocks)) stop("Empty BED input; missing data is not a negative result")
  reg_intervals(rbindlist(blocks), "BED")
}
reg_resource_defaults <- function(cfg, z) {
  row <- function(id, url, format, build = z$build, md5_url = NA_character_) {
    data.table(resource_id = id, url = url, path = file.path(z$cache_dir, basename(url)),
      format = format, build = build, version = basename(dirname(url)),
      md5_url = md5_url, required = TRUE, selection_role = "annotation",
      scope = "full_reference")
  }
  if (cfg$species %in% c("human", "mouse")) {
    build <- if (cfg$species == "human") "hg38" else "mm10"
    valid <- if (build == "hg38") c("hg38", "GRCh38") else c("mm10", "GRCm38")
    if (!z$build %in% valid) stop("No default repeat URL for build ", z$build)
    base <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/", build, "/database/")
    d <- rbind(row("repeats", paste0(base, "rmsk.txt.gz"), "rmsk"),
      row("rmsk_schema", paste0(base, "rmsk.sql"), "sql"))
    if (cfg$species == "human") {
      base <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200761/suppl/"
      for (mark in c("H3K27ac", "H3K4me1", "H3K4me3"))
        d <- rbind(d, row(paste0("hTSC_", mark),
          paste0(base, "GSE200761_hTSC_CutTag_", mark, "_peaks.bed.gz"), "bed"))
    }
  } else if (cfg$species == "maize") {
    if (!z$build %in% c("B73_RefGen_v4", "Zm-B73-REFERENCE-GRAMENE-4.0"))
      stop("Native maize defaults require B73v4")
    d <- row("repeats", paste0("https://download.maizegdb.org/",
      "Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0.TEs.1.0.gff3.gz"), "gff3")
  } else {
    if (z$build != "GCF_027474105.1") stop("Fish default requires GCF_027474105.1")
    b <- paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/474/105/",
      "GCF_027474105.1_Auburn_PproV2/")
    d <- rbind(row("fish_genome", paste0(b,
      "GCF_027474105.1_Auburn_PproV2_genomic.fna.gz"), "fasta_gz",
      md5_url = paste0(b, "md5checksums.txt")),
      row("fish_report", paste0(b,
        "GCF_027474105.1_Auburn_PproV2_assembly_report.txt"), "text",
        md5_url = paste0(b, "md5checksums.txt")))
  }
  # Explicit manifest replaces defaults by resource_id, allowing local-first
  # reuse without silently falling back when a requested local file is absent.
  if (!is.null(z$resources)) {
    local <- reg_read(reg_path(z$resources, cfg))
    require_columns(local, c("resource_id", "path", "format", "build"))
    if ("repeats" %in% local$resource_id &&
        (!"scope" %in% names(local) || anyNA(local[resource_id == "repeats"]$scope)))
      stop("Override repeat resources require scope=full_reference or subset")
    local[, path := vapply(path, reg_path, "", cfg = cfg)]
    d <- rbind(d[!resource_id %in% local$resource_id], local, fill = TRUE)
  }
  if (!is.null(z$repeat_file)) {
    d <- d[resource_id != "repeats"]
    d <- rbind(d, data.table(resource_id = "repeats", path = z$repeat_file,
      format = reg_default(z$repeat_format, "rmsk"), build = z$build,
      url = NA_character_, required = TRUE,
      scope = reg_default(z$repeat_scope, "subset")), fill = TRUE)
  }
  assert_unique(d, "resource_id", "Regulatory resource IDs")
  assert_unique(d, "path", "Regulatory resource destinations")
  for (nm in c("url", "md5_url")) if (!nm %in% names(d)) d[, (nm) := NA_character_]
  if (any(!d[resource_id == "repeats"]$scope %in% c("full_reference", "subset")))
    stop("Repeat scope must be full_reference or subset")
  if (anyNA(d[, .(path, format, build)])) stop("Incomplete resource manifest")
  if (any(d$format %in% c("rmsk", "gff3", "rmout", "bed") & d$build != z$build))
    stop("Resource/reference build mismatch")
  d
}
reg_prepare_resources <- function(cfg, z, root, job) {
  d <- reg_resource_defaults(cfg, z)
  dir.create(job, recursive = TRUE, showWarnings = FALSE)
  reg_write(d, file.path(job, "resources.tsv"))
  commands <- c("#!/bin/bash", "set -eo pipefail", "umask 077")
  for (i in seq_len(nrow(d))) {
    if (is.na(d$url[i]) || !nzchar(d$url[i])) {
      if (!file.exists(d$path[i])) stop("Configured local resource is missing: ", d$path[i])
      next
    }
    commands <- c(commands, paste("bash", shQuote(file.path(root, "bin/fetch_regulatory.sh")),
      shQuote(d$url[i]), shQuote(d$path[i]), shQuote(d$format[i]),
      shQuote(if (is.na(d$md5_url[i])) "-" else d$md5_url[i])))
  }
  commands <- c(commands, paste(shQuote(file.path(root, "bin/Rscript_sif")),
    shQuote(file.path(root, "poe_regulatory.R")), shQuote(attr(cfg, "file")),
    "worker-resources", "--manifest", shQuote(file.path(job, "resources.tsv"))))
  writeLines(commands, file.path(job, "run.sh"))
  invisible(d)
}
reg_finalize_resources <- function(cfg, z, root, manifest) {
  d <- reg_read(manifest)
  reg_resource_defaults(cfg, z) # revalidate declared builds
  for (i in seq_len(nrow(d))) {
    if (!file.exists(d$path[i]) || file.info(d$path[i])$size <= 0)
      stop("Missing/empty acquired resource: ", d$resource_id[i])
    if (d$format[i] == "bed") reg_bed(d$path[i])
    if (d$format[i] == "sql" &&
        !any(grepl("genoStart", readLines(d$path[i], warn = FALSE))))
      stop("Invalid rmsk SQL schema")
  }
  # The R image need not see /risapps. Runtime/image checksums are host receipts.
  dest <- file.path(z$output_dir, "resource_manifest.tsv")
  norm <- file.path(z$output_dir, "repeat_intervals.tsv.gz")
  outputs <- if (cfg$species == "fish" && !"repeats" %in% d$resource_id) dest else c(dest, norm)
  if (!reg_begin(z, "regulatory-resources", outputs)) return(invisible(TRUE))
  d[, `:=`(bytes = file.info(path)$size,
    sha256 = vapply(path, reg_sha, ""), retrieved_utc = format(Sys.time(), tz = "UTC"),
    source_status = ifelse(is.na(url), "configured_local", "local_or_verified_download"))]
  if ("repeats" %in% d$resource_id) {
    r <- d[resource_id == "repeats"]
    f <- switch(r$format, rmsk = reg_rmsk, gff3 = reg_maize_gff,
      rmout = reg_repeatmasker_out, stop("Unknown repeat format"))
    x <- f(r$path); reg_write(x, norm)
  }
  reg_write(d, dest)
  reg_finish(z, "regulatory-resources", outputs)
}
