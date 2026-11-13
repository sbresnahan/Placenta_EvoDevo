reg_fish_cpus <- function(z) {
  n <- z$fish$cpus
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 4L || n != floor(n))
    stop("Set regulatory.fish.cpus explicitly to allocated slots (>=4 for RMBlast)")
  as.integer(n)
}
reg_fish_paths <- function(z) {
  d <- file.path(z$work_dir, "fish_repeats")
  list(work = d, genome = file.path(z$cache_dir, "fish_genome.fa"),
    library = file.path(d, "Pprolifica-families.fa"),
    stockholm = file.path(d, "Pprolifica-families.stk"),
    log = file.path(d, "Pprolifica-rmod.log"),
    mask = file.path(d, "mask"),
    out = file.path(d, "mask", "fish_genome.fa.out"),
    tbl = file.path(d, "mask", "fish_genome.fa.tbl"),
    gff = file.path(d, "mask", "fish_genome.fa.out.gff"),
    intervals = file.path(z$output_dir, "fish_repeat_intervals.tsv.gz"),
    provenance = file.path(z$output_dir, "fish_repeat_provenance.tsv"))
}
reg_fasta_signature <- function(path) {
  con <- if (endsWith(path, ".gz")) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con)); name <- NULL; parts <- character(); ans <- list()
  flush <- function() {
    if (is.null(name)) return()
    seq <- toupper(paste(parts, collapse = ""))
    if (!nzchar(seq) || grepl("[^ACGTRYSWKMBDHVN]", seq))
      stop("Empty/invalid genomic FASTA sequence: ", name)
    ans[[length(ans) + 1L]] <<- data.table(chrom = name, length = nchar(seq),
      sequence_sha256 = digest::digest(seq, serialize = FALSE, algo = "sha256"))
  }
  repeat {
    lines <- readLines(con, n = 50000L, warn = FALSE)
    if (!length(lines)) break
    heads <- which(startsWith(lines, ">"))
    pos <- 1L
    for (i in heads) {
      if (i > pos) parts <- c(parts, lines[pos:(i - 1L)])
      flush()
      name <- sub("[ \t].*", "", substring(lines[i], 2L))
      parts <- character(); pos <- i + 1L
    }
    if (pos <= length(lines)) parts <- c(parts, lines[pos:length(lines)])
  }
  flush()
  x <- rbindlist(ans)
  if (!nrow(x)) stop("No sequences in genomic FASTA")
  assert_unique(x, "chrom", "FASTA contig IDs"); x
}
reg_fish_reference <- function(cfg, z) {
  p <- reg_fish_paths(z)
  d <- reg_read(file.path(z$output_dir, "resource_manifest.tsv"))
  d <- d[resource_id %in% c("fish_genome", "fish_report")]
  if (!all(c("fish_genome", "fish_report") %in% d$resource_id))
    stop("Run regulatory-resources for the fish genome/report first")
  dir.create(z$cache_dir, recursive = TRUE, showWarnings = FALSE)
  # The copied FASTA is never the frozen primary genome.
  if (!file.exists(p$genome)) {
    tmp <- paste0(p$genome, ".partial")
    reg_run(cfg$tools$gzip, c("-dc", d[resource_id == "fish_genome"]$path),
      stdout = tmp)
    if (!file.rename(tmp, p$genome)) stop("Cannot finalize fish FASTA")
  }
  official <- reg_fasta_signature(d[resource_id == "fish_genome"]$path)
  staged <- reg_fasta_signature(p$genome)
  if (!identical(official, staged)) stop("Staged fish FASTA differs from NCBI reference")
  if (sum(staged$length) != 674152735 || nrow(staged) != 415L)
    stop("Unexpected GCF_027474105.1 reference lengths")
  if (file.exists(cfg$paths$genome)) {
    primary <- reg_fasta_signature(cfg$paths$genome)
    a <- if (!is.null(z$contig_aliases)) reg_read(z$contig_aliases) else NULL
    primary[, chrom := reg_alias(chrom, official$chrom, a)]
    setorder(primary, chrom); setorder(official, chrom)
    if (!identical(primary, official))
      stop("Frozen primary fish genome is not sequence-identical to discovery reference")
  }
  invisible(staged)
}
reg_fish_runtime <- function(cfg, z, tool, args) {
    c(z$fish$executable, "run", "--cleanenv",
      "--bind", paste0(z$work_dir, ":", z$work_dir),
      "--bind", paste0(z$cache_dir, ":", z$cache_dir),
      "--env", "BLAST_USAGE_REPORT=false,OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1",
      z$fish$image, tool, args)
}
reg_command <- function(words) paste(vapply(words, shQuote, ""), collapse = " ")
reg_fish_begin <- function(cfg, z, root, phase, resume = NULL) {
  if (cfg$species != "fish") stop("fish-repeats requires a fish profile")
  reg_fish_cpus(z)
  if (!reg_valid(z, "regulatory-resources")) stop("Validated fish resources required")
  p <- reg_fish_paths(z)
  outputs <- if (phase == "discovery") c(p$library, p$stockholm, p$log) else
    c(p$out, p$tbl, p$gff, p$intervals, p$provenance)
  stage <- paste0("fish-", phase)
  if (reg_valid(z, stage)) return(FALSE)
  if (phase == "masking") {
    if (!reg_valid(z, "fish-discovery"))
      stop("Complete fish discovery is required before masking")
    if (!file.exists(p$genome) || file.info(p$genome)$size <= 0)
      stop("Missing staged fish reference for masking")
  } else {
    # One-time build/identity validation on discovery, not after every stage.
    reg_fish_reference(cfg, z)
  }
  if (!is.null(resume)) {
    if (phase != "discovery") stop("Recovery applies only to discovery")
    if (!dir.exists(resume) ||
        dirname(normalizePath(resume)) != normalizePath(p$work) ||
        !grepl("^RM_", basename(normalizePath(resume))))
      stop("Recovery directory must be an RM_* directory in this run's work directory")
    old <- reg_receipt_path(z, stage)
    if (!file.exists(old) || !identical(readRDS(old)$state, "running"))
      stop("Recovery requires a previous running receipt")
  } else if (phase == "discovery" && dir.exists(p$work) &&
      length(list.files(p$work, all.files = FALSE)))
    stop("Fish discovery work exists; choose a new work_dir or explicit recovery directory")
  dir.create(p$work, recursive = TRUE, showWarnings = FALSE)
  if (phase == "masking") dir.create(p$mask, recursive = TRUE, showWarnings = FALSE)
  reg_begin(z, stage, outputs, recover = !is.null(resume))
}
reg_fish_commands <- function(cfg, z, phase, resume = NULL) {
  p <- reg_fish_paths(z)
  reg_fish_cpus(z)
  base <- c("#!/bin/bash", "set -eo pipefail", "umask 077",
    "export BLAST_USAGE_REPORT=false", paste("cd", shQuote(p$work)))
  if (phase == "discovery") {
    if (is.null(resume)) base <- c(base, reg_command(reg_fish_runtime(cfg, z,
      "BuildDatabase", c("-name", "Pprolifica", p$genome))))
    args <- c("-database", "Pprolifica", "-threads", z$fish$cpus,
      "-srand", z$fish$seed, "-LTRStruct")
    if (!is.null(resume)) args <- c(args, "-recoverDir", resume)
    base <- c(base, paste(reg_command(reg_fish_runtime(cfg, z, "RepeatModeler", args)),
      "> discovery.stdout.log 2> discovery.stderr.log"))
  } else {
    pa <- max(1L, z$fish$cpus %/% 4L) # RMBlast uses four cores per parallel worker
    base <- c(base, paste(reg_command(reg_fish_runtime(cfg, z, "RepeatMasker",
      c("-engine", "rmblast", "-pa", pa, "-lib", p$library,
        "-gff", "-xsmall", "-dir", p$mask, p$genome))),
      "> masking.stdout.log 2> masking.stderr.log"))
  }
  base
}
reg_fish_finish <- function(cfg, z, root, phase) {
  p <- reg_fish_paths(z)
  if (phase == "discovery") {
    for (f in c(p$library, p$stockholm, p$log)) reg_sha(f)
    heads <- readLines(p$library, warn = FALSE)
    if (!any(startsWith(heads, ">"))) stop("Discovery did not produce a FASTA family library")
    if (!any(grepl("Random Number Seed|seed|srand", readLines(p$log, warn = FALSE),
        ignore.case = TRUE))) warning("Seed not found in summarized log; retained in run command")
    reg_finish(z, "fish-discovery", c(p$library, p$stockholm, p$log))
  } else {
    for (f in c(p$out, p$tbl, p$gff)) reg_sha(f)
    x <- reg_repeatmasker_out(p$out)
    reg_write(x, p$intervals)
    reg_write(data.table(build = z$build, genome_sha256 = reg_sha(p$genome),
      family_library_sha256 = reg_sha(p$library),
      image = z$fish$image,
      seed = z$fish$seed, cpus = z$fish$cpus,
      rmblast_parallel_workers = z$fish$cpus %/% 4L,
      n_repeat_hits = nrow(x), classified_te_hits = sum(x$is_te),
      unknown_hits = sum(x$is_unknown),
      limitation = "reference_insertions_not_Line_B_PAV_or_TE_expression"), p$provenance)
    reg_finish(z, "fish-masking",
      c(p$out, p$tbl, p$gff, p$intervals, p$provenance))
  }
}
