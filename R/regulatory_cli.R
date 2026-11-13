reg_actions <- c("regulatory-resources", "fish-repeats", "te-tss",
  "human-context", "regulatory-status")

reg_lsf_resource <- function(cfg, action, z) {
  if (action == "fish-repeats") {
    key <- "fish_repeats"
    r <- cfg$lsf$stages[[key]]
    if (is.null(r$queue)) stop("Set lsf.stages.", key, ".queue to a site-approved queue")
    for (k in c("cpus", "memory_value", "wall_time"))
      if (is.null(r[[k]])) stop("Explicit lsf.stages.", key, ".", k, " is required")
    r$reservation_value <- reg_default(r$reservation_value, r$memory_value)
    if (r$cpus < reg_fish_cpus(z)) stop("Fish CPU allocation below regulatory.fish.cpus")
  } else {
    key <- gsub("-", "_", action)
    r <- cfg$lsf$stages[[key]]
    if (is.null(r)) return(stage_resource(cfg,
      if (action == "regulatory-resources") "retrieve" else "models"))
    r$queue <- reg_default(r$queue, cfg$lsf$queue)
    r$reservation_value <- reg_default(r$reservation_value, r$memory_value)
  }
  for (k in c("cpus", "memory_value", "reservation_value")) positive(r[[k]], k)
  if (r$cpus != floor(r$cpus)) stop("LSF CPUs must be a whole allocated slot count")
  if (!grepl("^[A-Za-z0-9_.-]+$", r$queue) ||
      !grepl("^[0-9]+:[0-5][0-9]$", r$wall_time)) stop("Unsafe regulatory LSF resource")
  r
}
reg_plan_jobs <- function(cfg, z, root, action, submit, after = "none", resume = NULL) {
  if (!grepl("^(none|[0-9]+)$", after)) stop("--after requires a numeric LSF job ID")
  if (!is.null(resume) && action != "fish-repeats") stop("--resume requires fish-repeats")
  if (action == "fish-repeats" && cfg$species != "fish") stop("Fish profile required")
  if (action == "human-context" && cfg$species != "human") stop("Human profile required")
  r <- reg_lsf_resource(cfg, action, z)
  job <- file.path(z$output_dir, "jobs",
    paste0(action, "_", format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid()))
  dir.create(job, recursive = TRUE, showWarnings = FALSE)
  if (any(grepl("[\r\n]", c(job, root, attr(cfg, "file"))))) stop("Unsafe path")
  entry <- c(file.path(root, "bin/Rscript_sif"), file.path(root, "poe_regulatory.R"),
    attr(cfg, "file"))
  if (action == "regulatory-resources") {
    reg_prepare_resources(cfg, z, root, job)
    scripts <- file.path(job, "run.sh")
  } else if (action == "fish-repeats") {
    scripts <- file.path(job, c("discovery.sh", "masking.sh"))
    for (i in seq_along(scripts)) {
      phase <- c("discovery", "masking")[i]
      res <- if (phase == "discovery") resume else NULL
      args <- c("--phase", phase, if (!is.null(res)) c("--resume", res))
      lines <- c("#!/bin/bash", "set -eo pipefail", "umask 077",
        paste0("export SCRIPTS_DIR=", shQuote(root)),
        paste0("export POE_SITE_FILE=", shQuote(cfg$execution$site_file)),
        paste0("result=$(", reg_command(c(entry, "worker-fish-begin", args)), ")"),
        'printf "%s\\n" "$result"',
        'if [[ "${result##*$\'\\n\'}" == RUN=0 ]]; then exit 0; fi',
        '[[ "${result##*$\'\\n\'}" == RUN=1 ]] || exit 1',
        reg_fish_commands(cfg, z, phase, res)[-1L],
        reg_command(c(entry, "worker-fish-finish", "--phase", phase)))
      writeLines(lines, scripts[i])
    }
  } else {
    scripts <- file.path(job, "run.sh")
    writeLines(c("#!/bin/bash", "set -eo pipefail", "umask 077",
      reg_command(c(entry, paste0("worker-", action)))), scripts)
  }
  host <- file.path(job, "host.sh")
  env <- c(paste0("export SCRIPTS_DIR=", shQuote(root)),
    paste0("export POE_SITE_FILE=", shQuote(cfg$execution$site_file)))
  lines <- c("#!/bin/bash", "set -eo pipefail", "umask 077", env,
    paste("mkdir", shQuote(file.path(job, "submitted.guard"))),
    "declare -A seen", "previous=''")
  for (i in seq_along(scripts)) {
    name <- paste(cfg$cohort_id, action, i, basename(job), sep = "_")
    args <- c("-q", r$queue, "-n", r$cpus, "-M", r$memory_value,
      "-R", paste0("rusage[mem=", r$reservation_value, "] span[hosts=1]"),
      "-W", r$wall_time, "-J", name, "-o", file.path(job, paste0(i, ".%J.out")),
      "-e", file.path(job, paste0(i, ".%J.err")))
    dep <- if (i > 1L) ' -w "done($previous)"' else
      if (after != "none") paste0(" -w 'done(", after, ")'") else ""
    cmd <- paste0(reg_command(c(cfg$tools$bsub, args)), dep, " < ", shQuote(scripts[i]))
    lines <- c(lines, paste0("output=$(", cmd, ")"),
      'printf "%s\\n" "$output"',
      '[[ "$output" =~ Job\\ \\<([0-9]+)\\> ]] || { printf "Missing job ID; inspect LSF before retry\\n" >&2; exit 1; }',
      'previous="${BASH_REMATCH[1]}"',
      '[[ -z "${seen[$previous]:-}" ]] || exit 1', 'seen[$previous]=yes',
      paste0('printf "%s\\n" "$previous" > ', shQuote(file.path(job, paste0(i, ".jobid")))))
  }
  writeLines(lines, host)
  cat(if (submit) "Prepared host-only plan:" else "DRY RUN ONLY; prepared plan:", host, "\n")
  if (submit) cat("HOST_PLAN=", host, "\n", sep = "")
  invisible(host)
}
reg_main <- function(args, root) {
  if (length(args) < 2L) stop("Usage: poe_regulatory.R CONFIG ACTION [--submit] [--after JOB] [--resume DIR]")
  for (module in c("config", "common", "workflow", "parental", "phase_origin", "gusto",
      "models", "exports", "acquisition", "public", "reference", "pipeline", "stages",
      "additional", "scheduler", "seadragon", "lsf", "annotation",
      "regulatory_common", "regulatory_resources", "regulatory_fish",
      "regulatory_te", "regulatory_context"))
    source(file.path(root, "R", paste0(module, ".R")), local = .GlobalEnv)
  cfg <- read_config_v2(args[1])
  if (dir.exists(cfg$tools$r_library)) .libPaths(cfg$tools$r_library, include.site = FALSE)
  z <- reg_settings(cfg); action <- args[2]; flags <- args[-c(1, 2)]
  if (anyDuplicated(flags[startsWith(flags, "--")])) stop("Duplicate flag")
  values <- list()
  while (length(flags)) {
    key <- flags[1]
    if (key == "--submit") {values[[key]] <- TRUE; flags <- flags[-1L]; next}
    if (!key %in% c("--after", "--resume", "--phase", "--manifest") ||
        length(flags) < 2L || startsWith(flags[2], "--")) stop("Invalid/missing flag value: ", key)
    values[[key]] <- flags[2]; flags <- flags[-c(1, 2)]
  }
  if (action %in% reg_actions) {
    if (any(c("--phase", "--manifest") %in% names(values))) stop("Worker-only flag")
    if (action == "regulatory-status") {
      if (length(values)) stop("Status accepts no flags")
      print(reg_status(z)); return(invisible(TRUE))
    }
    return(reg_plan_jobs(cfg, z, root, action, isTRUE(values[["--submit"]]),
      reg_default(values[["--after"]], "none"), values[["--resume"]]))
  }
  if (any(c("--submit", "--after") %in% names(values))) stop("Submission flag on worker")
  if (action == "worker-resources") {
    if (is.null(values[["--manifest"]])) stop("--manifest required")
    reg_finalize_resources(cfg, z, root, values[["--manifest"]])
  } else if (action %in% c("worker-fish-begin", "worker-fish-finish")) {
    phase <- values[["--phase"]]
    if (is.null(phase) || !phase %in% c("discovery", "masking")) stop("--phase required")
    if (action == "worker-fish-begin")
      cat("RUN=", as.integer(reg_fish_begin(cfg, z, root, phase, values[["--resume"]])),
        "\n", sep = "")
    else reg_fish_finish(cfg, z, root, phase)
  } else if (action == "worker-te-tss") reg_te_worker(cfg, z, root)
  else if (action == "worker-human-context") reg_context_worker(cfg, z, root)
  else stop("Unknown regulatory action")
}
