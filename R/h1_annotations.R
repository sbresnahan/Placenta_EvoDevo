# Optional descriptive annotations; no count or model transformation.
h1_normalize_annotations <- function(x, label = "h1_bayes") {
  if (is.null(x)) return(list())
  if (!is.list(x) || (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x)))))
    stop(label, " must be a YAML mapping")
  shipped <- list(count_provenance = "REVIEW_REQUIRED_COUNT_PROVENANCE",
    producing_pipeline = c("REVIEW_REQUIRED_COUNT_PRODUCER", "DECLARE_REVIEWED_COUNT_PRODUCER"))
  for (key in names(shipped)) {
    value <- x[[key]]
    absent <- is.null(value) || !length(value) ||
      (length(value) == 1L && is.atomic(value) && is.na(value))
    if (!absent && is.character(value) && length(value) == 1L)
      absent <- !nzchar(trimws(value)) || value %in% shipped[[key]]
    if (absent) {
      x[[key]] <- NULL
      next
    }
    if (!is.character(value) || length(value) != 1L || is.na(value))
      stop(label, ".", key, " must be one optional text annotation")
    if (key == "count_provenance" && !value %in% c("fragment_unique_assignment",
        "upstream_duplicate_corrected", "assignment_covariance_unknown"))
      stop(label, ".count_provenance must be one of: fragment_unique_assignment, ",
        "upstream_duplicate_corrected, assignment_covariance_unknown; or absent")
    if (key == "producing_pipeline" &&
        grepl("REVIEW_REQUIRED|DECLARE|TODO|UNSPECIFIED", value, ignore.case = TRUE))
      stop(label, ".producing_pipeline must identify a recorded producer or be absent")
  }
  x
}
h1_validate_controls <- function(h, label = "h1_bayes") {
  integers <- c("chains", "warmup", "iterations", "parallel_chains", "max_treedepth", "seed")
  for (key in intersect(names(h), integers)) {
    value <- h[[key]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value))
      stop("Invalid ", label, ".", key, ": expected a positive integer")
  }
  chains <- if (is.null(h$chains)) 4L else h$chains
  parallel <- if (is.null(h$parallel_chains)) 4L else h$parallel_chains
  if (chains < 4L || parallel > chains)
    stop(label, " requires chains >= 4 and parallel_chains <= chains")
  if ("adapt_delta" %in% names(h) && (!is.numeric(h$adapt_delta) ||
      length(h$adapt_delta) != 1L || !is.finite(h$adapt_delta) ||
      h$adapt_delta < .95 || h$adapt_delta >= 1))
    stop(label, ".adapt_delta must be >= 0.95 and < 1")
  invisible(h)
}
h1_merge_overrides <- function(inherited, operator, label = "h1_bayes") {
  out <- h1_normalize_annotations(inherited, label)
  overlay <- h1_normalize_annotations(operator, "operator.h1_bayes")
  for (key in names(overlay)) out[[key]] <- overlay[[key]]
  h1_validate_controls(out, label)
  out
}
