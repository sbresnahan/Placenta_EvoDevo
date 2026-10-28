# Base-R-only startup: run before parsing YAML or loading any other dependency.
# POE_R_LIBRARY is passed explicitly by Rscript_sif; direct Rscript calls use
# the configured deployment default unless the caller selects a library.
local({
  library_path <- Sys.getenv("POE_R_LIBRARY", unset = NA_character_)
  if (is.na(library_path)) {
    # Production entry points already establish root from their --file metadata.
    if (!exists("root", inherits = TRUE))
      stop("R bootstrap requires the repository root or POE_R_LIBRARY.", call. = FALSE)
    default_file <- file.path(root, "config/r_library.path")
    if (!file.exists(default_file))
      stop("Missing direct-Rscript library setting: ", default_file, call. = FALSE)
    library_path <- readLines(default_file, warn = FALSE)
    if (length(library_path) != 1L)
      stop("config/r_library.path must contain one approved library path.", call. = FALSE)
  }
  if (!nzchar(library_path) || !dir.exists(library_path) ||
      file.access(library_path, mode = 5L) != 0L)
    stop("R library is missing or not readable/searchable: ", library_path,
         ". Set POE_R_LIBRARY to the approved R library.", call. = FALSE)
  library_path <- normalizePath(library_path, mustWork = TRUE)
  .libPaths(c(library_path, .libPaths()))
  if (!identical(.libPaths()[1L], library_path))
    stop("Failed to prepend the approved R library: ", library_path, call. = FALSE)
  invisible(library_path)
})
