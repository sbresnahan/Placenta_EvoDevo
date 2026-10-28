suppressPackageStartupMessages(library(data.table))
Sys.umask("0077")
options(stringsAsFactors = FALSE)
data.table::setDTthreads(1L)

require_columns <- function(x, required, label = "input") {
  if (anyDuplicated(names(x))) stop(label, ": duplicate column names")
  miss <- setdiff(required, names(x))
  if (length(miss)) stop(label, ": missing columns: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}
assert_unique <- function(x, cols, label = "input") {
  if (anyNA(x[, ..cols]) || anyDuplicated(x[, ..cols]))
    stop(label, ": missing or duplicate key: ", paste(cols, collapse = ", "))
  invisible(TRUE)
}
sha256 <- function(path) digest::digest(file = path, algo = "sha256")
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", na = "NA")
}
atomic_tsv <- function(x, path) {
  tmp <- tempfile(pattern = ".poe-", tmpdir = dirname(path))
  data.table::fwrite(x, tmp, sep = "\t", na = "NA")
  if (!file.rename(tmp, path)) stop("Atomic rename failed; use a POSIX work directory")
}
script_root <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  f <- gsub("~+~"," ",f,fixed=TRUE)
  normalizePath(file.path(dirname(f), ".."), mustWork = TRUE)
}
cli_args <- function(required) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) != length(required))
    stop("Arguments: ", paste(required, collapse = " "))
  setNames(as.list(a), required)
}
provenance <- function(paths) {
  data.table(file = basename(paths), bytes = file.info(paths)$size,
             sha256 = vapply(paths, sha256, character(1)))
}
