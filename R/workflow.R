read_config <- function(path) {
  read_config_v2(path)
}
run_checked <- function(command, args, stdout = "", stderr = "") {
  status <- system2(command, vapply(args, shQuote, character(1)), stdout = stdout, stderr = stderr)
  if (!identical(as.integer(status), 0L)) stop("Command failed: ", basename(command))
  invisible(TRUE)
}
stage_signature <- function(inputs, settings, source_files = character()) {
  paths <- c(inputs, source_files)
  if (any(!file.exists(paths)) || any(file.info(paths)$size <= 0))
    stop("Missing/empty stage inputs")
  digest::digest(list(hashes = vapply(paths, sha256, character(1)), settings = settings),
                 algo = "sha256")
}
stage_valid <- function(receipt, signature, outputs) {
  if (!file.exists(receipt) || any(!file.exists(outputs))) return(FALSE)
  r <- tryCatch(readRDS(receipt), error = function(e) NULL)
  !is.null(r) && identical(r$signature, signature) &&
    identical(r$output_hashes, setNames(vapply(unname(outputs), sha256, character(1)),unname(outputs))) &&
    all(file.info(outputs)$size > 0)
}
stage_finish <- function(receipt, signature, outputs) {
  if (any(!file.exists(outputs)) || any(file.info(outputs)$size <= 0))
    stop("Cannot commit missing/empty stage outputs")
  dir.create(dirname(receipt),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile(".receipt-",tmpdir=dirname(receipt))
  saveRDS(list(signature = signature,
    state = "complete", outputs=unname(outputs),
    output_hashes = setNames(vapply(unname(outputs), sha256, character(1)),unname(outputs))), tmp)
  if(!file.rename(tmp,receipt))stop("Atomic receipt commit failed")
}
check_linkage <- function(samples, trios) {
  require_columns(samples, c("sample_id","child_id","layout","strand_mode","fastq_1","fastq_2",
                            "identity_verified"))
  require_columns(trios, c("child_id","mother_id","father_id","family_id"))
  assert_unique(samples, "sample_id", "RNA manifest")
  assert_unique(samples, "child_id", "One primary RNA library per child")
  assert_unique(trios, "child_id", "Trio manifest")
  if (anyNA(trios) || any(trios$child_id == trios$mother_id |
      trios$child_id == trios$father_id | trios$mother_id == trios$father_id))
    stop("Missing/conflicting pedigree identities")
  if (any(!samples$child_id %in% trios$child_id)) stop("RNA/trio linkage incomplete")
  if (anyNA(samples$identity_verified) || !all(samples$identity_verified))
    stop("Independent RNA/DNA identity QC must be signed off locally")
  if (any(!samples$layout %in% c("PE","SE"))) stop("Unknown layout")
  if (anyNA(samples$strand_mode) || any(!samples$strand_mode %in% c("unstranded","fr","rf")))
    stop("Explicit independently verified Ornaments library strand_mode required")
  se<-samples[layout=="SE"]
  if(nrow(se) && (!all(c("fragment_mean","fragment_sd") %in% names(se)) ||
     any(!is.finite(c(se$fragment_mean,se$fragment_sd))) ||
     any(se$fragment_mean<=0 | se$fragment_sd<=0)))
    stop("SE library fragment mean/SD missing")
  if (any(!file.exists(samples$fastq_1)) ||
      any(!file.exists(samples[layout == "PE"]$fastq_2))) stop("FASTQs missing")
  merge(samples, trios, by = "child_id", all.x = TRUE, sort = FALSE)
}
