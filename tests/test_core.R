#!/usr/bin/env Rscript
# Logical/coordinate truth tables only; no simulated biological observations.
f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
f <- gsub("~+~", " ", f, fixed=TRUE)
root <- normalizePath(file.path(dirname(f), ".."))
source(file.path(root, "R/common.R"))
source(file.path(root, "R/parental.R"))
source(file.path(root, "R/models.R"))
n <- 0L
check <- function(x) {stopifnot(isTRUE(x)); n <<- n + 1L}
fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
g <- c("0/0","0/1","1/1")
for (mother in g) for (father in g) for (child in g) {
  possibilities <- expand.grid(m = unique(parse_gt(mother)), p = unique(parse_gt(father)))
  ok <- vapply(seq_len(nrow(possibilities)), function(i)
    paste(sort(unlist(possibilities[i,])), collapse = "") ==
      paste(sort(parse_gt(child)), collapse = ""), logical(1))
  actual <- transmitted_maternal_allele(child, mother, father)
  expected <- if (!any(ok)) "mendelian_conflict" else if (child != "0/1")
    "child_not_heterozygous" else if (length(unique(possibilities$m[ok])) > 1)
      "ambiguous_transmission" else "typed_trio"
  check(identical(actual$status, expected))
  if (expected == "typed_trio")
    check(actual$maternal == unique(possibilities$m[ok]))
}
check(transmitted_maternal_allele("0|1","0/0","1/1")$maternal == 0L)
check(transmitted_maternal_allele("./.","0/0","1/1")$status == "missing_or_nonbiallelic")
# Coordinate truth table including exon boundaries and negative strand.
ex <- data.table(transcript_id = rep(c("plus","minus"), each = 2),
  chrom = "chr1", start = c(10L,30L,30L,10L), end = c(19L,39L,39L,19L),
  strand = rep(c("+","-"), each = 2), tx_exon_order = rep(1:2,2))
counts <- data.table(transcript_id = rep(c("plus","minus"), each = 3),
  tx_pos = rep(c(1L,10L,11L),2), tx_ref = "A", tx_alt = "G",
  ref_count = 0, alt_count = 2.5, other_count = 0)
m <- map_tx_loci(counts, ex)
setorder(m, transcript_id, tx_pos)
check(identical(as.integer(m$genomic_pos), c(39L,30L,19L,10L,19L,30L)))
check(all(m[strand == "-"]$ref == "T") && all(m[strand == "-"]$alt == "C"))
trio <- unique(m[, .(chrom, genomic_pos, ref, alt)])
trio[, `:=`(child_gt = "0/1", mother_gt = "0/0", father_gt = "1/1")]
oriented <- orient_counts(m, trio)
check(all(oriented$maternal_count == 0) && all(oriented$paternal_count == 2.5))
trio[1, child_gt := "./."]
missing <- orient_counts(m, trio)
check(anyNA(missing$maternal_count))
check(fails(orient_counts(m, rbind(trio,trio[1]))))
tx <- data.table(transcript_id = c("plus","minus"), event_id = c("P","M"), gene_id = c("G","G"))
ev <- event_locus_counts(oriented, tx, "coordinate_test", "truth_table")
check(nrow(ev) == 6L) # NOT collapsed across loci.
check(sum(ev$maternal_count) == 0)
check(fails(read_ornament_counts(tempfile(), "left_right")))
check(fails(assert_unique(data.table(id = c("a","a")), "id")))
check(equivalence_call(0, .01, 20)$equivalent)
check(!equivalence_call(0, 1, 20)$equivalent)
check(is.null(hc3_intercept(1:10, cbind(1,1))))
cat("PASS:", n, "logical/coordinate assertions\n")
