#!/usr/bin/env Rscript
# Input test for one cohort: every dataset the pipeline reads after retrieval
# is present (and, where the manifests record it, has the expected byte size).
# Run through tests/test_inputs.sh CONFIG, which starts the R container.
f<-gsub("~+~"," ",sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]),fixed=TRUE)
root<-normalizePath(file.path(dirname(f),".."),mustWork=TRUE)
source(file.path(root,"R/library_paths.R"))
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=1L)stop("Usage: tests/test_inputs.sh CONFIG")
source(file.path(root,"R/config.R"));c<-read_config_v2(a[1])
if(dir.exists(c$tools$r_library)).libPaths(c$tools$r_library,include.site=FALSE)
for(k in c("common","workflow","parental","phase_origin","gusto","models","exports","acquisition","public",
           "reference","pipeline","stages","additional","scheduler","seadragon","lsf"))
  source(file.path(root,"R",paste0(k,".R")))

problems<-character();checked<-0L
problem<-function(...)problems<<-c(problems,paste0(...))
attempt<-function(label,expr)tryCatch(expr,error=function(e){problem(label,": ",conditionMessage(e));NULL})
size_of<-function(p)if(file.exists(p))file.info(p)$size else NA_real_
need<-function(path,label,bytes=NA_real_) {
  checked<<-checked+1L
  if(is.null(path)||is.na(path)||!nzchar(path))return(invisible(problem(label,": no path configured")))
  n<-size_of(path)
  why<-if(is.na(n))"missing" else if(n<=0)"empty" else
    if(is.finite(bytes)&&n!=bytes)sprintf("%s bytes, expected %s",
      format(n,scientific=FALSE),format(bytes,scientific=FALSE)) else ""
  if(nzchar(why))problem(label,": ",path," (",why,")")
  invisible(!nzchar(why))
}

m<-attempt("manifests",load_manifests(c))

if(c$profile=="public_cross") {
  # Raw reads: ENA byte sizes, or the SRA receipt written by bin/fetch_sra.sh.
  if(!is.null(m))for(i in seq_len(nrow(m$sources)))
    for(mate in seq_len(if(m$sources$layout[i]=="PE")2L else 1L)) {
      s<-m$sources[i];p<-cfg_path(c,"fastq",run=s$run_id,mate=mate);checked<-checked+1L
      if(!isTRUE(verify_public_source(c,s,mate)))
        problem("FASTQ ",s$run_id,"_",mate,": ",p,
          if(file.exists(p))" (size does not match manifest or SRA receipt)" else " (missing)")
    }
  # Reference resources listed as required in the resources manifest.
  items<-attempt("resources manifest",host_retrieval_items(c,resources_only=TRUE))
  if(!is.null(items))for(i in seq_len(nrow(items))) {
    x<-items[i]
    bytes<-if(x$method=="cyverse") {
      r<-paste0(x$destination,".cyverse.tsv")
      d<-if(file.exists(r))tryCatch(fread(r,colClasses="character"),error=function(e)NULL)
      if(is.null(d)||!"bytes" %in% names(d)) {problem("resource ",x$id,": no CyVerse receipt ",r);NA_real_}
      else as.numeric(d$bytes[1])
    } else x$bytes
    need(x$destination,paste("resource",x$id),bytes)
    if(!is.na(x$decompress_to))need(x$decompress_to,paste("resource",x$id,"(decompressed)"))
  }
  need(c$paths$resource_receipt,"retrieval receipt")
  # Parental marker sources (AGPv3 sources are lifted over in the markers stage).
  ms<-attempt("marker sources",resolve_manifest_paths(fread(c$paths$marker_sources),c("path","processed_path"),c))
  if(!is.null(ms)) {
    if(!nrow(ms))problem("marker sources: ",c$paths$marker_sources," has no rows")
    for(i in seq_len(nrow(ms)))need(ms$path[i],paste("marker source",ms$source_id[i]))
    if(any(ms$build=="AGPv3"))need(c$markers$liftover_chain,"AGPv3-to-v4 liftover chain")
  }
} else {
  # Human trio: private local FASTQs and genotype/origin inputs (never retrieved).
  if(!is.null(m)) {
    attempt("trio linkage",human_linkage(c,m))
    for(i in seq_len(nrow(m$sources)))for(mate in seq_len(if(m$sources$layout[i]=="PE")2L else 1L)) {
      s<-m$sources[i]
      p<-attempt(paste0("FASTQ ",s$run_id,"_",mate),source_fastq(c,s,mate))
      if(!is.null(p))need(p,paste0("FASTQ ",s$run_id,"_",mate))
    }
  }
  if(gusto_enabled(c)) {
    for(k in c("vcf","origin_matrix","fastq_id_map"))
      need(c$gusto_precomputed[[k]],paste0("gusto_precomputed.",k))
  } else for(k in c("typed_vcf","personalization_vcf"))need(c$paths[[k]],paste0("paths.",k))
  need(c$paths$trios,"trio manifest")
  need(c$paths$covariates,"covariates manifest")
}

# Reference genome/annotation read by the reference stage.
if(file.exists(c$paths$genome)||is.null(c$reference$genome_source))need(c$paths$genome,"genome FASTA") else
  need(c$reference$genome_source,"reference.genome_source")
need(c$paths$annotation,"annotation")
if(!is.null(c$reference$transcript_source))need(c$reference$transcript_source,"reference.transcript_source")
if(!is.null(c$additional_analysis$te_rmsk))need(c$additional_analysis$te_rmsk,"additional_analysis.te_rmsk")

if(length(problems)) {
  cat(sprintf("%s: %d of %d input checks FAILED\n",c$cohort_id,length(problems),checked))
  cat(paste0("  - ",problems),sep="\n")
  quit(status=1L)
}
cat(sprintf("%s: all %d input checks passed; inputs are ready.\n",c$cohort_id,checked))
