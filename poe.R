#!/usr/bin/env Rscript
# Locate the scripts and initialize the R library.
f<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
f<-gsub("~+~"," ",f,fixed=TRUE)
root<-normalizePath(dirname(f),mustWork=TRUE)
source(file.path(root,"R/library_paths.R"))
a<-commandArgs(trailingOnly=TRUE)
if(length(a)<2)stop("Usage: Rscript poe.R CONFIG STAGE [--submit] [--execute] [--after JOB]")
source(file.path(root,"R/config.R"))
stage<-a[2];c<-read_config_v2(a[1],setup=stage=="preflight");flags<-a[-c(1,2)]
if(dir.exists(c$tools$r_library)).libPaths(c$tools$r_library,include.site=FALSE)
for(module in c("common","workflow","parental","models","exports","acquisition",
                 "public","reference","pipeline","stages","scheduler","seadragon","lsf"))
  source(file.path(root,"R",paste0(module,".R")))
allowed<-c("--submit","--execute","--after")
value<-function(flag) {
  k<-match(flag,flags);if(is.na(k))return(NULL)
  if(k==length(flags)||startsWith(flags[k+1],"--"))stop("Missing value for ",flag)
  flags[k+1]
}
extras<-flags
for(flag in "--after")if(flag %in% extras) {
  k<-match(flag,extras);if(k==length(extras))stop("Missing flag value")
  extras<-extras[-c(k,k+1)]
}
if(any(!extras %in% c("--submit","--execute")))stop("Unknown CLI flag")
if(anyDuplicated(flags[startsWith(flags,"--")]))stop("Duplicate CLI flag")
data.table::setDTthreads(if(stage=="preflight")1L else as.integer(c$ornaments$threads))
tryCatch({
  if(stage=="preflight")preflight_v2(c,root,"--execute" %in% flags)
  else if(stage=="retrieve") {
    validate_inputs_v2(c,root)
    if(is_seadragon(c)) {
      if("--execute" %in% flags)stop("Use bin/submit_pipeline.sh CONFIG retrieve --submit")
      submit_graph(c,root,"retrieve","--submit" %in% flags,
        if(is.null(value("--after")))"none" else value("--after"))
    } else if(!"--execute" %in% flags)cat("DRY RUN ONLY: retrieve requires --execute\n")
    else {snapshot(c);retrieve_inputs(c)}
  } else if(stage %in% c("reference","pilot","scale","models")) {
    validate_inputs_v2(c,root)
    if(stage=="scale"&&"--submit" %in% flags)request_scale(c,root)
    submit_stage(c,root,stage,"--submit" %in% flags,
      if(is.null(value("--after")))"none" else value("--after"))
  } else if(startsWith(stage,"worker-")) {
    if(is_seadragon(c))stop("Use bin/submit_pipeline.sh CONFIG STAGE --submit")
    # Each worker runs an explicitly requested stage.
    suppressPackageStartupMessages(library(GenomicRanges))
    source(file.path(root,"R/annotation.R"))
    action<-sub("^worker-","",stage)
    if(!action %in% c("reference","pilot","scale","models"))stop("Unknown worker stage")
    get(paste0("worker_",action))(c,root)
  } else if(stage=="review")review_v2(c)
  else if(stage=="export")export_v2(c)
  else if(stage=="status") {
    status_v2(c,root)
  } else stop("Unknown stage")
},error=function(e) {
  if(!stage %in% c("preflight","export"))write_tsv(data.table(stage=stage,
    state=if(startsWith(stage,"worker-"))"failed" else "blocked",reason=conditionMessage(e)),c$paths$status)
  stop(conditionMessage(e),call.=FALSE)
})
