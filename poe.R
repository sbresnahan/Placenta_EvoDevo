#!/usr/bin/env Rscript
# Repository location is bootstrap metadata, never a data-path fallback.
f<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
f<-gsub("~+~"," ",f,fixed=TRUE)
root<-normalizePath(dirname(f),mustWork=TRUE)
a<-commandArgs(trailingOnly=TRUE)
if(length(a)<2)stop("Usage: Rscript poe.R CONFIG STAGE [--submit] [--execute] [--after JOB] [--accept SHA256]")
source(file.path(root,"R/config.R"))
c<-read_config_v2(a[1]);stage<-a[2];flags<-a[-c(1,2)]
if(dir.exists(c$tools$r_library)).libPaths(c$tools$r_library,include.site=FALSE)
for(module in c("common","workflow","parental","models","exports","acquisition",
                 "public","reference","pipeline","stages","scheduler","seadragon","lsf"))
  source(file.path(root,"R",paste0(module,".R")))
allowed<-c("--submit","--execute","--after","--accept")
value<-function(flag) {
  k<-match(flag,flags);if(is.na(k))return(NULL)
  if(k==length(flags)||startsWith(flags[k+1],"--"))stop("Missing value for ",flag)
  flags[k+1]
}
extras<-flags
for(flag in c("--after","--accept"))if(flag %in% extras) {
  k<-match(flag,extras);if(k==length(extras))stop("Missing flag value")
  extras<-extras[-c(k,k+1)]
}
if(any(!extras %in% c("--submit","--execute")))stop("Unknown CLI flag")
if(anyDuplicated(flags[startsWith(flags,"--")]))stop("Duplicate CLI flag")
data.table::setDTthreads(as.integer(c$ornaments$threads))
tryCatch({
  if(stage=="preflight")preflight_v2(c,root,"--execute" %in% flags)
  else if(stage=="retrieve") {
    preflight_v2(c,root)
    if(is_seadragon(c)) {
      if("--execute" %in% flags)stop("Run 01_retrieval/01_retrieve_public.lsf directly on the approved host; R cannot launch native host tools")
      submit_graph(c,root,"retrieve","--submit" %in% flags,
        if(is.null(value("--after")))"none" else value("--after"))
    } else if(!"--execute" %in% flags)cat("DRY RUN ONLY: retrieve requires --execute on the approved download host\n")
    else {snapshot(c);retrieve_inputs(c)}
  } else if(stage %in% c("reference","pilot","scale","models")) {
    preflight_v2(c,root)
    if(stage=="scale"&&"--submit" %in% flags)authorize_scale(c,root,value("--accept"))
    submit_stage(c,root,stage,"--submit" %in% flags,
      if(is.null(value("--after")))"none" else value("--after"))
  } else if(startsWith(stage,"worker-")) {
    if(is_seadragon(c))stop("Use checked-in LSF workers; legacy all-library native execution is disabled")
    # Workers are reached through explicit stage commands; scale always checks approval.
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
  if(stage!="export")write_tsv(data.table(stage=stage,
    state=if(startsWith(stage,"worker-"))"failed" else "blocked",reason=conditionMessage(e)),c$paths$status)
  stop(conditionMessage(e),call.=FALSE)
})
