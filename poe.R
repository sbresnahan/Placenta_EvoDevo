#!/usr/bin/env Rscript
# Locate the scripts and initialize the R library.
f<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
f<-gsub("~+~"," ",f,fixed=TRUE)
root<-normalizePath(dirname(f),mustWork=TRUE)
source(file.path(root,"R/library_paths.R"))
a<-commandArgs(trailingOnly=TRUE)
if(length(a)<2)stop("Usage: Rscript poe.R CONFIG STAGE [--submit] [--execute] [--after JOB] [--resume] [--sample ID]")
if(a[2] %in% c("regulatory-resources","fish-repeats","te-tss","te-genes","human-context","regulatory-status")) {
  source(file.path(root,"R/regulatory_cli.R"))
  reg_main(a,root)
  quit(status=0L)
}
source(file.path(root,"R/config.R"))
stage<-a[2];c<-read_config_v2(a[1]);flags<-a[-c(1,2)]
if(dir.exists(c$tools$r_library)).libPaths(c$tools$r_library,include.site=FALSE)
for(module in c("common","workflow","parental","phase_origin","gusto","models","exports","acquisition",
                 "public","reference","pipeline","stages","additional","scheduler","seadragon","lsf"))
  source(file.path(root,"R",paste0(module,".R")))
allowed<-c("--submit","--execute","--after","--resume","--sample")
value<-function(flag) {
  k<-match(flag,flags);if(is.na(k))return(NULL)
  if(k==length(flags)||startsWith(flags[k+1],"--"))stop("Missing value for ",flag)
  flags[k+1]
}
extras<-flags
for(flag in c("--after","--sample"))if(flag %in% extras) {
  k<-match(flag,extras);if(k==length(extras))stop("Missing flag value")
  extras<-extras[-c(k,k+1)]
}
if(any(!extras %in% c("--submit","--execute","--resume")))stop("Unknown CLI flag")
if(anyDuplicated(flags[startsWith(flags,"--")]))stop("Duplicate CLI flag")
if("--resume" %in% flags && !stage %in% c("reference","pilot"))
  stop("--resume requires reference or pilot stage")
if(!is.null(value("--sample")) && !"--resume" %in% flags)
  stop("--sample requires reference --resume")
data.table::setDTthreads(as.integer(c$ornaments$threads))
tryCatch({
  if(stage=="retrieve") {
    validate_inputs_v2(c,root)
    if(is_seadragon(c)) {
      if("--execute" %in% flags)stop("Use bin/submit_pipeline.sh CONFIG retrieve --submit")
      submit_graph(c,root,"retrieve","--submit" %in% flags,
        if(is.null(value("--after")))"none" else value("--after"))
    } else if(!"--execute" %in% flags)cat("DRY RUN ONLY: retrieve requires --execute\n")
    else {snapshot(c);retrieve_inputs(c)}
  } else if(stage %in% c("reference","pilot","scale","models","additional")) {
    validate_inputs_v2(c,root)
    if(stage=="scale"&&"--submit" %in% flags)request_scale(c,root)
    after<-if(is.null(value("--after")))"none" else value("--after")
    if("--resume" %in% flags)
      submit_graph(c,root,stage,"--submit" %in% flags,after,resume=TRUE,sample=value("--sample"))
    else submit_stage(c,root,stage,"--submit" %in% flags,after)
  } else if(startsWith(stage,"worker-")) {
    if(is_seadragon(c))stop("Use bin/submit_pipeline.sh CONFIG STAGE --submit")
    # Each worker runs an explicitly requested stage.
    suppressPackageStartupMessages(library(GenomicRanges))
    source(file.path(root,"R/annotation.R"))
    action<-sub("^worker-","",stage)
    if(!action %in% c("reference","pilot","scale","models","additional"))stop("Unknown worker stage")
    get(paste0("worker_",action))(c,root)
  } else if(stage=="review")review_v2(c)
  else if(stage=="export")export_v2(c)
  else if(stage=="status") {
    status_v2(c,root)
  } else if(stage=="additional-status") {
    paths<-additional_paths(c)
    ok<-stage_valid(additional_receipt(c),paths)
    print(data.table(stage="additional_v1",state=if(ok)"complete" else
      if(file.exists(additional_receipt(c)))"incomplete" else "not_requested",
      output_dir=additional_dir(c)))
  } else stop("Unknown stage")
},error=function(e) {
  if(stage!="export")write_tsv(data.table(stage=stage,
    state=if(startsWith(stage,"worker-"))"failed" else "blocked",reason=conditionMessage(e)),c$paths$status)
  stop(conditionMessage(e),call.=FALSE)
})
