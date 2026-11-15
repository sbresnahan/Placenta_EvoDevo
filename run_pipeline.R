#!/usr/bin/env Rscript
# R counterpart to xQTL's driver; prepares checked-in-worker graph for HOST bsub.
f<-gsub("~+~"," ",sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]),fixed=TRUE)
root<-normalizePath(dirname(f),mustWork=TRUE)
source(file.path(root,"R/library_paths.R"))
a<-commandArgs(trailingOnly=TRUE)
if(length(a)<2L)stop("Usage: run_pipeline.R CONFIG retrieve|reference|pilot|scale|models|additional [--submit] [--resources-only] [--after JOB] [--resume] [--sample ID]")
if(a[2] %in% c("regulatory-resources","fish-repeats","te-tss","te-genes","human-context","regulatory-status")) {
  source(file.path(root,"R/regulatory_cli.R"))
  reg_main(a,root)
  quit(status=0L)
}
source(file.path(root,"R/config.R"));c<-read_config_v2(a[1])
if(dir.exists(c$tools$r_library)).libPaths(c$tools$r_library,include.site=FALSE)
for(k in c("common","workflow","parental","phase_origin","gusto","models","exports","acquisition","public","reference",
           "pipeline","stages","additional","scheduler","seadragon","lsf"))source(file.path(root,"R",paste0(k,".R")))
flags<-a[-c(1,2)]
if(anyDuplicated(flags[startsWith(flags,"--")]))stop("Duplicate flag")
take<-function(k) {
  i<-match(k,flags);if(is.na(i))return(NULL)
  if(i==length(flags)||startsWith(flags[i+1],"--"))stop("Missing value: ",k)
  flags[i+1]
}
after<-take("--after")
sample<-take("--sample")
extra<-flags
for(k in c("--after","--sample"))if(k %in% extra) {
  i<-match(k,extra);extra<-extra[-c(i,i+1)]
}
if(any(!extra %in% c("--submit","--resume","--resources-only")))stop("Unknown flag")
resources_only<-"--resources-only" %in% flags
if(resources_only && a[2]!="retrieve")stop("--resources-only requires retrieve stage")
resume<-"--resume" %in% flags
if(resume && !a[2] %in% c("reference","pilot"))
  stop("--resume requires reference or pilot stage")
if(!is.null(sample) && !resume)stop("--sample requires reference --resume")
if(!is_seadragon(c))stop("run_pipeline.R requires the Seadragon backend")
validate_inputs_v2(c,root,resources_only=resources_only)
if(a[2]=="scale"&&"--submit" %in% flags)request_scale(c,root)
submit_graph(c,root,a[2],"--submit" %in% flags,if(is.null(after))"none" else after,
             resume=resume,sample=sample,resources_only=resources_only)
