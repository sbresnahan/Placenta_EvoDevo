#!/usr/bin/env Rscript
# R counterpart to xQTL's driver; prepares checked-in-worker graph for HOST bsub.
f<-gsub("~+~"," ",sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]),fixed=TRUE)
root<-normalizePath(dirname(f),mustWork=TRUE)
a<-commandArgs(trailingOnly=TRUE)
if(length(a)<2L)stop("Usage: run_pipeline.R CONFIG retrieve|reference|pilot|scale|models [--submit] [--accept SHA] [--after JOB]")
source(file.path(root,"R/config.R"));c<-read_config_v2(a[1])
if(dir.exists(c$tools$r_library)).libPaths(c$tools$r_library,include.site=FALSE)
for(k in c("common","workflow","parental","models","exports","acquisition","public","reference",
           "pipeline","stages","scheduler","seadragon","lsf"))source(file.path(root,"R",paste0(k,".R")))
flags<-a[-c(1,2)]
if(anyDuplicated(flags[startsWith(flags,"--")]))stop("Duplicate flag")
take<-function(k) {
  i<-match(k,flags);if(is.na(i))return(NULL)
  if(i==length(flags)||startsWith(flags[i+1],"--"))stop("Missing value: ",k)
  flags[i+1]
}
after<-take("--after");accept<-take("--accept")
extra<-flags
for(k in c("--after","--accept"))if(k %in% extra) {
  i<-match(k,extra);extra<-extra[-c(i,i+1)]
}
if(any(!extra %in% "--submit"))stop("Unknown flag")
if(!is_seadragon(c))stop("run_pipeline.R requires the Seadragon backend")
preflight_v2(c,root)
if(a[2]=="scale"&&"--submit" %in% flags)authorize_scale(c,root,accept)
submit_graph(c,root,a[2],"--submit" %in% flags,if(is.null(after))"none" else after)
