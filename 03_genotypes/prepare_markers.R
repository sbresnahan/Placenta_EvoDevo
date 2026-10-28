#!/usr/bin/env Rscript
# Stage phase entry point. Native operations are emitted for the host, never run in R.
f<-gsub("~+~"," ",sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]),fixed=TRUE)
root<-normalizePath(file.path(dirname(f),".."),mustWork=TRUE)
source(file.path(root,"R/library_paths.R"))
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=4L)stop("Usage: stage.R CONFIG PHASE LIBRARY JOB_DIRECTORY; use the checked-in LSF wrapper")
a<-c(a[1],"markers",a[-1])
# runtime.R reads commandArgs; source its loader and invoke the same tested phase dispatcher.
source(file.path(root,"R/config.R"));c<-read_config_v2(a[1])
if(dir.exists(c$tools$r_library)).libPaths(c$tools$r_library,include.site=FALSE)
for(k in c("common","workflow","parental","models","exports","acquisition","public","reference","pipeline","stages","scheduler","seadragon","lsf"))source(file.path(root,"R",paste0(k,".R")))
suppressPackageStartupMessages(library(GenomicRanges));source(file.path(root,"R/annotation.R"))
runtime_action(c,root,a[2],a[3],a[4],a[5])
