#!/usr/bin/env Rscript
# Internal worker phase interface, always called through the R-container bridge.
f<-gsub("~+~"," ",sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]),fixed=TRUE)
root<-normalizePath(dirname(f),mustWork=TRUE)
source(file.path(root,"R/library_paths.R"))
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=5L)stop("Usage: runtime.R CONFIG STAGE PHASE LIBRARY JOB_DIRECTORY")
source(file.path(root,"R/config.R"));c<-read_config_v2(a[1],setup=a[2]=="preflight")
if(dir.exists(c$tools$r_library)).libPaths(c$tools$r_library,include.site=FALSE)
for(k in c("common","workflow","parental","models","exports","acquisition","public","reference",
           "pipeline","stages","scheduler","seadragon","lsf"))source(file.path(root,"R",paste0(k,".R")))
suppressPackageStartupMessages(library(GenomicRanges))
source(file.path(root,"R/annotation.R"))
data.table::setDTthreads(if(a[2]=="preflight")1L else as.integer(c$ornaments$threads))
if(!is_seadragon(c))stop("This worker requires execution.backend: seadragon")
runtime_action(c,root,a[2],a[3],a[4],a[5])
