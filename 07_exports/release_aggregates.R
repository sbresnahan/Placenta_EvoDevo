#!/usr/bin/env Rscript
f<-gsub("~+~"," ",sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]),fixed=TRUE)
root<-normalizePath(file.path(dirname(f),".."),mustWork=TRUE)
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=1L)stop("One explicit private CONFIG required")
source(file.path(root,"R/config.R"));c<-read_config_v2(a[1])
if(dir.exists(c$tools$r_library)).libPaths(c$tools$r_library,include.site=FALSE)
for(k in c("common","workflow","exports","stages"))source(file.path(root,"R",paste0(k,".R")))
export_v2(c)
