#!/usr/bin/env Rscript
f<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
f<-gsub("~+~"," ",f,fixed=TRUE)
root<-normalizePath(file.path(dirname(f),".."))
for(t in c("test_libpaths.R","test_core.R","test_workflow.R","test_v2.R","test_cli.R","test_seadragon.R")) {
  rc<-system2(file.path(R.home("bin"),"Rscript"),shQuote(file.path(root,"tests",t)))
  if(rc!=0L)stop("Test failed: ",t)
}
cat("ALL LOCAL LOGICAL TEST SUITES PASSED. No cohort execution or live LSF validation.\n")
