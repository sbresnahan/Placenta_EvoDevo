#!/usr/bin/env Rscript
# Relocate configuration/manifests and change current directory. Never invoke bsub.
f<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
f<-gsub("~+~"," ",f,fixed=TRUE)
root<-normalizePath(file.path(dirname(f),".."))
d<-tempfile("relocated config with spaces ");dir.create(d)
dir.create(file.path(d,"config"));dir.create(file.path(d,"manifests"))
cfg<-file.path(d,"config","mouse.example.yaml")
stopifnot(file.copy(file.path(root,"config/mouse_Babr_CAST.example.yaml"),cfg))
stopifnot(all(file.copy(list.files(file.path(root,"manifests"),full.names=TRUE),file.path(d,"manifests"))))
gusto<-file.path(d,"config","human.example.yaml")
stopifnot(file.copy(file.path(root,"config/GUSTO.example.yaml"),gusto))
for(p in c(cfg,gusto)) {
  x<-yaml::read_yaml(p)
  # Deliberate test-only path relocation; never write to production site roots.
  for(section in c("paths","templates"))
    x[[section]]<-lapply(x[[section]],function(v)
      if(startsWith(v,"/rsrch"))file.path(d,"local",sub("^/","",v)) else v)
  yaml::write_yaml(x,p)
}
old<-getwd();setwd(d)
rs<-file.path(R.home("bin"),"Rscript")
run<-function(config,stage,extra=character())suppressWarnings(system2(rs,
 c(shQuote(file.path(root,"poe.R")),shQuote(config),stage,vapply(extra,shQuote,character(1))),
 stdout=TRUE,stderr=TRUE))
n<-0L;check<-function(x){stopifnot(isTRUE(x));n<<-n+1L}
x<-run(cfg,"preflight");check(is.null(attr(x,"status"))&&any(grepl("valid;",x,fixed=TRUE)))
x<-run(cfg,"pilot",c("--after","123"))
if(!is.null(attr(x,"status")))cat(paste(x,collapse="\n"),"\n")
check(is.null(attr(x,"status"))&&any(grepl("DRY RUN ONLY",x,fixed=TRUE)))
check(any(grepl("done(123)",x,fixed=TRUE)))
jobs<-list.files(file.path(d,"local"),"jobs_.*\\.tsv$",recursive=TRUE,full.names=TRUE)
check(length(jobs)==1L)
check(nrow(data.table::fread(jobs[1]))==6L)
x<-run(cfg,"scale",c("--submit","--accept",strrep("0",64)))
check(!is.null(attr(x,"status")))
x<-run(gusto,"export")
check(!is.null(attr(x,"status"))&&any(grepl("Release disabled",x,fixed=TRUE)))
x<-run(cfg,"pilot",c("--after","123;bad"))
check(!is.null(attr(x,"status"))&&any(grepl("Numeric LSF dependency",x,fixed=TRUE)))
setwd(old);unlink(d,recursive=TRUE)
cat("PASS:",n,"relocated-config/CLI/quoting/fail-closed assertions\n")
