job_command <- function(c,root,stage) {
  cmd<-c(c$tools$Rscript,file.path(root,"poe.R"),attr(c,"file"),paste0("worker-",stage))
  if(isTRUE(c$execution$container_enabled)) {
    binds<-unlist(lapply(c$execution$binds,function(b)c("--bind",paste(b$host,b$container,b$mode,sep=":"))))
    cmd<-c(c$tools$singularity,"exec","--cleanenv","--env",
      paste0("R_LIBS_USER=",c$tools$r_library),binds,c$execution$image,cmd)
  }
  paste(vapply(cmd,shQuote,character(1)),collapse=" ")
}
submit_stage <- function(c,root,stage,submit=FALSE,dependency="none") {
  if(identical(c$execution$backend,"seadragon"))
    return(submit_graph(c,root,stage,submit,dependency))
  if(!stage %in% c("reference","pilot","scale","models"))stop("Invalid scheduler stage")
  if(!grepl("^(none|[0-9]+)$",dependency))stop("Numeric LSF dependency required")
  id<-paste0(format(Sys.time(),"%Y%m%dT%H%M%S"),"_",Sys.getpid())
  p<-function(k)cfg_path(c,k,stage=stage,job=id)
  script<-p("job_script");ensure_parent(script)
  lines<-c("#!/bin/bash","set -euo pipefail","umask 077",
    paste0("export R_LIBS_USER=",shQuote(c$tools$r_library)),job_command(c,root,stage))
  writeLines(lines,script)
  r<-c$lsf$stages[[stage]]
  args<-c("-q",c$lsf$queue,"-n",r$cpus,"-W",r$wall_time,"-M",r$memory_value,
    "-R",paste0("rusage[mem=",r$memory_value,"]"),"-J",paste(c$cohort_id,stage,id,sep="_"),
    "-oo",p("job_stdout"),"-eo",p("job_stderr"))
  if(dependency!="none")args<-c(args,"-w",paste0("done(",dependency,")"))
  if(!submit) {
    cat("DRY RUN ONLY\n",paste(shQuote(c$tools$bsub),paste(vapply(args,shQuote,character(1)),collapse=" "),
      "<",shQuote(script)),"\n");return(invisible(script))
  }
  if(!isTRUE(c$lsf$memory_units_verified)||!isTRUE(c$lsf$memory_per_job_verified)||
    !isTRUE(c$lsf$execution_environment_verified))stop("LSF memory/environment attestation required")
  preflight_v2(c,root,execution=TRUE)
  if(!file.exists(c$tools$bsub))stop("Configured bsub is not visible")
  if(isTRUE(c$execution$container_enabled)&&!file.exists(c$execution$image))stop("Container image missing")
  for(k in c("job_stdout","job_stderr"))ensure_parent(p(k))
  output<-system2(c$tools$bsub,vapply(args,shQuote,character(1)),stdin=script,stdout=TRUE,stderr=TRUE)
  if(!is.null(attr(output,"status")))stop("LSF submission failed: ",paste(output,collapse=" "))
  job<-regmatches(paste(output,collapse=" "),regexec("Job <([0-9]+)>",paste(output,collapse=" ")))[[1]]
  if(length(job)!=2)stop("LSF returned no numeric job ID; inspect scheduler before retrying")
  cat("Submitted LSF job ",job[2],"; dependency uses done(), so failed jobs cannot unlock downstream stages.\n",sep="")
  invisible(job[2])
}
