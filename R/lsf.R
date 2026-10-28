# LSF workers are the Seadragon submission interface.
lsf_scripts <- c(retrieve="01_retrieval/01_retrieve_public.lsf",
 reference="02_reference/01_prepare_reference.lsf",markers="03_genotypes/01_prepare_markers.lsf",
 child="03_genotypes/02_prepare_child.lsf",index="04_ornaments/01_build_index.lsf",
 pilot="04_ornaments/02_run_pilot.lsf",quantify="04_ornaments/03_quantify_library.lsf",
 `collect-pilot`="04_ornaments/10_collect_pilot.lsf",
 `collect-quantification`="04_ornaments/11_collect_quantification.lsf",
 assign="05_parental_expression/01_assign_parents.lsf",models="06_promoter_models/01_run_models.lsf")
graph_spec <- function(c,stage) {
  m<-load_manifests(c);g<-list()
  add<-function(key,action,sample="shared",deps=character()) {
    g[[length(g)+1L]]<<-list(key=key,action=action,sample=sample,deps=deps)
  }
  if(stage=="retrieve") {
    if(c$profile!="public_cross")stop("Public retrieval is not a human-data access route")
    add("retrieve","retrieve")
  } else if(stage=="reference") {
    add("reference","reference")
    if(c$profile=="public_cross") {
      add("markers","markers",deps="reference");add("index","index",deps="markers")
    } else for(s in m$libraries$sample_id) {
      add(paste0("child_",s),"child",s,"reference")
      add(paste0("index_",s),"index",s,paste0("child_",s))
    }
  } else if(stage=="pilot") {
    for(s in m$libraries$sample_id)add(paste0("pilot_",s),"pilot",s)
    add("collect-pilot","collect-pilot",deps=vapply(g,`[[`,character(1),"key"))
  } else if(stage=="scale") {
    for(s in m$libraries$sample_id) {
      add(paste0("quantify_",s),"quantify",s)
      add(paste0("assign_",s),"assign",s,paste0("quantify_",s))
    }
    add("collect-quantification","collect-quantification",
      deps=vapply(g[vapply(g,function(x)x$action=="assign",logical(1))],`[[`,character(1),"key"))
  } else if(stage=="models")add("models","models")
  else stop("Unknown submission stage")
  g
}
stage_resource <- function(c,action) {
  key<-if(action %in% c("quantify","assign","collect-quantification"))"scale" else
    if(action %in% c("markers","child","index"))"reference" else
    if(action=="collect-pilot")"pilot" else action
  r<-c$lsf$stages[[key]]
  if(is.null(r))stop("Missing resources for ",action)
  if(is.null(r$reservation_value))r$reservation_value<-r$memory_value
  if(is.null(r$queue))r$queue<-c$lsf$queue
  for(k in c("memory_value","reservation_value","cpus"))positive(r[[k]],paste(action,k))
  if(!grepl("^[A-Za-z0-9_.-]+$",r$queue)||!grepl("^[0-9]+:[0-5][0-9]$",r$wall_time))
    stop("Unsafe LSF resource value")
  if(action %in% c("pilot","quantify")&&r$cpus<c$ornaments$threads)stop("LSF CPUs below Ornaments threads")
  r
}
submit_graph <- function(c,root,stage,submit=FALSE,dependency="none") {
  if(!grepl("^(none|[0-9]+)$",dependency))stop("Numeric LSF dependency required")
  g<-graph_spec(c,stage)
  if(submit) {
    validate_inputs_v2(c,root,execution=TRUE)
    if(stage=="pilot")reference_gate(c,root)
    if(stage %in% c("scale","models"))scale_gate(c,root,load_manifests(c))
  }
  id<-paste0(format(Sys.time(),"%Y%m%dT%H%M%OS6"),"_",Sys.getpid(),"_",basename(tempfile("r")))
  id<-gsub("[^A-Za-z0-9_.-]","_",id)
  registry<-sub("\\.tsv$",paste0("_",id,".tsv"),c$paths$job_registry)
  if(identical(registry,c$paths$job_registry))registry<-paste0(registry,"_",id,".tsv")
  jobs<-character();rows<-list();hostlines<-character()
  hostplan<-paste0(registry,".host.sh");submitted<-paste0(registry,".submitted.tsv")
  if(submit) {
    ensure_parent(hostplan)
    inventory<-paste0(registry,".native.tsv")
    hostlines<-c("#!/bin/bash","set -eo pipefail","umask 077",
      paste("export SCRIPTS_DIR=",shQuote(root),sep=""),
      paste("export POE_SITE_FILE=",shQuote(c$execution$site_file),sep=""),
      paste("export POE_NATIVE_RECORD=",shQuote(inventory),sep=""),
      paste("mkdir",shQuote(paste0(hostplan,".started"))),
      paste(shQuote(file.path(root,"bin/host_inventory.sh")),shQuote(inventory)),
      plan_line(file.path(root,"bin/Rscript_sif"),c(file.path(root,"runtime.R"),attr(c,"file"),
        "preflight","finish","shared","unused")),
      "declare -A jobs seen",
      "rows=$'node_key\\taction\\tsample\\tjob_id\\tdepends_on\\tstate\\n'")
  }
  for(i in seq_along(g)) {
    node<-g[[i]];r<-stage_resource(c,node$action)
    script<-file.path(root,lsf_scripts[[node$action]])
    if(!file.exists(script))stop("Missing LSF worker: ",script)
    deps<-unname(jobs[node$deps])
    if(!length(node$deps)&&dependency!="none")deps<-dependency
    jobname<-paste(c$cohort_id,node$action,node$sample,id,sep="_")
    if(nchar(jobname)>180)stop("LSF job name too long")
    logs<-file.path(c$paths$logs,"lsf");dir.create(logs,recursive=TRUE,showWarnings=FALSE)
    # LSF -env parses comma separators; reject rather than misparse values.
    env<-c(CONFIG=attr(c,"file"),SCRIPTS_DIR=root,POE_SAMPLE=node$sample,
      POE_SITE_FILE=c$execution$site_file)
    if(any(grepl("[,\r\n\"']",env)))stop("Path cannot be represented safely in LSF -env")
    # -env is already shell-quoted below. Adding quotes around individual
    # values passes literal quote characters into the LSF job environment.
    envarg<-paste0("all,",paste(paste0(names(env),"=",env),collapse=","))
    args<-c("-q",r$queue,"-n",r$cpus,"-M",r$memory_value,
      "-R",paste0("rusage[mem=",r$reservation_value,"]"),"-W",r$wall_time,
      "-J",jobname,"-o",file.path(logs,paste0(jobname,".%J.out")),
      "-e",file.path(logs,paste0(jobname,".%J.err")),"-env",envarg)
    baseargs<-args
    if(length(deps))args<-c(args,"-w",paste(paste0("done(",deps,")"),collapse=" && "))
    command<-paste(shQuote(c$tools$bsub),paste(vapply(as.character(args),shQuote,character(1)),collapse=" "),
      "<",shQuote(script))
    jid<-as.character(900000L+i) # dry-run/planning labels ONLY; never scheduler results
    if(!submit) {
      cat("DRY RUN ONLY:",command,"\n")
    } else {
      dep_expr<-if(length(node$deps))paste(paste0("done(${jobs[",node$deps,"]})"),collapse=" && ") else
        if(dependency!="none")paste0("done(",dependency,")") else ""
      hostcommand<-paste(shQuote(c$tools$bsub),
        paste(vapply(as.character(baseargs),shQuote,character(1)),collapse=" "),
        if(nzchar(dep_expr))paste0('-w "',dep_expr,'"') else "",
        "<",shQuote(script))
      hostlines<-c(hostlines,
        paste0("if ! output=$(",hostcommand,"); then printf 'LSF submission failed; inspect already submitted jobs\\n' >&2; exit 1; fi"),
        'if [[ "$output" =~ Job\\ \\<([0-9]+)\\> ]]; then jid="${BASH_REMATCH[1]}"; else printf "Missing numeric LSF ID; inspect scheduler before retrying\\n" >&2; exit 1; fi',
        '[[ -z "${seen[$jid]:-}" ]] || { printf "Duplicate scheduler ID\\n" >&2; exit 1; }',
        'seen[$jid]=yes',paste0('jobs[',node$key,']="$jid"'),
        paste0("rows+=$(printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s' ",shQuote(node$key)," ",
          shQuote(node$action)," ",shQuote(node$sample),' "$jid" "',dep_expr,'" submitted)'),
        "rows+=$'\\n'",
        paste0('printf "%s" "$rows" > ',shQuote(paste0(submitted,".tmp"))),
        paste("mv -f",shQuote(paste0(submitted,".tmp")),shQuote(submitted)),
        paste0('printf "Submitted %s as %s\\n" ',shQuote(node$key),' "$jid"'))
    }
    jobs[node$key]<-jid
    rows[[i]]<-data.table(node_key=node$key,action=node$action,sample=node$sample,
      job_id=paste0(if(submit)"HOST_PLANNED_" else "DRYRUN_",jid),depends_on=paste(deps,collapse=","),
      script=lsf_scripts[[node$action]],state=if(submit)"host_planned" else "dry_run",command=command)
    write_tsv(rbindlist(rows),registry)
  }
  cat("Job registry:",registry,"\n")
  if(submit) {
    writeLines(hostlines,hostplan)
    cat("Prepared host-only submission plan; R has not called bsub.\n")
    cat("HOST_PLAN=",hostplan,"\n",sep="")
  }
  invisible(rbindlist(rows))
}
