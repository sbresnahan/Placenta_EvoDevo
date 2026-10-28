#!/usr/bin/env Rscript
# No biological observations. Quantifier, bcftools, container and scheduler are MOCKS.
f<-gsub("~+~"," ",sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]),fixed=TRUE)
root<-normalizePath(file.path(dirname(f),".."),mustWork=TRUE)
for(k in c("common","workflow","config","parental","models","exports","acquisition","public",
           "reference","pipeline","stages","scheduler","seadragon","lsf"))
  source(file.path(root,"R",paste0(k,".R")))
source(file.path(root,"tests/helpers_seadragon.R"))
n<-0L;check<-function(x){stopifnot(isTRUE(x));n<<-n+1L}
fails<-function(x)inherits(tryCatch(force(x),error=identity),"error")
d<-tempfile("poe host boundary fixtures ");dir.create(d)
cat("Deterministic mocked-host fixture directory:",d,"\n")
host<-make_host_fixture(root,file.path(d,"mock host"))
oldenv<-Sys.getenv(c("SCRIPTS_DIR","POE_SITE_FILE","CONFIG","POE_SAMPLE","POE_NATIVE_RECORD"),unset=NA)
restore<-function()for(k in names(oldenv))if(is.na(oldenv[[k]]))Sys.unsetenv(k) else
  do.call(Sys.setenv,setNames(list(oldenv[[k]]),k))
Sys.setenv(SCRIPTS_DIR=root,POE_SITE_FILE=host$site)
Sys.unsetenv("POE_NATIVE_RECORD")
run_shell<-function(file,args=character(),okay=TRUE) {
  out<-suppressWarnings(system2("bash",c(shQuote(file),vapply(args,shQuote,character(1))),stdout=TRUE,stderr=TRUE))
  rc<-attr(out,"status");if(is.null(rc))rc<-0L
  if(okay&&rc!=0L)stop("Worker failed: ",file,"\n",paste(out,collapse="\n"))
  list(rc=rc,out=out)
}
allfiles<-list.files(root,recursive=TRUE,full.names=TRUE)
shell<-allfiles[grepl("\\.(sh|lsf)$|/bin/",allfiles)]
for(p in shell)check(system2("bash",c("-n",shQuote(p)))==0L)
check(length(lsf_scripts)==11L)
site<-readLines(file.path(root,"config/seadragon.sh"))
check(any(grepl("/software/ornaments/bin/activate",site,fixed=TRUE)))
check(any(grepl("bcftools-1.16",site,fixed=TRUE)))
check(!any(grepl("conda activate ornaments",site,fixed=TRUE)))
q<-plan_line("/bin/printf",c("%s","literal ; $HOME ' quoted"))
testplan<-write_plan(q,d)
out<-run_shell(testplan)
check(identical(out$out,"literal ; $HOME ' quoted"))
check(!file.exists(file.path(d,"quoted")))
# Scope one tiny public and human profile, two libraries each.
profiles<-lapply(c(FALSE,TRUE),function(h)make_profile_fixture(root,
  file.path(d,if(h)"human" else "public"),host,h))
for(p in profiles) {
  c<-p$c;Sys.setenv(CONFIG=p$config,POE_SAMPLE="shared")
  check(fails(native_gate(c)))
  out<-run_shell(file.path(root,"00_setup/preflight.sh"),c(p$config,"--execute"))
  check(file.exists(c$paths$native_inventory))
  check(is.data.table(native_gate(c)))
  check(fails(scale_gate(c,root,load_manifests(c))))
  # A changed native executable is rejected, not silently accepted/rebuilt.
  original<-readLines(host$ornaments);writeLines(c(original,"# changed mock executable"),host$ornaments)
  out<-run_shell(file.path(root,"bin/host_inventory.sh"),file.path(d,"changed.tsv"))
  Sys.setenv(POE_NATIVE_RECORD=file.path(d,"changed.tsv"))
  check(fails(native_gate(c)))
  writeLines(original,host$ornaments);Sys.unsetenv("POE_NATIVE_RECORD")
  # Inspect ordinary fan-out/dependencies without invoking a real scheduler.
  g<-graph_spec(c,"pilot")
  check(length(g)==3L&&g[[3]]$action=="collect-pilot")
  check(!any(vapply(g,function(x)x$action=="quantify",logical(1))))
  text<-capture.output(reg<-submit_graph(c,root,"pilot",FALSE,"123"))
  check(nrow(reg)==3L&&all(reg$state=="dry_run"))
  check(any(grepl("done(123)",text,fixed=TRUE)))
  check(grepl("done(900001) && done(900002)",reg$command[3],fixed=TRUE))
  check(all(grepl("-M' '8",reg$command,fixed=TRUE)))
  check(fails(submit_graph(c,root,"pilot",FALSE,"12;bad")))
  for(action in c("reference",if(c$profile=="public_cross")c("markers","index"))) {
    Sys.setenv(POE_SAMPLE="shared")
    out<-run_shell(file.path(root,lsf_scripts[[action]]))
  }
  if(c$profile=="human_trio")for(s in p$samples)for(action in c("child","index")) {
    Sys.setenv(POE_SAMPLE=s);out<-run_shell(file.path(root,lsf_scripts[[action]]))
  }
  check(!fails(reference_gate(c,root)))
  if(c$profile=="human_trio") {
    check(!fails(child_index_gate(c,root,p$samples[1])))
    typed<-fread(cfg_path(c,"child_query",sample=p$samples[1]),header=FALSE)
    check(ncol(typed)==7L&&typed$V6[1]=="0/0"&&typed$V7[1]=="1/1")
  }
  for(s in p$samples) {
    Sys.setenv(POE_SAMPLE=s)
    # Manual scale worker cannot bypass the missing report authorization.
    out<-run_shell(file.path(root,lsf_scripts[["quantify"]]),okay=FALSE)
    check(out$rc!=0L&&any(grepl("Explicit scale authorization missing",out$out,fixed=TRUE)))
    before<-proc.time()[["elapsed"]]
    out<-run_shell(file.path(root,lsf_scripts[["pilot"]]))
    cat(c$cohort_id,s,"mock pilot wall seconds:",round(proc.time()[["elapsed"]]-before,2),"\n")
    check(file.exists(library_pilot_report(c,s)))
  }
  check(!file.exists(c$paths$scale_receipt))
  Sys.setenv(POE_SAMPLE="shared")
  out<-run_shell(file.path(root,lsf_scripts[["collect-pilot"]]))
  check(file.exists(c$paths$pilot_report))
  report<-fread(c$paths$pilot_report)
  check(nrow(report)==8L&&all(report$elapsed_seconds==1))
  check(fails(authorize_scale(c,root,strrep("0",64))))
  authorize_scale(c,root,sha256(c$paths$pilot_report))
  check(!fails(scale_gate(c,root,load_manifests(c))))
  Sys.setenv(POE_SAMPLE=p$samples[1],POE_FIXTURE_FAIL="yes")
  out<-run_shell(file.path(root,lsf_scripts[["quantify"]]),okay=FALSE)
  check(out$rc!=0L)
  check(!file.exists(cfg_path(c,"quant_receipt",sample=p$samples[1],mode="complete",strand="fr")))
  Sys.unsetenv("POE_FIXTURE_FAIL")
  for(s in p$samples) {
    Sys.setenv(POE_SAMPLE=s)
    out<-run_shell(file.path(root,lsf_scripts[["quantify"]]))
    out<-run_shell(file.path(root,lsf_scripts[["assign"]]))
  }
  Sys.setenv(POE_SAMPLE="shared")
  out<-run_shell(file.path(root,lsf_scripts[["collect-quantification"]]))
  out<-run_shell(file.path(root,lsf_scripts[["models"]]))
  check(file.exists(c$paths$model_promoter)&&file.exists(c$paths$model_receipt))
  counts<-completed_counts(c,root,load_manifests(c))
  check(uniqueN(counts$event$sample_id)==2L)
  check(all(counts$event$maternal_count+counts$event$paternal_count==20))
  if(c$profile=="human_trio") {
    check(all(counts$event[sample_id=="fixtureA"]$maternal_count==12))
    check(all(counts$event[sample_id=="fixtureB"]$maternal_count==8))
    check(all(fread(c$paths$model_promoter)$status!="eligible_conditional_inference"))
    check(fails(export_v2(c)))
    typed<-cfg_path(c,"child_query",sample=p$samples[1]);original<-readLines(typed)
    writeLines("CHANGED_TYPED_ANCHOR_FIXTURE",typed)
    check(fails(completed_counts(c,root,load_manifests(c))))
    writeLines(original,typed)
  }
  # Duplicate/stale lock blocks the real checked-in worker before any native work.
  lock<-file.path(c$paths$receipts,"locks","models_shared")
  dir.create(lock,recursive=TRUE);writeLines("DETERMINISTIC_LOCK_FIXTURE",file.path(lock,"owner"))
  out<-run_shell(file.path(root,lsf_scripts[["models"]]),okay=FALSE)
  check(out$rc!=0L&&any(grepl("Stage/library locked",out$out,fixed=TRUE)))
  check(file.exists(file.path(lock,"owner")));unlink(lock,recursive=TRUE)
  # An edited output invalidates completion; no automatic reauthorization.
  sample<-p$samples[1]
  output<-cfg_path(c,"event_counts",sample=sample,mode="complete",strand="fr")
  writeLines("CHANGED_COORDINATE_FIXTURE",output)
  check(fails(completed_counts(c,root,load_manifests(c))))
}
# Test scheduler success/failure with a dedicated executable mock; NO bsub installed/called.
p<-profiles[[1]];c<-p$c
fake<-file.path(d,"fake_bsub")
counter<-file.path(d,"scheduler_counter")
writeLines(c("#!/bin/bash","set -euo pipefail",
  "read -r first; [[ \"$first\" == '#!/bin/bash' ]] || exit 4",
  paste0("counter=",shQuote(counter)),
  'n=12344; [[ ! -f "$counter" ]] || read -r n < "$counter"; n=$((n+1))',
  'printf "%s\\n" "$n" > "$counter"',
  'printf "%s\\n" "$@" > "${counter}.${n}.args"',
  "printf 'Job <%s> is submitted to queue <fixture>.\\n' \"$n\""),fake);Sys.chmod(fake,"0755")
c$tools$bsub<-fake
# Only acquisition needs no completed-reference gate; remains a no-download scheduler mock.
c$retrieval$internet_host_approved<-TRUE
text<-capture.output(reg<-submit_graph(c,root,"retrieve",TRUE))
check(reg$state=="host_planned"&&!file.exists(counter))
plan<-sub("^HOST_PLAN=","",grep("^HOST_PLAN=",text,value=TRUE))
check(length(plan)==1L&&file.exists(plan))
out<-run_shell(plan)
submitted<-paste0(sub("\\.host\\.sh$","",plan),".submitted.tsv")
check(fread(submitted)$job_id==12345L)
check(run_shell(plan,okay=FALSE)$rc!=0L) # one-shot plan cannot duplicate a submission
text<-capture.output(reg<-submit_graph(c,root,"reference",TRUE))
plan<-sub("^HOST_PLAN=","",grep("^HOST_PLAN=",text,value=TRUE))
out<-run_shell(plan)
args_second<-readLines(paste0(counter,".12347.args"))
check("done(12346)" %in% args_second)
args_third<-readLines(paste0(counter,".12348.args"))
check("done(12347)" %in% args_third)
writeLines(c("#!/bin/bash","printf 'invalid scheduler response\\n'"),fake)
text<-capture.output(reg<-submit_graph(c,root,"retrieve",TRUE))
plan<-sub("^HOST_PLAN=","",grep("^HOST_PLAN=",text,value=TRUE))
check(run_shell(plan,okay=FALSE)$rc!=0L)
writeLines(c("#!/bin/bash","exit 5"),fake)
Sys.sleep(1) # ensure a distinct single-use plan filename under the same PID
text<-capture.output(reg<-submit_graph(c,root,"retrieve",TRUE))
plan<-sub("^HOST_PLAN=","",grep("^HOST_PLAN=",text,value=TRUE))
check(run_shell(plan,okay=FALSE)$rc!=0L)
restore()
unlink(d,recursive=TRUE)
cat("PASS:",n,"Seadragon shell/host-boundary/fan-out/mock integration assertions\n")
cat("All native counts, resource metrics, bcftools operations, scheduler and container behavior above are MOCK fixtures.\n")
