#!/usr/bin/env Rscript
f<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
f<-gsub("~+~"," ",f,fixed=TRUE)
root<-normalizePath(file.path(dirname(f),".."))
source(file.path(root,"R/common.R"));source(file.path(root,"R/workflow.R"))
source(file.path(root,"R/config.R"))
source(file.path(root,"R/exports.R"));source(file.path(root,"R/models.R"))
n<-0L
check<-function(x){stopifnot(isTRUE(x));n<<-n+1L}
fails<-function(expr)inherits(tryCatch(force(expr),error=identity),"error")
d<-tempfile();dir.create(d)
i<-file.path(d,"input");o<-file.path(d,"output");r<-file.path(d,"receipt")
writeLines("logical fixture",i);writeLines("completed",o)
check(!stage_valid(r,o))
stage_finish(r,o);check(stage_valid(r,o))
writeLines("changed",o);check(stage_valid(r,o))
writeLines("changed input",i);check(stage_valid(r,o))
stage_begin(r);check(!stage_valid(r,o))
stage_finish(r,o);unlink(o);check(!stage_valid(r,o))
writeLines("read fixture",file.path(d,"fq"))
s<-data.table(sample_id="logic",child_id="C",layout="SE",strand_mode="unstranded",
 fragment_mean=1,fragment_sd=1,fastq_1=file.path(d,"fq"),
 fastq_2=NA_character_)
t<-data.table(child_id="C",mother_id="M",father_id="F",family_id="family")
check(nrow(check_linkage(s,t))==1)
check(fails(check_linkage(rbind(s,s),t)))
check(nrow(check_linkage(copy(s)[,identity_verified:=FALSE],t))==1) # Legacy extra column is ignored.
check(fails(check_linkage(s,copy(t)[,mother_id:="C"])))
check(fails(check_linkage(copy(s)[,strand_mode:="unknown"],t)))
check(fails(check_linkage(copy(s)[,fragment_sd:=NA_real_],t)))
check(fails(validate_aggregate(data.table(participant_id="not_exportable"),
 data.table(),10)))
check(fails(read_config(file.path(root,"config/maize_B73_Ki11.example.yaml")))==FALSE)
check(nrow(fit_all_pairs(data.table()))==0L)
# Algebraic HC3 test, not a biological data generator.
X<-cbind(1,seq(-1,1,length.out=20));y<-(1:20)^2
check(abs(hc3_intercept(y,X)$estimate-mean(y))<1e-10)
# Nonbiological schema constants: no donor observations or discovery estimates.
public<-data.table(gene_id="REFERENCE_GENE",gene_name="REFERENCE_SYMBOL",
 event_a="REFERENCE_EVENT_A",event_b="REFERENCE_EVENT_B",pair_type="TE_nonTE")
safe<-copy(public)
safe[,`:=`(n_families=10L,status="eligible_conditional_inference",estimate=0,
 se=1,p_value=.5,q_value=.5,ci_low=-1,ci_high=1,residual_df=5L,
 bootstrap_valid_fraction=1,cis_basis="two_cis_terms",
 uncertainty="between_families_conditional_on_Ornaments",discordant=FALSE)]
check(validate_aggregate(safe,public,10))
check(fails(validate_aggregate(copy(safe)[,participant_id:="FORBIDDEN"],public,10)))
check(fails(validate_aggregate(safe,public,11)))
check(fails(validate_aggregate(copy(safe)[,uncertainty:="FORBIDDEN_FREE_TEXT"],public,10)))
check(fails(validate_aggregate(copy(safe)[,event_a:="NOT_IN_PUBLIC_REFERENCE"],public,10)))
unlink(d,recursive=TRUE)
cat("PASS:",n,"workflow/privacy/algebra assertions\n")
