#!/usr/bin/env Rscript
# All observations here are deterministic schema/coordinate/algebra fixtures, not biological data.
f<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
f<-gsub("~+~"," ",f,fixed=TRUE)
root<-normalizePath(file.path(dirname(f),".."))
for(k in c("common","workflow","config","annotation","parental","models","exports","acquisition",
  "public","reference","pipeline","stages","scheduler"))source(file.path(root,"R",paste0(k,".R")))
n<-0L
check<-function(x){stopifnot(isTRUE(x));n<<-n+1L}
fails<-function(x)inherits(tryCatch(force(x),error=identity),"error")
temp<-tempfile("poe logical fixtures ");dir.create(temp)
c<-read_config_v2(file.path(root,"config/mouse_Babr_CAST.example.yaml"))
check(c$markers$proxy&&c$parents$focal=="C57BL_6J_Babr")
m<-load_manifests(c)
check(nrow(m$libraries)==5&&all(m$libraries$fragment_mean_nt==250)&&all(m$libraries$fragment_sd_nt==50))
check(!inherits(m$sources$bytes_1,"integer64"))
check(!any(grepl("grid",m$libraries$fragment_parameter_origin)))
check(fails(resolve_value("${UNKNOWN}/x",temp,character())))
check(fails(resolve_value("$HOME/x",temp,character())))
check(resolve_value("relative space/x",temp)==file.path(temp,"relative space/x"))
check(fails(cfg_path(c,"fastq",run="bad;command",mate=1)))
check(fails(cfg_path(c,"fastq",run="okay")))
for(id in c("maize_B73_Ki11","maize_B73_W22","mouse_FVB_CAST_E12")) {
  cc<-read_config_v2(file.path(root,"config",paste0(id,".example.yaml")))
  mm<-load_manifests(cc)
  check(all(mm$libraries$layout=="PE"))
  check(all(is.na(mm$libraries$fragment_mean_nt)))
  check(is.list(validate_inputs_v2(cc,root)))
}
cc<-read_config_v2(file.path(root,"config/GUSTO.example.yaml"))
check(fails(load_manifests(cc)))
check(fails(export_v2(cc)))
# Config relocation, strict profile/threshold types and template tokens.
raw<-yaml::read_yaml(file.path(root,"config/mouse_Babr_CAST.example.yaml"))
for(change in list(list(schema_version=1L),list(profile="human_trio"),list(species="unknown"))) {
  bad<-modifyList(raw,change);p<-file.path(temp,"bad.yaml");yaml::write_yaml(bad,p)
  check(fails(read_config_v2(p)))
}
bad<-raw;bad$templates$prefix<-"bad/{unknown}";yaml::write_yaml(bad,p)
check(fails(read_config_v2(p)))
bad<-raw;bad$models$min_depth<-19;yaml::write_yaml(bad,p)
check(fails(read_config_v2(p)))
bad<-raw;bad$fragments$fragment_sd_nt<-0;yaml::write_yaml(bad,p)
check(fails(read_config_v2(p)))
raw$paths$work<-"../local/relocated-work"
yaml::write_yaml(raw,p);relocated<-read_config_v2(p)
check(startsWith(relocated$paths$work,temp))
# HTTP response truth table. No network traffic or mock download success claims.
header<-charToRaw("HTTP/1.1 206 Partial Content\r\nContent-Range: bytes 10-19/20\r\n")
check(range_ok(206L,header,10,10,20))
check(!range_ok(200L,header,10,10,20))
check(!range_ok(206L,header,9,10,20))
check(!range_ok(206L,charToRaw("Content-Range: bytes 10-19/21"),10,10,20))
check(range_ok(200L,raw(),0,20,20))
empty<-file.path(temp,"empty");file.create(empty)
check(!verify_file(empty,1,strrep("0",32)))
dat<-file.path(temp,"checksum");writeLines("logical",dat)
check(verify_file(dat,file.info(dat)$size,unname(tools::md5sum(dat))))
check(!verify_file(dat,file.info(dat)$size+1,unname(tools::md5sum(dat))))
dest<-file.path(temp,"download");part<-file.path(temp,"partial");resp<-file.path(temp,"response")
writeBin(charToRaw("abcdef"),dest);md5<-unname(tools::md5sum(dest));unlink(dest)
writeBin(charToRaw("abc"),part)
download_cfg<-list(retrieval=list(attempts=1,timeout_seconds=1,connect_timeout_seconds=1,retry_seconds=0))
transport206<-function(url,path,handle) {
  writeBin(charToRaw("def"),path)
  list(status_code=206L,headers=charToRaw("Content-Range: bytes 3-5/6\r\n"))
}
fetch_one("https://example.invalid/fixture",dest,part,resp,6,md5,download_cfg,transport206)
check(verify_file(dest,6,md5))
unlink(dest);writeBin(charToRaw("abc"),part)
transport200<-function(url,path,handle) {
  writeBin(charToRaw("abcdef"),path);list(status_code=200L,headers=raw())
}
check(fails(fetch_one("https://example.invalid/fixture",dest,part,resp,6,md5,download_cfg,transport200)))
check(file.info(part)$size==3&&!file.exists(dest))
unlink(part)
check(fails(fetch_one("https://example.invalid/fixture",dest,part,resp,6,strrep("0",32),
  download_cfg,transport200)))
# The first maize FASTQ exceeds the 32-bit CURLOPT_MAXFILESIZE limit.
# Confirm curl handle setup reaches a mock transport for a >4 GiB file.
large_seen<-FALSE
large_transport<-function(url,path,handle) {
  large_seen<<-TRUE
  stop("Intentional offline transport failure")
}
check(fails(fetch_one("https://example.invalid/large",dest,part,resp,4359512380,
  strrep("0",32),download_cfg,large_transport)))
check(large_seen)
# FASTQ synchronization, nonrandom prefix and malformed/empty rejection.
fq<-function(path,id,seq="AC") {con<-gzfile(path,"wt");writeLines(c(paste0("@",id),seq,"+",
  strrep("!",nchar(seq))),con);close(con)}
f1<-file.path(temp,"1.gz");f2<-file.path(temp,"2.gz")
fq(f1,"truth/1");fq(f2,"truth/2")
check(validate_fastq(f1,f2)==1)
pref<-file.path(temp,"prefix.gz");check(validate_fastq(f1,limit=1,output1=pref)==1)
fq(f2,"mismatch/2");check(fails(validate_fastq(f1,f2)))
con<-gzfile(f2,"wt");writeLines(c("@x","AC","+","!"),con);close(con)
check(fails(validate_fastq(f2)))
check(fails(validate_fastq(empty)))
badgzip<-file.path(temp,"corrupt.gz");writeBin(charToRaw("not a gzip stream"),badgzip)
check(fails(run_checked(Sys.which("gzip"),c("-t",badgzip),stdout=FALSE,stderr=FALSE)))
# GT confidence: another strain's site FILTER never substitutes for parental FI.
check(parent_call(c(GT="0/0",FI="1"),"mgp")==0L)
check(is.na(parent_call(c(GT="0/0",FI="0"),"mgp")))
check(is.na(parent_call(c(GT="./.",FI="1"),"mgp")))
check(parent_call(c(GT="1/1",AD="0,8"),"panzea")==1L)
check(is.na(parent_call(c(GT="1/1",AD="3,5"),"panzea")))
check(is.na(parent_call(c(GT="0/1",AD="5,5"),"panzea")))
check(identical(parent_calls(c("0/0:1","1/1:0","./.:1"),rep("GT:FI",3),"mgp"),c(0L,NA_integer_,NA_integer_)))
check(identical(parent_calls(c("0/0:8,0","1/1:0,8","1/1:2,6"),rep("GT:AD",3),"panzea"),
  c(0L,1L,NA_integer_)))
vpath<-file.path(temp,"markers.vcf.gz");con<-gzfile(vpath,"wt")
writeLines(c("##fileformat=VCFv4.2",
 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFVB_NJ\tCAST_EiJ",
 "1\t10\t.\tA\tG\t.\tFAIL_OTHER_STRAIN\t.\tGT:FI\t0/0:1\t1/1:1",
 "1\t20\t.\tA\tG\t.\tPASS\t.\tGT:FI\t./.:1\t1/1:1"),con);close(con)
cc<-c;cc$markers$proxy<-FALSE;cc$markers$focal_sample<-"FVB_NJ"
vm<-vcf_markers(vpath,cc)
check(nrow(vm)==1&&vm$genomic_pos==10)
# Generic parental labels, conservation and missing-not-zero.
c$parents<-list(focal="FOCAL",alternate="OTHER")
x<-data.table(chrom="chrT",genomic_pos=c(1L,2L),ref="A",alt="G",ref_count=c(2.5,0),alt_count=c(7.5,1))
mk<-data.table(chrom="chrT",genomic_pos=1L,ref="A",alt="G",focal_index=1L,alternate_index=0L)
y<-orient_public(x,mk,data.table(maternal_parent="FOCAL"),c)
check(y$maternal_count[1]==7.5&&y$paternal_count[1]==2.5&&is.na(y$maternal_count[2]))
yy<-orient_public(x,mk,data.table(maternal_parent="OTHER"),c)
check(yy$maternal_count[1]==2.5&&sum(y[1,.(maternal_count,paternal_count)])==10)
# Fixed-anchor grouping, no 0->40->80 chaining.
tx<-data.table(transcript_id=c("a","b","c"),gene_id="T",chrom="chrT",strand="+",
 n_exons=2L,first_donor=100L,tss=c(1L,41L,81L))
gr<-group_promoters(tx,50)
check(uniqueN(gr$event_id)==2L)
# Algebraic rank/df behavior: not biological observations or a power simulation.
check(public_fit(c(1,2),c(1,-1))$status=="nonpositive_df_or_rank_deficient")
check(public_fit(1:4,rep(1,4))$status=="single_direction_cis_POE_confounded")
check(public_fit(1:6,c(1,1,1,-1,-1,-1))$residual_df==4)
check(public_fit(1:6,c(1,1,1,-1,-1,-1),FALSE)$status=="independence_unresolved")
# Minimum-total-depth sentinel, never maternal fraction; no missing units.
q<-data.table(unit_id=rep(c("u1","u2"),2),gene_id="T",locus_id=rep(c("A","B"),each=2),
 maternal_count=c(0,20,15,15),paternal_count=c(20,0,15,15))
frozen<-freeze_sentinels(q,"gene_id",c("u1","u2"),20)
check(all(frozen$locus_id=="B"))
check(nrow(freeze_sentinels(q[unit_id=="u1"],"gene_id",c("u1","u2")))==0)
check(fails(freeze_sentinels(rbind(q,q[1]),"gene_id",c("u1","u2"))))
u<-data.table(unit_id=c("u1","u2"),cluster_id=c("u1","u2"),
  maternal_parent=c("FOCAL","OTHER"))
e<-rbind(copy(q)[,event_id:="E1"],copy(q)[,event_id:="E2"])
endpoints<-public_endpoints(q,e,u,c)
check(nrow(endpoints$paired)==2L)
check(endpoints$gene$status=="nonpositive_df_or_rank_deficient")
u[,cluster_id:="same"];endpoints<-public_endpoints(q,e,u,c)
check(endpoints$promoter$status=="independence_unresolved")
# Negative-strand REF mapping, checked in existing core truth table; forward command contracts here.
l<-m$libraries[1];qargs<-quant_args(c,l,"INDEX","VCF","GENOTYPE","READ","OUT")
check("--single" %in% qargs && qargs[match("-l",qargs)+1]=="250")
l$layout<-"PE";qargs<-quant_args(c,l,"INDEX","VCF","GENOTYPE",c("R1","R2"),"OUT")
check(!"--single" %in% qargs && !"-l" %in% qargs)
# Completion records do not determine whether inputs or outputs are current.
i<-file.path(temp,"i");o<-file.path(temp,"o");r<-file.path(temp,"r.rds")
writeLines("truth",i);writeLines("truth",o)
stage_finish(r,o);check(stage_valid(r,o))
writeLines("changed",o);check(stage_valid(r,o))
writeLines("changed input",i);check(stage_valid(r,o))
stage_begin(r);check(!stage_valid(r,o))
# Scheduler quoting, dependencies, no actual bsub.
# Legacy adapter remains regression-tested independently of the new host backend.
c$execution$backend<-"legacy"
c$templates$job_script<-file.path(temp,"jobs/{stage}_{job}.sh")
c$templates$job_stdout<-file.path(temp,"logs/{stage}_{job}.out")
c$templates$job_stderr<-file.path(temp,"logs/{stage}_{job}.err")
txt<-capture.output(script<-submit_stage(c,root,"pilot",FALSE,"123"))
check(any(grepl("done(123)",txt,fixed=TRUE))&&any(grepl("DRY RUN ONLY",txt,fixed=TRUE)))
check(any(grepl("set -euo pipefail",readLines(script),fixed=TRUE)))
check(fails(submit_stage(c,root,"pilot",FALSE,"12;exit")))
check(fails(submit_stage(c,root,"pilot",TRUE,"none")))
check(fails(scale_gate(c,root,m)))
check(fails(completed_counts(c,root,m)))
# Human direct entry cannot count the same family twice.
human<-data.table(sample_id=c("one","two"),family_id="same")
check(fails(fit_promoter_pair(human)))
check(is.null(hc3_intercept(1:6,cbind(1,1:6))))
# Deterministic algebra fixtures exercise the human gates, not biological observations.
hd<-data.table(sample_id=paste0("truth",1:12),family_id=paste0("family",1:12),
  gene_id="TRUTH",gene_name="TRUTH",event_a="A",event_b="B",pair_type="structural_promoter_pair",
  maternal_count_a=11:22,paternal_count_a=20,maternal_count_b=20,paternal_count_b=12:23,
  ref_maternal_a=rep(c(TRUE,FALSE),6),ref_maternal_b=rep(c(TRUE,FALSE),6))
check(fit_promoter_pair(hd[1:9])$status=="insufficient_donors")
check(fit_promoter_pair(copy(hd)[,ref_maternal_a:=TRUE])$status=="insufficient_cis_orientation")
check(fit_promoter_pair(copy(hd)[,covar:=NA_real_],"covar")$status=="missing_covariates")
check(fit_promoter_pair(copy(hd)[,covar:=factor("one")],"covar")$status=="unidentifiable_covariates")
fit<-fit_promoter_pair(hd)
check(fit$cis_basis=="one_redundant_cis_term_removed"&&fit$residual_df==10)
check(fit$status=="eligible_conditional_inference"&&fit$bootstrap_valid_fraction>=.9)
rare<-copy(hd[1:10])
rare[,`:=`(ref_maternal_a=c(rep(TRUE,5),rep(FALSE,5)),
  ref_maternal_b=c(rep(TRUE,7),rep(FALSE,3)))]
rare_fit<-fit_promoter_pair(rare)
check(rare_fit$status=="unstable_bootstrap"&&is.na(rare_fit$p_value))
# Real samtools on a tiny coordinate-only FASTA/GTF fixture; no public cohort or quantification.
if(nzchar(Sys.which("samtools"))) {
  rcfg<-read_config_v2(file.path(root,"config/GUSTO.example.yaml"))
  keys<-c("genome","annotation","gtf","transcripts","exons","tx_map","event_map","events","pairs",
    "sizes","annotation_audit","deferred")
  for(k in keys)rcfg$paths[[k]]<-file.path(temp,paste0("reference_",k))
  rcfg$paths$genome_fai<-paste0(rcfg$paths$genome,".fai")
  rcfg$paths$transcripts_fai<-paste0(rcfg$paths$transcripts,".fai")
  rcfg$tools$samtools<-Sys.which("samtools")
  writeLines(c(">chrT",paste(rep("ACGT",20),collapse="")),rcfg$paths$genome)
  records<-list()
  addtx<-function(id,gene,strand,starts,ends) {
    attr<-paste0('gene_id "',gene,'"; gene_name "',gene,'"; transcript_id "',id,'";')
    records[[length(records)+1L]]<<-paste("chrT","truth","transcript",min(starts),max(ends),".",strand,".",attr,sep="\t")
    for(i in seq_along(starts))records[[length(records)+1L]]<<-paste("chrT","truth","exon",starts[i],
      ends[i],".",strand,".",paste0(attr,' exon_number "',i,'";'),sep="\t")
  }
  addtx("plus1","G","+",c(1,21),c(8,28));addtx("plus2","G","+",c(11,21),c(18,28))
  addtx("minus","H","-",c(31,51),c(38,58))
  writeLines(unlist(records),rcfg$paths$annotation)
  build_reference(rcfg)
  check(nrow(fread(rcfg$paths$tx_map))==3)
  check(all(Biostrings::width(Biostrings::readDNAStringSet(rcfg$paths$transcripts))==16))
  vv<-data.table(chrom="chrT",genomic_pos=58L,ref="C",alt="A",gt="0/1")
  mapped<-tx_variant_map(vv,fread(rcfg$paths$exons),rcfg,file.path(temp,"truth.vcf"),
    file.path(temp,"truth_map.tsv"),"TRUTH")
  check(mapped$tx_pos==1L&&mapped$tx_ref=="G"&&mapped$tx_alt=="T")
  check(nrow(fread(rcfg$paths$pairs))==1)
}
# Runtime script scan: no site/session paths or arbitrary environment_prefix.
scripts<-c(file.path(root,"poe.R"),list.files(file.path(root,"R"),"\\.R$",full.names=TRUE))
text<-unlist(lapply(scripts,readLines))
check(!any(grepl("/workspace|/mnt/|/rsrch[0-9]|environment_prefix|eval\\(parse",text)))
unlink(temp,recursive=TRUE)
cat("PASS:",n,"schema/acquisition/public-model/scheduler assertions\n")
