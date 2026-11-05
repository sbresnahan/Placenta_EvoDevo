#!/usr/bin/env Rscript
# Run on Seadragon via bin/Rscript_sif tests/test_additional.R
f<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
root<-normalizePath(file.path(dirname(f),".."))
source(file.path(root,"R/library_paths.R"))
for(k in c("common","workflow","parental","phase_origin","models","public","additional"))
  source(file.path(root,"R",paste0(k,".R")))
library(data.table)

# Chromosome vote from homozygous parental transmissions; detect a phase flip.
trio<-data.table(chrom=rep("1",50),genomic_pos=seq(1000000,50000000,by=1000000),
  ref="A",alt="G",child_gt="0|1",mother_gt="0/0",father_gt="1/1")
s<-list(min_informative=20L,min_agreement=0.95,window_bp=5000000L,
  min_window_votes=5L,min_covered_windows=3L)
p<-infer_chromosome_origin(trio,s)
stopifnot(nrow(p)==1L,p$status=="callable",p$maternal_haplotype==1L)
one<-data.table(chrom="1",genomic_pos=1000000L,ref="A",alt="G",
  transcript_id="tx",ref_count=7,alt_count=3)
oriented<-orient_counts_chromosome(one,trio[1],p)
stopifnot(oriented$parent_status=="phased_chromosome_trio",
  oriented$maternal_count==7,oriented$paternal_count==3)
trio[genomic_pos>25000000,child_gt:="1|0"]
switched<-infer_chromosome_origin(trio,s)
stopifnot(switched$status!="callable")

# TSS and donor are mutually exclusive, conditional structural comparisons.
ev<-data.table(window=50L,gene_id="g",gene_name="g",chrom="1",strand="+",
  event_id=c("a","b","c"),tss_min=c(100L,200L,100L),
  first_donor=c(300L,300L,400L),distinguishing_bp=30L)
q<-additional_event_pairs(ev,list())
stopifnot(nrow(q)==2L,
  setequal(q$effect_type,c("TSS_conditional_on_first_donor","first_donor_conditional_on_TSS")))

units<-data.table(unit_id=paste0("u",1:4),cluster_id=paste0("cluster",1:4),
  maternal_parent=c("B73","B73","Ki11","Ki11"))
counts<-data.table(gene_id="g",event_id=rep(c("a","b"),each=4),
  unit_id=rep(units$unit_id,2),locus_id="1:100:A:G",
  maternal_count=c(14,12,10,11,9,10,12,11),
  paternal_count=c(9,10,12,11,14,13,9,10),
  distinguishing_coordinate=FALSE)
cfg<-list(parents=list(focal="B73"),dosage=list(maternal=2,paternal=1))
settings<-list(min_depth=10L,min_per_direction=2L,max_pairs_per_gene=10L)
result<-additional_public_pairs(counts,ev,units,cfg,settings)
tss<-result[effect_type=="TSS_conditional_on_first_donor"]
stopifnot(nrow(result)==2L,nrow(tss)==1L,tss$n_units==4L,
  tss$n_focal_mother==2L,tss$n_alternate_mother==2L)

cat("PASS: chromosome phase and distinct biological event hypotheses\n")
