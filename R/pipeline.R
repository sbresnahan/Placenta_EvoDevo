ensure_parent <- function(p)dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
source_fastq <- function(c,s,mate) {
  col<-paste0("local_",mate)
  if(c$profile=="human_trio") {
    if(!col %in% names(s)||is.na(s[[col]])||!nzchar(s[[col]]))stop("Private local FASTQ path required")
    return(resolve_value(s[[col]],dirname(attr(c,"file")),unlist(c$environment_allowlist)))
  }
  cfg_path(c,"fastq",run=s$run_id,mate=mate)
}
tool_gate <- function(c) {
  if(identical(c$execution$backend,"seadragon"))return(native_gate(c))
  if(!dir.exists(c$tools$r_library))stop("Configured R library directory is not visible")
  for(p in c("data.table","yaml","curl","GenomicRanges","IRanges","S4Vectors","Rsamtools","Biostrings"))
    if(!requireNamespace(p,quietly=TRUE))stop("Missing configured R dependency: ",p)
  for(k in c("Rscript","bash","samtools","bcftools","gzip","ornaments","python","personalizer","time"))
    if(!file.exists(c$tools[[k]]))stop("Configured tool missing: ",k)
  if(getRversion()<"4.3.1")stop("R >=4.3.1 required")
  version_output<-function(tool) {
    x<-system2(tool,"--version",stdout=TRUE,stderr=TRUE)
    if(!is.null(attr(x,"status")))stop("Tool version check failed: ",basename(tool))
    paste(x,collapse="; ")
  }
  versions<-data.table(component=c("R","samtools","bcftools"),
    value=c(as.character(getRversion()),
      version_output(c$tools$samtools),version_output(c$tools$bcftools)))
  write_tsv(versions,c$paths$tool_versions)
}
snapshot <- function(c) {
  for(k in c("work","scratch","logs","receipts"))dir.create(c$paths[[k]],recursive=TRUE,showWarnings=FALSE)
  for(k in c("work","scratch","logs","receipts"))if(file.access(c$paths[[k]],2)!=0)stop("Unwritable root: ",k)
  ensure_parent(c$paths$snapshot)
  tmp<-tempfile(".snapshot-",tmpdir=dirname(c$paths$snapshot));yaml::write_yaml(c,tmp)
  if(!file.rename(tmp,c$paths$snapshot))stop("Atomic configuration snapshot failed")
}
human_linkage <- function(c,m) {
  s<-copy(m$libraries)
  require_columns(s,c("child_id"))
  t<-fread(c$paths$trios)
  require_columns(t,c("child_id","mother_id","father_id","family_id"))
  assert_unique(s,"child_id");assert_unique(t,"child_id")
  if(anyNA(t)||any(!s$child_id %in% t$child_id)||
    any(t$child_id==t$mother_id|t$child_id==t$father_id|t$mother_id==t$father_id))
    stop("Missing or conflicting sample/trio linkage")
  merge(s,t,by="child_id",all.x=TRUE)
}
preflight_v2 <- function(c,root,execution=FALSE) {
  tool_gate(c)
  cat("Required tools and R packages are accessible.\n")
  invisible(TRUE)
}
validate_inputs_v2 <- function(c,root,execution=FALSE) {
  m<-load_manifests(c)
  r<-fread(c$paths$resources)
  require_columns(r,c("resource_id","url","destination","partial","response","bytes","md5","required","build"))
  assert_unique(r,"resource_id")
  if(c$profile=="human_trio")human_linkage(c,m)
  if(execution) {tool_gate(c);snapshot(c)}
  cat("Configuration/manifests valid.\n")
  invisible(m)
}
make_index <- function(c,vcf,sample,fasta,index) {
  ensure_parent(fasta);ensure_parent(index)
  run_checked(c$tools$python,c(c$tools$personalizer,"-f",c$paths$transcripts,
    "-v",vcf,"-t","ornament","-o",fasta,"-s",sample))
  run_checked(c$tools$ornaments,c("index","-i",index,fasta))
}
human_prepare <- function(c,m,sample,root) {
  s<-human_linkage(c,m)[sample_id==sample]
  if(nrow(s)!=1)stop("Unique human sample required")
  outputs<-vapply(c("child_query","personal_query","child_vcf","personal_vcf",
    "child_tx_vcf","child_tx_map","child_fasta","child_index"),
    function(k)cfg_path(c,k,sample=sample),character(1))
  receipt<-cfg_path(c,"human_receipt",sample=sample)
  stage_begin(receipt)
  for(mode in c("typed","personalization")) {
    ids<-if(mode=="typed")c(s$child_id,s$mother_id,s$father_id) else s$child_id
    vcf<-cfg_path(c,if(mode=="typed")"child_vcf" else "personal_vcf",sample=sample)
    query<-cfg_path(c,if(mode=="typed")"child_query" else "personal_query",sample=sample)
    ensure_parent(vcf)
    # norm performs REF agreement; view keeps biallelic SNPs only.
    # Two explicit calls avoid shell pipelines and arbitrary shell interpolation.
    run_checked(c$tools$bcftools,c("view","-s",paste(ids,collapse=","),"-m2","-M2","-v","snps",
      "-f","PASS,.","-Ov","-o",vcf,c$paths[[paste0(mode,"_vcf")]]))
    run_checked(c$tools$bcftools,c("norm","-f",c$paths$genome,"-c","e","-Ov",vcf),
      stdout=query)
    if(!file.rename(query,vcf))stop("Normalized VCF staging failed")
    run_checked(c$tools$bcftools,c("query","-s",paste(ids,collapse=","),
      "-f","%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n","-o",query,vcf))
  }
  d<-fread(cfg_path(c,"personal_query",sample=sample),header=FALSE,
    col.names=c("chrom","genomic_pos","ref","alt","gt"))
  d<-d[gt %in% c("0/1","1/0","1/1","0|1","1|0","1|1")]
  tx_variant_map(d,fread(c$paths$exons),c,cfg_path(c,"child_tx_vcf",sample=sample),
    cfg_path(c,"child_tx_map",sample=sample),s$child_id)
  make_index(c,cfg_path(c,"child_tx_vcf",sample=sample),s$child_id,
    cfg_path(c,"child_fasta",sample=sample),cfg_path(c,"child_index",sample=sample))
  ensure_parent(receipt);stage_finish(receipt,outputs)
  s
}
quant_args <- function(c,l,index,vcf,genotype,files,out,strand=l$strand_mode) {
  a<-c("quant","-i",index,"-o",out,"--vcf",vcf,"--sample",genotype,
    "--threads",as.character(c$ornaments$threads))
  if(strand=="fr")a<-c(a,"--fr-stranded")
  if(strand=="rf")a<-c(a,"--rf-stranded")
  if(l$layout=="SE")a<-c(a,"--single","-l",l$fragment_mean_nt,"-s",l$fragment_sd_nt)
  c(a,files)
}
quant_sample <- function(c,root,m,sample,mode="complete",strand=NULL,limit=Inf) {
  l<-m$libraries[sample_id==sample];if(nrow(l)!=1)stop("Unknown sample")
  if(is.null(strand))strand<-l$strand_mode
  s<-m$sources[sample_id==sample]
  files<-character();n<-0
  for(i in seq_len(nrow(s))) {
    a<-source_fastq(c,s[i],1);b<-if(l$layout=="PE")source_fastq(c,s[i],2) else NULL
    if(!file.exists(a)||(!is.null(b)&&!file.exists(b)))stop("Raw FASTQ missing")
    if(c$profile=="public_cross")for(mate in seq_len(if(l$layout=="PE")2 else 1)) {
      f<-if(mate==1)a else b
      if(!verify_file(f,s[[paste0("bytes_",mate)]][i],s[[paste0("md5_",mate)]][i]))
        stop("Public raw source size/checksum mismatch")
    }
    if(is.finite(limit)) {
      if(n>=limit)break
      pa<-cfg_path(c,"prefix",sample=sample,run=s$run_id[i],mate=1,mode=mode,strand=strand)
      pb<-if(l$layout=="PE")cfg_path(c,"prefix",sample=sample,run=s$run_id[i],mate=2,mode=mode,strand=strand) else NULL
      ensure_parent(pa)
      k<-validate_fastq(a,b,limit-n,pa,pb);n<-n+k
      files<-c(files,pa,pb)
    } else {
      # Completeness, CRC and mate identity must be checked, including private FASTQs.
      run_checked(c$tools$gzip,c("-t",a));if(!is.null(b))run_checked(c$tools$gzip,c("-t",b))
      n<-n+validate_fastq(a,b);files<-c(files,a,b)
    }
  }
  if(c$profile=="human_trio") {
    hs<-human_prepare(c,m,sample,root);genotype<-hs$child_id
    index<-cfg_path(c,"child_index",sample=sample);vcf<-cfg_path(c,"child_tx_vcf",sample=sample)
  } else {genotype<-c$cohort_id;index<-c$paths$index;vcf<-c$paths$tx_vcf}
  args<-list(sample=sample,mode=mode,strand=strand)
  p<-function(k)do.call(cfg_path,c(list(c=c,key=k),args))
  outputs<-vapply(c("allele_counts","thetas","tpms","quant_log","resource_log"),p,character(1))
  rec<-p("quant_receipt")
  stage_begin(rec)
  for(o in outputs)ensure_parent(o)
  dir.create(p("quant_dir"),recursive=TRUE,showWarnings=FALSE)
  t<-proc.time()[["elapsed"]]
  qa<-quant_args(c,l,index,vcf,genotype,files,p("quant_dir"),strand)
  run_checked(c$tools$time,c("-v","-o",p("resource_log"),c$tools$ornaments,qa),
    stdout=p("quant_log"),stderr=p("quant_log"))
  elapsed<-proc.time()[["elapsed"]]-t
  x<-read_ornament_counts(p("allele_counts"),c$ornaments$semantics)
  if(!nrow(x))stop("Empty allelic output is failed callability, not successful completion")
  log<-readLines(p("quant_log"),warn=FALSE)
  q<-grep("\\[quant\\] processed .* reads, .* reads pseudoaligned",log,value=TRUE)
  if(length(q)!=1)stop("Cannot verify Ornaments processed/pseudoaligned diagnostics")
  processed<-as.numeric(gsub(",","",sub(".*processed ([0-9,]+) reads,.*","\\1",q)))
  mapped<-as.numeric(gsub(",","",sub(".*reads, ([0-9,]+) reads pseudoaligned.*","\\1",q)))
  if(processed!=n)stop("Quantifier did not process all configured fragments")
  rssline<-grep("Maximum resident set size",readLines(p("resource_log")),value=TRUE)
  rss<-as.numeric(sub(".*: *","",rssline))
  disk<-sum(file.info(list.files(c$paths$work,recursive=TRUE,full.names=TRUE))$size,na.rm=TRUE)
  metrics<-data.table(sample_id=sample,mode=mode,strand=strand,fragments=n,
    pseudoalignment_fraction=mapped/processed,elapsed_seconds=elapsed,peak_rss_kb=rss,
    work_bytes_endpoint_lower_bound=disk,index_bytes=file.info(index)$size,
    fragment_mean_nt=if(l$layout=="SE")l$fragment_mean_nt else NA_real_,
    fragment_sd_nt=if(l$layout=="SE")l$fragment_sd_nt else NA_real_,
    fragment_parameter_origin=if(l$layout=="SE")l$fragment_parameter_origin else "paired_end_quantifier",
    marker_proxy=isTRUE(c$markers$proxy),proxy_reason=if(isTRUE(c$markers$proxy))c$markers$proxy_reason else "")
  ensure_parent(rec);stage_finish(rec,outputs)
  r<-readRDS(rec);r$metrics<-metrics;saveRDS(r,rec);metrics
}
assign_sample <- function(c,m,sample,mode="complete") {
  l<-m$libraries[sample_id==sample];strand<-l$strand_mode
  p<-function(k)cfg_path(c,k,sample=sample,mode=mode,strand=strand)
  x<-read_ornament_counts(p("allele_counts"),c$ornaments$semantics)
  ex<-fread(c$paths$exons);tx<-fread(c$paths$tx_map)
  if(any(!x$transcript_id %in% tx$transcript_id))stop("Counts/reference transcript mismatch")
  x<-map_tx_loci(x,ex);x[,locus_id:=paste(chrom,genomic_pos,ref,alt,sep=":")]
  map<-fread(c$paths$event_map)[window==50]
  ev<-fread(c$paths$events)[window==50]
  intervals<-rbindlist(lapply(seq_len(nrow(ev)),function(i) {
    if(is.na(ev$distinguishing_intervals[i])||!nzchar(ev$distinguishing_intervals[i]))return(NULL)
    a<-tstrsplit(strsplit(ev$distinguishing_intervals[i],";",fixed=TRUE)[[1]],"-",fixed=TRUE)
    data.table(event_id=ev$event_id[i],chrom=ev$chrom[i],start=as.integer(a[[1]]),end=as.integer(a[[2]]))
  }))
  if(!ncol(intervals))intervals<-data.table(event_id=character(),chrom=character(),start=integer(),end=integer())
  x<-mark_direct_coordinates(x,intervals,map)
  if(c$profile=="public_cross") {
    unit<-m$units[unit_id==l$unit_id];x<-orient_public(x,fread(c$paths$markers),unit,c);family<-unit$unit_id
  } else {
    hs<-human_linkage(c,m)[sample_id==sample];family<-hs$family_id
    anchors<-fread(cfg_path(c,"child_query",sample=sample),header=FALSE,
      col.names=c("chrom","genomic_pos","ref","alt","child_gt","mother_gt","father_gt"))
    x<-orient_counts(x,anchors)
  }
  # Cross-gene overlapping exons cannot provide a unique gene anchor.
  x[,gene_id:=tx$gene_id[match(transcript_id,tx$transcript_id)]]
  valid<-which(!is.na(x$genomic_pos))
  h<-GenomicRanges::findOverlaps(GenomicRanges::GRanges(x$chrom[valid],
    IRanges::IRanges(x$genomic_pos[valid],width=1)),GenomicRanges::GRanges(ex$chrom,
    IRanges::IRanges(ex$start,ex$end)),ignore.strand=TRUE)
  bad<-unique(valid[S4Vectors::queryHits(h)[
    x$gene_id[valid[S4Vectors::queryHits(h)]]!=ex$gene_id[S4Vectors::subjectHits(h)]]])
  x[bad,`:=`(maternal_count=NA_real_,paternal_count=NA_real_,ref_maternal=NA,
    parent_status="other_gene_overlap")]
  gene<-x[,.(maternal_count=if(anyNA(maternal_count))NA_real_ else sum(maternal_count),
    paternal_count=if(anyNA(paternal_count))NA_real_ else sum(paternal_count)),
    by=.(gene_id,locus_id)]
  x[,gene_id:=NULL]
  event<-event_locus_counts(x[transcript_id %in% map$transcript_id],map,sample,family)
  for(z in c("gene","event")) {
    d<-get(z);d[,`:=`(unit_id=l$unit_id,marker_proxy=isTRUE(c$markers$proxy),
      evidence_tier="model_based_expected_counts_not_joint_molecules")]
    write_tsv(d,p(paste0(z,"_counts")))
  }
  write_tsv(x[,.N,by=.(mapping_status,parent_status)],p("assignment_qc"))
}
