# Seadragon host/R bridge. R never launches a native bioinformatics executable.
is_seadragon <- function(c)identical(c$execution$backend,"seadragon")
runtime_record <- function(c) {
  x<-Sys.getenv("POE_NATIVE_RECORD")
  if(nzchar(x))x else c$paths$native_inventory
}
native_gate <- function(c,establish=FALSE) {
  for(p in c("data.table","yaml","curl","GenomicRanges","IRanges","S4Vectors","Rsamtools","Biostrings"))
    if(!requireNamespace(p,quietly=TRUE))stop("Missing R-container package: ",p)
  if(getRversion()<"4.3.1")stop("R >=4.3.1 required")
  f<-runtime_record(c)
  if(!file.exists(f))stop("Host inventory missing; run 00_setup/preflight.sh CONFIG --execute")
  d<-fread(f,colClasses="character",na.strings=NULL)
  require_columns(d,c("component","path","version"))
  assert_unique(d,"component")
  need<-c("ornaments","personalizer","samtools","bcftools","python","python_dependencies","site","r_image")
  if(!all(need %in% d$component)||anyNA(d$path)||any(!nzchar(d$path)))
    stop("Incomplete host inventory")
  val<-function(k,field)d[[field]][match(k,d$component)]
  if(!dir.exists(c$tools$r_library))stop("Declared R library is not visible inside the R container")
  if(!identical(normalizePath(val("site","path"),mustWork=TRUE),
                normalizePath(c$execution$site_file,mustWork=TRUE)))
    stop("Executed site file differs from the explicitly configured site file")
  if(!identical(val("r_image","path"),c$execution$image))
    stop("Executed R image differs from the configured image")
  invisible(d)
}
catalog_outputs <- function(c)unlist(c$paths[c("genome_fai","transcripts_fai","gtf","transcripts",
  "exons","tx_map","event_map","events","pairs","annotation_audit","deferred")])
catalog_gate <- function(c,root) {
  if(!stage_valid(c$paths$catalog_receipt,catalog_outputs(c)))
    stop("Missing/incomplete reference catalog")
}
marker_outputs <- function(c)unlist(c$paths[c("markers","tx_vcf","tx_marker_map")])
marker_gate <- function(c,root) {
  catalog_gate(c,root)
  if(!stage_valid(c$paths$marker_receipt,marker_outputs(c)))
    stop("Missing/incomplete parental marker/transcript VCF receipt")
}
child_outputs <- function(c,sample)vapply(c("child_query","personal_query","child_vcf","personal_vcf",
  "child_tx_vcf","child_tx_map"),function(k)cfg_path(c,k,sample=sample),character(1))
child_gate <- function(c,root,sample) {
  reference_gate(c,root)
  if(!stage_valid(cfg_path(c,"child_receipt",sample=sample),child_outputs(c,sample)))
    stop("Missing/incomplete typed child preparation")
}
child_index_gate <- function(c,root,sample) {
  child_gate(c,root,sample)
  outputs<-vapply(c("child_fasta","child_index"),function(k)cfg_path(c,k,sample=sample),character(1))
  if(!stage_valid(cfg_path(c,"human_receipt",sample=sample),outputs))
    stop("Missing/incomplete child personalized index")
}
plan_line <- function(tool,args,stdout=NULL,stderr=NULL) {
  cmd<-paste(vapply(c(tool,as.character(args)),shQuote,character(1)),collapse=" ")
  if(!is.null(stdout))cmd<-paste(cmd,">",shQuote(stdout))
  if(!is.null(stderr))cmd<-if(identical(stdout,stderr))paste(cmd,"2>&1") else paste(cmd,"2>",shQuote(stderr))
  cmd
}
write_plan <- function(lines,job,name="native") {
  p<-file.path(job,paste0(name,".sh"))
  writeLines(c("#!/bin/bash","set -eo pipefail","umask 077",lines),p)
  p
}
sample_manifest <- function(c,sample) {
  m<-load_manifests(c)
  if(length(sample)!=1L||!sample %in% m$libraries$sample_id)stop("Exact manifest library required")
  m
}
prepare_quant <- function(c,root,m,sample,mode,strand,limit,job) {
  l<-m$libraries[sample_id==sample]
  s<-m$sources[sample_id==sample]
  if(c$profile=="human_trio")child_index_gate(c,root,sample) else reference_gate(c,root)
  files<-character();n<-0;raw<-character()
  for(i in seq_len(nrow(s))) {
    a<-source_fastq(c,s[i],1);b<-if(l$layout=="PE")source_fastq(c,s[i],2) else NULL
    raw<-c(raw,a,b)
    if(!all(file.exists(c(a,b))))stop("Raw FASTQ missing")
    if(c$profile=="public_cross")for(k in seq_len(if(l$layout=="PE")2 else 1)) {
      if(!verify_file(c(a,b)[k],s[[paste0("bytes_",k)]][i],s[[paste0("md5_",k)]][i]))
        stop("Public source size/checksum mismatch")
    }
    if(is.finite(limit)) {
      if(n>=limit)next
      pa<-cfg_path(c,"prefix",sample=sample,run=s$run_id[i],mate=1,mode=mode,strand=strand)
      pb<-if(l$layout=="PE")cfg_path(c,"prefix",sample=sample,run=s$run_id[i],mate=2,mode=mode,strand=strand) else NULL
      ensure_parent(pa);n<-n+validate_fastq(a,b,limit-n,pa,pb);files<-c(files,pa,pb)
    } else {n<-n+validate_fastq(a,b);files<-c(files,a,b)}
  }
  genotype<-if(c$profile=="human_trio")human_linkage(c,m)[sample_id==sample]$child_id else c$cohort_id
  index<-if(c$profile=="human_trio")cfg_path(c,"child_index",sample=sample) else c$paths$index
  vcf<-if(c$profile=="human_trio")cfg_path(c,"child_tx_vcf",sample=sample) else c$paths$tx_vcf
  p<-function(k)cfg_path(c,k,sample=sample,mode=mode,strand=strand)
  outputs<-vapply(c("allele_counts","thetas","tpms","quant_log","resource_log"),p,character(1))
  stage_begin(p("quant_receipt"))
  for(o in outputs)ensure_parent(o)
  dir.create(p("quant_dir"),recursive=TRUE,showWarnings=FALSE)
  args<-quant_args(c,l,index,vcf,genotype,files,p("quant_dir"),strand)
  lines<-vapply(raw,function(f)plan_line(c$tools$gzip,c("-t",f)),character(1))
  lines<-c(lines,plan_line(c$tools$time,c("-v","-o",p("resource_log"),c$tools$ornaments,args),
    p("quant_log"),p("quant_log")))
  list(sample=sample,mode=mode,strand=strand,fragments=n,outputs=outputs,
    raw=raw,files=files,index=index,vcf=vcf,receipt=p("quant_receipt"),lines=lines)
}
finish_quant <- function(c,root,q) {
  x<-read_ornament_counts(q$outputs["allele_counts"],c$ornaments$semantics)
  if(!nrow(x))stop("Empty allelic output is failed callability")
  log<-readLines(q$outputs["quant_log"],warn=FALSE)
  line<-grep("\\[quant\\] processed .* reads, .* reads pseudoaligned",log,value=TRUE)
  if(length(line)!=1L)stop("Missing/ambiguous quantification diagnostics")
  processed<-as.numeric(gsub(",","",sub(".*processed ([0-9,]+) reads,.*","\\1",line)))
  mapped<-as.numeric(gsub(",","",sub(".*reads, ([0-9,]+) reads pseudoaligned.*","\\1",line)))
  if(!is.finite(processed)||processed!=q$fragments||!is.finite(mapped)||mapped<0||mapped>processed)
    stop("Quantifier processed/pseudoaligned count mismatch")
  rl<-readLines(q$outputs["resource_log"],warn=FALSE)
  number<-function(label) {
    line<-grep(label,rl,value=TRUE,fixed=TRUE)
    if(length(line)!=1L)stop("Missing native resource metric: ",label)
    as.numeric(sub(".*: *","",line))
  }
  rss<-number("Maximum resident set size")
  # Native wall time is measured independently of R/container startup.
  elapsed_line<-grep("Elapsed (wall clock)",rl,value=TRUE,fixed=TRUE)
  if(length(elapsed_line)!=1L)stop("Missing elapsed wall time")
  parts<-as.numeric(strsplit(sub("^.*\\): *","",elapsed_line),":",fixed=TRUE)[[1]])
  elapsed<-sum(rev(parts)*60^(seq_along(parts)-1L))
  if(!is.finite(elapsed)||elapsed<=0||!is.finite(rss)||rss<=0)stop("Invalid native resource metrics")
  l<-load_manifests(c)$libraries[sample_id==q$sample]
  metrics<-data.table(sample_id=q$sample,mode=q$mode,strand=q$strand,fragments=q$fragments,
    pseudoalignment_fraction=mapped/processed,elapsed_seconds=elapsed,peak_rss_kb=rss,
    work_bytes_endpoint_lower_bound=sum(file.info(list.files(c$paths$work,recursive=TRUE,full.names=TRUE))$size,na.rm=TRUE),
    index_bytes=file.info(q$index)$size,fragment_mean_nt=l$fragment_mean_nt,
    fragment_sd_nt=l$fragment_sd_nt,fragment_parameter_origin=l$fragment_parameter_origin,
    marker_proxy=isTRUE(c$markers$proxy),proxy_reason=if(isTRUE(c$markers$proxy))c$markers$proxy_reason else "")
  stage_finish(q$receipt,q$outputs)
  r<-readRDS(q$receipt);r$metrics<-metrics;saveRDS(r,q$receipt)
  metrics
}
library_pilot_receipt <- function(c,s)cfg_path(c,"library_pilot_receipt",sample=s)
library_pilot_report <- function(c,s)cfg_path(c,"library_pilot_report",sample=s)
collect_pilot <- function(c,root,m) {
  stage_begin(c$paths$pilot_receipt)
  reports<-lapply(m$libraries$sample_id,function(s) {
    p<-library_pilot_report(c,s)
    if(!stage_valid(library_pilot_receipt(c,s),p))stop("Incomplete pilot library: ",s)
    fread(p)
  })
  report<-rbindlist(reports)
  total<-sum(m$sources$read_count/ifelse(m$sources$layout=="PE",2,1))
  if(!is.finite(total)||total<=0)stop("Complete-cohort read counts required")
  throughput<-report[mode=="pilot",sum(fragments)/sum(elapsed_seconds)]
  report[,`:=`(projected_quant_seconds=total/throughput,configured_total_fragments=total,
    prefix_sampling="nonrandom_prefix_not_whole_file_checksum_identity",
    actual_focal_parent=if(c$profile=="public_cross")c$parents$focal else "typed_trio",
    actual_alternate_parent=if(c$profile=="public_cross")c$parents$alternate else "typed_trio")]
  if(any(report$projected_quant_seconds>c$pilot$max_projected_seconds)||
     any(report$work_bytes_endpoint_lower_bound>c$pilot$max_work_bytes))
    stop("Pilot exceeds configured compute/storage budget")
  write_tsv(report,c$paths$pilot_report)
  stage_finish(c$paths$pilot_receipt,c$paths$pilot_report)
  cat("Pilot report complete. Review it before explicitly submitting scale.\n")
}
runtime_action <- function(c,root,stage,phase,sample,job) {
  if(phase=="bootstrap") {
    if(!grepl("^[A-Za-z0-9_.-]+$",stage)||!grepl("^[A-Za-z0-9_.-]+$",sample))stop("Unsafe stage/library")
    if(!isTRUE(c$execution$posix_work_approved))stop("POSIX work/scratch attestation required")
    id<-paste0(stage,"_",sample,"_",Sys.getpid(),"_",format(Sys.time(),"%Y%m%dT%H%M%S"))
    job<-file.path(c$paths$scratch,"jobs",id);dir.create(job,recursive=TRUE,showWarnings=FALSE)
    key<-if(stage %in% c("reference","markers")||(stage=="index"&&c$profile=="public_cross"))"reference" else
      if(stage %in% c("child","index","pilot","quantify","assign"))paste0("library_",sample) else paste(stage,sample,sep="_")
    lock<-file.path(c$paths$receipts,"locks",key);ensure_parent(lock)
    cat(job,lock,sep="\t");return(invisible())
  }
  if(phase=="inventory-path") {cat(c$paths$native_inventory);return(invisible())}
  if(phase=="accept-inventory") {
    native_gate(c,establish=TRUE);ensure_parent(c$paths$native_inventory)
    if(!identical(runtime_record(c),c$paths$native_inventory))
      writeLines(readLines(runtime_record(c)),c$paths$native_inventory)
    cat("Host inventory recorded; native and R environments checked separately.\n")
    return(invisible())
  }
  native_gate(c);snapshot(c);m<-load_manifests(c)
  if(stage=="preflight") {preflight_v2(c,root);return(invisible())}
  if(stage=="reference") {
    if(phase=="prepare") {
      stage_begin(c$paths$catalog_receipt);stage_begin(c$paths$reference_receipt)
      lines<-character()
      if(!is.null(c$reference$genome_source)&&!file.exists(c$paths$genome)) {
        ensure_parent(c$paths$genome)
        lines<-plan_line("/bin/cp",c("--reflink=auto","--",c$reference$genome_source,c$paths$genome))
      }
      write_plan(c(lines,plan_line(c$tools$samtools,c("faidx",c$paths$genome))),job,"pre")
    } else if(phase=="middle") {
      build_reference(c,native=FALSE,public=FALSE)
      write_plan(plan_line(c$tools$samtools,c("faidx",c$paths$transcripts)),job,"post")
    } else if(phase=="finish") {
      stage_finish(c$paths$catalog_receipt,catalog_outputs(c))
      if(c$profile=="human_trio")
        stage_finish(c$paths$reference_receipt,reference_outputs(c))
    }
  } else if(stage=="markers") {
    if(c$profile!="public_cross")stop("Public markers stage cannot run on human trios")
    stage_begin(c$paths$marker_receipt)
    catalog_gate(c,root);v<-prepare_markers(c);v[,gt:="0/1"]
    tx_variant_map(v,fread(c$paths$exons),c,c$paths$tx_vcf,c$paths$tx_marker_map,c$cohort_id)
    stage_finish(c$paths$marker_receipt,marker_outputs(c))
  } else if(stage=="child") {
    if(c$profile!="human_trio")stop("Child stage requires human trio profile")
    reference_gate(c,root);sample_manifest(c,sample)
    s<-human_linkage(c,m)[sample_id==sample];if(nrow(s)!=1L)stop("Unique child required")
    if(phase=="prepare") {
      stage_begin(cfg_path(c,"child_receipt",sample=sample))
      lines<-character()
      for(mode in c("typed","personalization")) {
        ids<-if(mode=="typed")c(s$child_id,s$mother_id,s$father_id) else s$child_id
        vcf<-cfg_path(c,if(mode=="typed")"child_vcf" else "personal_vcf",sample=sample)
        query<-cfg_path(c,if(mode=="typed")"child_query" else "personal_query",sample=sample)
        ensure_parent(vcf);tmp<-paste0(vcf,".normalized")
        lines<-c(lines,
          plan_line(c$tools$bcftools,c("view","-s",paste(ids,collapse=","),"-m2","-M2","-v","snps",
            "-f","PASS,.","-Ov","-o",vcf,c$paths[[paste0(mode,"_vcf")]])),
          plan_line(c$tools$bcftools,c("norm","-f",c$paths$genome,"-c","e","-Ov","-o",tmp,vcf)),
          plan_line("/bin/mv",c("--",tmp,vcf)),
          plan_line(c$tools$bcftools,c("query","-s",paste(ids,collapse=","),
            "-f","%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n","-o",query,vcf)))
      }
      write_plan(lines,job,"pre")
    } else {
      d<-fread(cfg_path(c,"personal_query",sample=sample),header=FALSE,
        col.names=c("chrom","genomic_pos","ref","alt","gt"))
      d<-d[gt %in% c("0/1","1/0","1/1","0|1","1|0","1|1")]
      tx_variant_map(d,fread(c$paths$exons),c,cfg_path(c,"child_tx_vcf",sample=sample),
        cfg_path(c,"child_tx_map",sample=sample),s$child_id)
      stage_finish(cfg_path(c,"child_receipt",sample=sample),child_outputs(c,sample))
    }
  } else if(stage=="index") {
    if(c$profile=="public_cross") {
      marker_gate(c,root);vcf<-c$paths$tx_vcf;genotype<-c$cohort_id
      fasta<-c$paths$personal_fasta;index<-c$paths$index
    } else {
      child_gate(c,root,sample);vcf<-cfg_path(c,"child_tx_vcf",sample=sample)
      genotype<-human_linkage(c,m)[sample_id==sample]$child_id
      fasta<-cfg_path(c,"child_fasta",sample=sample);index<-cfg_path(c,"child_index",sample=sample)
    }
    if(phase=="prepare") {
      stage_begin(if(c$profile=="public_cross")c$paths$reference_receipt else
        cfg_path(c,"human_receipt",sample=sample))
      ensure_parent(fasta);ensure_parent(index)
      write_plan(c(plan_line(c$tools$personalizer,c("-f",c$paths$transcripts,"-v",vcf,"-t","ornament",
        "-o",fasta,"-s",genotype)),plan_line(c$tools$ornaments,c("index","-i",index,fasta))),job,"pre")
    } else {
      if(!all(file.exists(c(fasta,index)))||any(file.info(c(fasta,index))$size<=0))stop("Missing/empty personalized index")
      if(c$profile=="public_cross")
        stage_finish(c$paths$reference_receipt,reference_outputs(c)) else
        stage_finish(cfg_path(c,"human_receipt",sample=sample),c(fasta,index))
    }
  } else if(stage %in% c("pilot","quantify")) {
    sample_manifest(c,sample);reference_gate(c,root)
    if(stage=="quantify")scale_gate(c,root,m)
    state<-file.path(job,"quant_state.rds")
    if(phase=="prepare") {
      if(stage=="pilot")stage_begin(library_pilot_receipt(c,sample))
      l<-m$libraries[sample_id==sample]
      specs<-if(stage=="pilot")c(lapply(c("fr","rf","unstranded"),function(s)
        list(mode="strand",strand=s,limit=c$pilot$strand_fragments)),
        list(list(mode="pilot",strand=l$strand_mode,limit=c$pilot$fragments))) else
        list(list(mode="complete",strand=l$strand_mode,limit=Inf))
      qs<-lapply(specs,function(s)prepare_quant(c,root,m,sample,s$mode,s$strand,s$limit,job))
      saveRDS(qs,state);write_plan(unlist(lapply(qs,`[[`,"lines")),job,"pre")
    } else {
      qs<-readRDS(state);report<-rbindlist(lapply(qs,function(q)finish_quant(c,root,q)))
      if(stage=="pilot") {
        l<-m$libraries[sample_id==sample];diag<-report[mode=="strand"]
        chosen<-diag[strand==l$strand_mode]$pseudoalignment_fraction
        if(!all(l$strand_verified %in% TRUE)||length(chosen)!=1L||
          chosen<c$pilot$min_pseudoalignment_fraction||
          (l$strand_mode!="unstranded"&&chosen<max(diag[strand!="unstranded"]$pseudoalignment_fraction)-c$pilot$strand_tolerance)||
          report[mode=="pilot",fragments]<c$pilot$fragments)stop("Pilot strand/protocol/chunk gate failed")
        assign_sample(c,m,sample,"pilot")
        write_tsv(report,library_pilot_report(c,sample))
        stage_finish(library_pilot_receipt(c,sample),library_pilot_report(c,sample))
      } else write_tsv(report,cfg_path(c,"library_metrics",sample=sample))
    }
  } else if(stage=="collect-pilot") {
    reference_gate(c,root);collect_pilot(c,root,m)
  } else if(stage=="assign") {
    scale_gate(c,root,m);sample_manifest(c,sample)
    if(c$profile=="human_trio")child_index_gate(c,root,sample) else reference_gate(c,root)
    l<-m$libraries[sample_id==sample]
    qr<-cfg_path(c,"quant_receipt",sample=sample,mode="complete",strand=l$strand_mode)
    if(!file.exists(qr))stop("Quantification receipt missing")
    quant_outputs<-vapply(c("allele_counts","thetas","tpms","quant_log","resource_log"),
      function(k)cfg_path(c,k,sample=sample,mode="complete",strand=l$strand_mode),character(1))
    if(!stage_valid(qr,quant_outputs))stop("Quantification receipt or outputs incomplete")
    ar<-cfg_path(c,"assign_receipt",sample=sample,mode="complete",strand=l$strand_mode)
    stage_begin(ar)
    assign_sample(c,m,sample)
    stage_finish(ar,assignment_outputs(c,sample,"complete",l))
  } else if(stage=="collect-quantification") {
    scale_gate(c,root,m);completed_counts(c,root,m)
    write_tsv(rbindlist(lapply(m$libraries$sample_id,function(s) {
      l<-m$libraries[sample_id==s]
      readRDS(cfg_path(c,"quant_receipt",sample=s,mode="complete",strand=l$strand_mode))$metrics
    })),c$paths$metrics)
  } else if(stage=="models")worker_models(c,root)
  else if(stage=="retrieve") {
    if(phase=="prepare") {
      cc<-c;cc$paths$resource_receipt<-file.path(job,"resource_receipt.tsv")
      retrieve_inputs(cc,native=FALSE)
      d<-fread(cc$paths$resource_receipt)
      lines<-vapply(d$path[endsWith(d$path,".gz")],function(f)plan_line(c$tools$gzip,c("-t",f)),character(1))
      write_plan(lines,job,"pre")
    } else write_tsv(fread(file.path(job,"resource_receipt.tsv")),c$paths$resource_receipt)
  } else stop("Unknown Seadragon stage: ",stage)
  invisible(TRUE)
}
