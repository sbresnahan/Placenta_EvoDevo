pin <- "2c4b1a2c18616f9de1cb307749688daaf3313d17"
scalar_text <- function(x, label) {
  if (!is.character(x)||length(x)!=1L||is.na(x)||!nzchar(x)) stop("Required scalar text: ",label)
  x
}
positive <- function(x,label) {
  if (!is.numeric(x)||length(x)!=1L||!is.finite(x)||x<=0) stop("Required positive number: ",label)
  invisible(x)
}
resolve_value <- function(x,base,allowed=character()) {
  scalar_text(x,"path")
  hits<-regmatches(x,gregexpr("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}",x))[[1]]
  for(h in hits) {
    k<-substring(h,3,nchar(h)-1L)
    if(!k %in% allowed||!nzchar(Sys.getenv(k))) stop("Unresolved/unapproved variable: ",k)
    x<-gsub(h,Sys.getenv(k),x,fixed=TRUE)
  }
  if(grepl("\\$",x))stop("Unsupported environment expression")
  if(!grepl("^/",x))x<-file.path(base,x)
  x
}
read_config_v2 <- function(path) {
  path<-normalizePath(path,mustWork=TRUE); c<-yaml::read_yaml(path,eval.expr=FALSE)
  if(!identical(c$schema_version,2L))stop("Schema 2 required; schema 1 is retired")
  if(!c$profile %in% c("public_cross","human_trio")||!c$species %in% c("maize","mouse","human"))
    stop("Unknown profile/species")
  if((c$species=="human")!=(c$profile=="human_trio"))stop("Mixed public/trio profile")
  scalar_text(c$cohort_id,"cohort_id")
  if(!grepl("^[A-Za-z0-9_.-]+$",c$cohort_id))stop("Unsafe cohort ID")
  if(!identical(c$ornaments$commit,pin)||!identical(c$ornaments$semantics,"ref_alt"))
    stop("Only audited REF/ALT Ornaments is supported")
  for(k in c("binary_sha256","personalizer_sha256")) {
    scalar_text(c$ornaments[[k]],paste0("ornaments.",k))
    if(!grepl("^[a-f0-9]{64}$",c$ornaments[[k]]))
      stop("Required installed-file SHA256 (64 lowercase hex characters): ornaments.",k)
  }
  required<-list(
    paths=c("work","scratch","logs","receipts","snapshot","sources","libraries","units","resources",
      "genome","annotation","gtf","transcripts","exons","tx_map","event_map","events","pairs",
      "marker_sources","markers","tx_vcf","tx_marker_map","personal_fasta","index","pilot_report",
      "pilot_receipt","scale_receipt","model_gene","model_promoter","model_pairs","model_receipt",
      "status","candidate","release","trios","covariates","typed_vcf","personalization_vcf","sizes",
      "annotation_audit","resource_receipt","reference_receipt","deferred","metrics","tool_versions",
      "genome_fai","transcripts_fai"),
    templates=c("fastq","partial","response","prefix","child_query","personal_query","child_vcf",
      "personal_vcf","child_tx_vcf","child_tx_map","child_fasta","child_index","quant_dir",
      "allele_counts","thetas","tpms","event_counts","gene_counts","assignment_qc","quant_receipt",
      "assign_receipt","human_receipt","quant_log","resource_log","job_script","job_stdout","job_stderr"),
    tools=c("Rscript","bash","bsub","git","samtools","bcftools","gzip","curl","ornaments",
      "python","personalizer","ornaments_source","singularity","r_library","time"))
  for(section in names(required)) {
    for(k in required[[section]])if(is.null(c[[section]][[k]]))stop("Missing ",section,".",k)
    c[[section]]<-lapply(c[[section]],resolve_value,base=dirname(path),allowed=unlist(c$environment_allowlist))
  }
  if(!is.logical(c$execution$container_enabled)||length(c$execution$container_enabled)!=1L)
    stop("container_enabled must be boolean")
  if(c$paths$genome_fai!=paste0(c$paths$genome,".fai")||
     c$paths$transcripts_fai!=paste0(c$paths$transcripts,".fai"))
    stop("Declare samtools/Rsamtools adjacent .fai paths explicitly")
  c$execution$image<-resolve_value(c$execution$image,dirname(path),unlist(c$environment_allowlist))
  if(!is.null(c$reference$genome_source))
    c$reference$genome_source<-resolve_value(c$reference$genome_source,dirname(path),unlist(c$environment_allowlist))
  if(identical(c$execution$backend,"seadragon")) {
    for(k in c("native_inventory","catalog_receipt","marker_receipt","job_registry"))
      if(is.null(c$paths[[k]]))stop("Missing Seadragon path: ",k)
    for(k in c("child_receipt","library_pilot_report","library_pilot_receipt","library_metrics"))
      if(is.null(c$templates[[k]]))stop("Missing Seadragon template: ",k)
    c$execution$site_file<-resolve_value(c$execution$site_file,dirname(path),unlist(c$environment_allowlist))
  }
  c$execution$binds<-lapply(c$execution$binds,function(b) {
    if(!is.list(b)||!b$mode %in% c("ro","rw"))stop("Structured bind required")
    b$host<-resolve_value(b$host,dirname(path),unlist(c$environment_allowlist))
    b$container<-resolve_value(b$container,dirname(path),unlist(c$environment_allowlist));b
  })
  for(s in c("pilot","scale","reference","models"))for(k in c("memory_value","cpus"))
    positive(c$lsf$stages[[s]][[k]],paste(s,k))
  for(k in c("fragments","strand_fragments","max_projected_seconds","max_work_bytes"))
    positive(c$pilot[[k]],paste("pilot",k))
  positive(c$models$min_depth,"models.min_depth")
  positive(c$ornaments$threads,"ornaments.threads")
  positive(c$fragments$fragment_mean_nt,"fragments.fragment_mean_nt")
  positive(c$fragments$fragment_sd_nt,"fragments.fragment_sd_nt")
  if(!c$fragments$fragment_parameter_origin %in% c("assumed_user_defined","measured"))
    stop("Explicit fragment-parameter origin required")
  if(c$ornaments$threads!=floor(c$ornaments$threads))stop("Integer thread count required")
  for(k in c("attempts","timeout_seconds","connect_timeout_seconds"))positive(c$retrieval[[k]],paste("retrieval",k))
  if(c$retrieval$attempts!=floor(c$retrieval$attempts))stop("Integer retry count required")
  if(c$models$min_depth<20)stop("Depth floor cannot be below 20")
  if(c$profile=="human_trio" && (c$models$min_families<10||c$models$min_orientation<3||
      c$models$bootstrap<2000||c$models$seed!=1147))stop("Human gates cannot be relaxed")
  for(x in c$templates) {
    tok<-regmatches(x,gregexpr("\\{[^{}]+\\}",x))[[1]]
    if(any(!tok %in% paste0("{",c("sample","run","mate","mode","stage","job","strand"),"}")))
      stop("Unknown output-template token")
  }
  if(c$profile=="public_cross") {
    scalar_text(c$parents$focal,"focal");scalar_text(c$parents$alternate,"alternate")
    if(c$parents$focal==c$parents$alternate)stop("Distinct parents required")
    positive(c$dosage$maternal,"maternal dosage");positive(c$dosage$paternal,"paternal dosage")
    if(!isTRUE(c$dosage$reviewed))stop("Tissue/ploidy-specific dosage review required")
    if(!c$markers$mode %in% c("panzea","mgp","prepared"))stop("Unknown marker adapter")
    positive(c$markers$chunk_rows,"markers.chunk_rows")
    if(c$markers$chunk_rows!=floor(c$markers$chunk_rows))stop("Integer VCF chunk size required")
    if(!is.finite(c$markers$min_ref_fraction)||c$markers$min_ref_fraction<.99||
       c$markers$min_ref_fraction>1)stop("Reference agreement floor is 99 percent")
    if(c$markers$mode=="panzea" &&
       (c$markers$min_depth<5||c$markers$min_allele_fraction<.9||
        c$markers$min_allele_fraction>1))stop("Panzea depth/allele-fraction gates cannot be relaxed")
    if(isTRUE(c$markers$proxy)&&(!isTRUE(c$markers$proxy_approved)||!nzchar(c$markers$proxy_reason)))
      stop("Explicit proxy approval/reason required")
    if(isTRUE(c$markers$proxy)&&c$parents$focal=="C57BL_6J_Babr"&&
       c$markers$mode=="mgp"&&c$markers$focal_sample!="REFERENCE_C57BL_6J")
      stop("The approved Babr MGP proxy is B6J reference, not B6N/another strain")
  }
  attr(c,"file")<-path;c
}
cfg_path <- function(c,key,...) {
  x<-c$templates[[key]]
  if(is.null(x))return(c$paths[[key]])
  vals<-list(...)
  for(k in names(vals)) {
    v<-as.character(vals[[k]])
    if(length(v)!=1L||is.na(v)||!grepl("^[A-Za-z0-9_.-]+$",v))stop("Unsafe template value")
    x<-gsub(paste0("{",k,"}"),v,x,fixed=TRUE)
  }
  if(grepl("[{}]",x))stop("Unresolved template: ",key)
  x
}
resolve_manifest_paths <- function(d,columns,c) {
  for(col in intersect(columns,names(d))) {
    ii<-which(!is.na(d[[col]]) & nzchar(d[[col]]))
    set(d,ii,col,vapply(d[[col]][ii],resolve_value,character(1),
       base=dirname(attr(c,"file")),allowed=unlist(c$environment_allowlist)))
  };d
}
load_manifests <- function(c,require_files=FALSE) {
  s<-fread(c$paths$sources,na.strings=c("","NA"));l<-fread(c$paths$libraries)
  for(k in intersect(c("bytes_1","bytes_2","read_count","base_count"),names(s)))
    set(s,j=k,value=as.numeric(s[[k]]))
  u<-fread(c$paths$units,na.strings=c("","NA"))
  require_columns(s,c("run_id","sample_id","layout","url_1","url_2","md5_1","md5_2","bytes_1","bytes_2"))
  require_columns(l,c("sample_id","unit_id","layout","strand_mode","strand_verified"))
  require_columns(u,c("unit_id","cluster_id","independence_verified","maternal_parent",
    "paternal_parent","stage","tissue"))
  assert_unique(s,"run_id");assert_unique(l,"sample_id");assert_unique(u,"unit_id")
  for(v in list(s$run_id,l$sample_id,u$unit_id))
    if(any(!grepl("^[A-Za-z0-9_.-]+$",v)))stop("Unsafe manifest identifier")
  assert_unique(l,"unit_id","One primary library per unit")
  if(!nrow(s)||!nrow(l)||!nrow(u)||any(!s$sample_id %in% l$sample_id)||
      any(!l$sample_id %in% s$sample_id)||!setequal(l$unit_id,u$unit_id))stop("Incomplete manifest linkage")
  if(any(!l$layout %in% c("SE","PE"))||any(!l$strand_mode %in% c("fr","rf","unstranded")))
    stop("Explicit layout/strand required")
  if(any(s$layout!=l$layout[match(s$sample_id,l$sample_id)]))stop("Run/library layout conflict")
  for(i in seq_len(nrow(l))) {
    vals<-list(sample=l$sample_id[i],mode="complete",strand=l$strand_mode[i])
    qp<-do.call(cfg_path,c(list(c=c,key="quant_dir"),vals))
    for(k in c("allele_counts","thetas","tpms")) {
      op<-do.call(cfg_path,c(list(c=c,key=k),vals))
      required_name<-c(allele_counts="allele_counts.txt",thetas="thetas.txt",tpms="tpms.txt")[[k]]
      if(dirname(op)!=qp||basename(op)!=required_name)
        stop("Pinned quantifier output contract mismatch: ",k)
    }
  }
  # YAML is the sole fragment-setting authority. PE libraries never receive SE settings.
  l[,`:=`(fragment_mean_nt=ifelse(layout=="SE",c$fragments$fragment_mean_nt,NA_real_),
    fragment_sd_nt=ifelse(layout=="SE",c$fragments$fragment_sd_nt,NA_real_),
    fragment_parameter_origin=ifelse(layout=="SE",c$fragments$fragment_parameter_origin,"paired_end_quantifier"))]
  se<-l[layout=="SE"]
  if(nrow(se)&&(any(!is.finite(c(se$fragment_mean_nt,se$fragment_sd_nt)))||
      any(se$fragment_mean_nt<=0|se$fragment_sd_nt<=0)||
      any(!se$fragment_parameter_origin %in% c("assumed_user_defined","measured"))))
    stop("One finite positive SE fragment mean/SD with origin required")
  if(c$profile=="public_cross") {
    if(uniqueN(u$stage)!=1L||uniqueN(u$tissue)!=1L)stop("Separate stage/tissue configs required")
    p<-c(c$parents$focal,c$parents$alternate)
    if(any(!u$maternal_parent %in% p)||any(!u$paternal_parent %in% p)||
       any(u$maternal_parent==u$paternal_parent))stop("Cross parent mismatch")
  }
  if(require_files)for(i in seq_len(nrow(s)))for(m in seq_len(if(s$layout[i]=="PE")2 else 1))
    if(!file.exists(cfg_path(c,"fastq",run=s$run_id[i],mate=m)))stop("FASTQ missing")
  list(sources=s,libraries=l,units=u)
}
