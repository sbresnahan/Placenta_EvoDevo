# Adapter for the already-phased GUSTO transcriptomic VCF. The source VCF,
# chromosome origin matrix and linked FASTQs are private, external inputs.
gusto_enabled <- function(c) c$profile == "human_trio" && isTRUE(c$gusto_precomputed$enabled)

# The upstream transcriptomic VCF also contains boundary-clipped deletion
# placeholders (e.g. POS=0, REF=., ALT=7). These cannot be represented as
# sequence-resolved VCF variants. The Seadragon bcftools 1.16 build evaluates
# REF/ALT regex predicates as false even for valid alleles, so use validated
# scalar predicates here. SNPs AND indels remain eligible. The separate source
# audit checks allele sequence content, and `bcftools norm -c e` checks REF.
# GT="alt" MUST be evaluated after a separate `bcftools view -s` invocation;
# combining the two options filters all children before subsetting to one.
gusto_valid_variant_filter <- function() {
  'GT="alt" && POS>=1 && REF!="." && ALT[0]!="."'
}

gusto_vcf_child_id <- function(column_name) {
  if(anyNA(column_name) || any(!grepl("^.+_B[0-9]{3}-[0-9]{5}$",column_name)))
    stop("GUSTO VCF sample columns must contain a delimited B infant ID")
  sub("^.*_","",column_name)
}

gusto_vcf_sample <- function(child_id,columns) {
  if(length(child_id)!=1L || is.na(child_id) ||
     !grepl("^B[0-9]{3}-[0-9]{5}$",child_id))stop("Invalid GUSTO infant ID")
  children<-gusto_vcf_child_id(columns)
  if(anyDuplicated(children))stop("Repeated infant identity in GUSTO VCF columns")
  matched<-columns[children==child_id]
  if(length(matched)!=1L)stop("Exact GUSTO infant ID absent from VCF: ",child_id)
  matched # bcftools requires the literal composite column, not its normalized identity.
}

gusto_vcf_samples <- function(path) {
  if(!file.exists(path) || file.info(path)$size<=0)stop("Missing precomputed GUSTO VCF: ",path)
  con<-if(grepl("\\.gz$",path))gzfile(path,"rt") else file(path,"rt")
  on.exit(close(con))
  repeat {
    line<-readLines(con,n=1L,warn=FALSE)
    if(!length(line))stop("GUSTO VCF has no #CHROM header")
    if(startsWith(line,"#CHROM"))break
  }
  h<-strsplit(line,"\t",fixed=TRUE)[[1]]
  if(length(h)<10L || !identical(h[1:9],c("#CHROM","POS","ID","REF","ALT",
    "QUAL","FILTER","INFO","FORMAT")))stop("Invalid GUSTO VCF header")
  if(anyDuplicated(h[-(1:9)]))stop("Duplicated GUSTO VCF samples")
  h[-(1:9)]
}

gusto_validate_linkage <- function(c,libraries) {
  require_columns(libraries,c("sample_id","child_id"))
  assert_unique(libraries,"child_id","GUSTO RNA-to-genotype linkage")
  m<-data.table::fread(c$gusto_precomputed$fastq_id_map,header=FALSE,
    col.names=c("infant_id","fastq_id"),colClasses="character")
  if(ncol(m)!=2L || !nrow(m))stop("Invalid infant FASTQ ID table")
  assert_unique(m,"infant_id");assert_unique(m,"fastq_id")
  m[,`:=`(child_id=paste0("B",infant_id),sample_id=paste0("J",fastq_id))]
  if(anyDuplicated(m[,.(child_id,sample_id)]) || any(!libraries$sample_id %in% m$sample_id))
    stop("GUSTO sample absent from FASTQ ID mapping")
  linked<-merge(libraries,m[,.(sample_id,child_id)],by="sample_id",all.x=TRUE,
    suffixes=c("","_mapping"),sort=FALSE)
  if(anyNA(linked$child_id_mapping) || any(linked$child_id!=linked$child_id_mapping))
    stop("GUSTO RNA-to-infant mapping disagrees with library manifest")
  vcf_samples<-gusto_vcf_samples(c$gusto_precomputed$vcf)
  invisible(vapply(libraries$child_id,gusto_vcf_sample,character(1),columns=vcf_samples))
  origins<-data.table::fread(c$gusto_precomputed$origin_matrix,nrows=1L)
  if(!"chr" %in% names(origins))stop("GUSTO origin matrix requires chr column")
  missing<-setdiff(c(paste0(libraries$child_id,"_L"),paste0(libraries$child_id,"_R")),names(origins))
  if(length(missing))stop("Missing GUSTO maternal/paternal origin columns: ",paste(missing,collapse=","))
  invisible(TRUE)
}

gusto_collision_audit <- function(vcf)paste0(vcf,".excluded_ambiguous.tsv")

gusto_child_map <- function(query,output,vcf=NULL) {
  z<-data.table::fread(query,header=FALSE,
    col.names=c("transcript_id","tx_pos","variant_id","tx_ref","tx_alt","gt"),
    colClasses=c("character","integer","character","character","character","character"))
  if(!nrow(z))stop("No alternate-bearing GUSTO transcript variants for this child")
  if(anyNA(z))stop("Missing values in GUSTO child transcript variant query")
  valid_gt<-c("0|1","1|0","1|1","0/1","1/0","1/1")
  if(any(!z$gt %in% valid_gt))
    stop("GUSTO child query contains missing or reference-only genotypes: ",
      paste(unique(z$gt[!z$gt %in% valid_gt]),collapse=","))
  if(any(!grepl("^chr([0-9]+|X|Y|M|MT):[0-9]+:[ACGT]+:[ACGT]+(v[0-9]+)?$",z$variant_id)))
    stop("Invalid genomic locus ID in GUSTO transcript variant query")
  # Identify conflicts in this child's actual ALT-bearing variants, including
  # indels. A site can have multiple genomic ALT alleles even if its transcript
  # representations are distributed across different isoforms. All copies of
  # such a genomic site must be excluded together for consistent personalizing
  # and RNA-count orientation. Repeated transcript coordinates also trigger a
  # site-wide exclusion, even if their IDs appear equivalent.
  locus<-data.table::tstrsplit(sub("v[0-9]+$","",z$variant_id),":",fixed=TRUE)
  genomic_site<-paste(locus[[1]],locus[[2]],sep=":")
  alleles<-paste(locus[[3]],locus[[4]],sep=">")
  multiallelic_sites<-data.table::data.table(site=genomic_site,allele=alleles)[
    ,.(n_alleles=data.table::uniqueN(allele)),by=site][n_alleles>1L,site]
  tx_key<-paste(z$transcript_id,z$tx_pos,sep="\t")
  repeated_tx<-z[,.(n=.N),by=.(transcript_id,tx_pos)][n>1L]
  repeated_key<-paste(repeated_tx$transcript_id,repeated_tx$tx_pos,sep="\t")
  ambiguous_sites<-union(multiallelic_sites,unique(genomic_site[tx_key %in% repeated_key]))
  if(length(ambiguous_sites) && is.null(vcf))
    stop("Ambiguous GUSTO genomic sites require the child VCF for exclusion")
  if(!is.null(vcf)) {
    key<-tx_key
    excluded<-genomic_site %in% ambiguous_sites
    lines<-readLines(vcf,warn=FALSE)
    is_header<-startsWith(lines,"#")
    records<-strsplit(lines[!is_header],"\t",fixed=TRUE)
    if(!any(startsWith(lines,"#CHROM")) || length(records)!=nrow(z) ||
       any(lengths(records)<10L))
      stop("GUSTO child VCF and query have inconsistent records")
    vcf_key<-vapply(records,function(x)paste(x[1:2],collapse="\t"),character(1))
    if(!identical(vcf_key,key))stop("GUSTO child VCF/query transcript positions differ")
    audit<-data.table::copy(z[excluded])
    audit[,`:=`(genomic_site=genomic_site[excluded],
      reason="conflicting_alt_or_repeated_transcript_position")]
    atomic_tsv(audit,gusto_collision_audit(vcf))
    if(any(excluded)) {
      # Both outputs are regenerated when a failed child is resumed. Use temp
      # files in the same directory so individual replacements are atomic.
      keep<-is_header
      keep[!is_header]<-!excluded
      tmp_vcf<-tempfile(".gusto-filter-",tmpdir=dirname(vcf))
      tmp_query<-tempfile(".gusto-filter-",tmpdir=dirname(query))
      writeLines(lines[keep],tmp_vcf)
      data.table::fwrite(z[!excluded],tmp_query,sep="\t",col.names=FALSE)
      if(!file.rename(tmp_vcf,vcf) || !file.rename(tmp_query,query))
        stop("Failed to replace filtered GUSTO child VCF/query")
      message("GUSTO: excluded ",length(ambiguous_sites)," ambiguous genomic sites (",
        sum(excluded)," VCF records); audit: ",gusto_collision_audit(vcf))
      z<-z[!excluded]
    }
  }
  # Keep indels in the personalized VCF/index, but use only SNVs as direct
  # genomic allelic-count anchors. Indel normalization across exon boundaries
  # requires a separate allele-aware mapping contract.
  z<-z[nchar(tx_ref)==1L & nchar(tx_alt)==1L &
    tx_ref %in% c("A","C","G","T") & tx_alt %in% c("A","C","G","T")]
  if(!nrow(z))stop("No SNV anchors in the personalized GUSTO transcript VCF")
  z[,locus_id:=sub("v[0-9]+$","",variant_id)]
  parts<-data.table::tstrsplit(z$locus_id,":",fixed=TRUE)
  z[,`:=`(chrom=parts[[1]],genomic_pos=as.integer(parts[[2]]))]
  reverse_base<-function(a)chartr("ACGT","TGCA",a)
  strand_compatible<-(parts[[3]]==z$tx_ref & parts[[4]]==z$tx_alt) |
    (parts[[3]]==reverse_base(z$tx_ref) & parts[[4]]==reverse_base(z$tx_alt))
  if(any(!strand_compatible))stop("VCF genomic ID alleles disagree with transcript REF/ALT")
  if(anyNA(z$genomic_pos)||any(z$genomic_pos<1L)||any(z$tx_pos<1L))
    stop("Invalid genomic or transcript position in GUSTO VCF")
  assert_unique(z,c("transcript_id","tx_pos"),"GUSTO transcript VCF variant positions")
  data.table::setorder(z,transcript_id,tx_pos)
  write_tsv(z,output)
}

gusto_orient_counts <- function(mapped,c,sample,child_id) {
  anchors<-data.table::fread(cfg_path(c,"child_tx_map",sample=sample))
  require_columns(anchors,c("transcript_id","tx_pos","tx_ref","tx_alt","gt",
    "locus_id","chrom","genomic_pos"))
  assert_unique(anchors,c("transcript_id","tx_pos"),"GUSTO child transcript SNPs")
  origins<-data.table::fread(c$gusto_precomputed$origin_matrix,
    select=c("chr",paste0(child_id,"_L"),paste0(child_id,"_R")),
    colClasses="character")
  data.table::setnames(origins,c("chr","left_origin","right_origin"))
  assert_unique(origins,"chr","GUSTO chromosome origin matrix")
  valid<-c("maternal","paternal","undefined")
  if(anyNA(origins)||any(!origins$left_origin %in% valid) ||
     any(!origins$right_origin %in% valid) ||
     any((origins$left_origin=="maternal")!=(origins$right_origin=="paternal")) ||
     any((origins$left_origin=="paternal")!=(origins$right_origin=="maternal")))
    stop("Noncomplementary or invalid precomputed GUSTO origin labels")
  a<-merge(anchors,origins,by.x="chrom",by.y="chr",all.x=TRUE,sort=FALSE)
  a[,`:=`(vcf_chrom=chrom,vcf_genomic_pos=genomic_pos)]
  allele_parts<-data.table::tstrsplit(a$locus_id,":",fixed=TRUE)
  a[,`:=`(vcf_genomic_ref=allele_parts[[3]],vcf_genomic_alt=allele_parts[[4]])]
  a[,c("chrom","genomic_pos"):=NULL]
  a[,left_allele:=substr(gt,1L,1L)]
  a[,right_allele:=substr(gt,3L,3L)]
  a[,phased_het:=grepl("^[01]\\|[01]$",gt) & left_allele!=right_allele]
  a[,maternal_index:=data.table::fifelse(phased_het & left_origin=="maternal",
    as.integer(left_allele),data.table::fifelse(phased_het & right_origin=="maternal",
      as.integer(right_allele),NA_integer_))]
  key<-c("transcript_id","tx_pos","tx_ref","tx_alt")
  a<-a[,c(key,"vcf_chrom","vcf_genomic_pos","vcf_genomic_ref","vcf_genomic_alt",
    "locus_id","maternal_index",
    "left_origin","phased_het"),with=FALSE]
  out<-merge(mapped,a,by=key,all.x=TRUE,sort=FALSE)
  if(any(!is.na(out$vcf_genomic_pos) & is.na(out$genomic_pos)))
    stop("GUSTO transcript SNP is absent from the configured GTF exon map")
  valid_coordinate<-!is.na(out$vcf_genomic_pos) & !is.na(out$genomic_pos)
  if(any(out$genomic_pos[valid_coordinate]!=out$vcf_genomic_pos[valid_coordinate] |
         out$chrom[valid_coordinate]!=out$vcf_chrom[valid_coordinate] |
         out$ref[valid_coordinate]!=out$vcf_genomic_ref[valid_coordinate] |
         out$alt[valid_coordinate]!=out$vcf_genomic_alt[valid_coordinate]))
    stop("GUSTO VCF genomic ID alleles/position do not match the transcript GTF")
  out[,parent_status:=data.table::fifelse(is.na(locus_id),"no_precomputed_variant",
    data.table::fifelse(!phased_het,"child_not_phased_heterozygous",
      data.table::fifelse(is.na(maternal_index),"origin_undefined","precomputed_origin")))]
  out[,`:=`(ref_maternal=ifelse(is.na(maternal_index),NA,maternal_index==0L),
    maternal_count=ifelse(is.na(maternal_index),NA_real_,
      ifelse(maternal_index==0L,ref_count,alt_count)),
    paternal_count=ifelse(is.na(maternal_index),NA_real_,
      ifelse(maternal_index==0L,alt_count,ref_count)))]
  out[is.na(locus_id),locus_id:=paste(chrom,genomic_pos,ref,alt,sep=":")]
  out[,c("vcf_chrom","vcf_genomic_pos","vcf_genomic_ref","vcf_genomic_alt",
    "maternal_index","left_origin","phased_het"):=NULL]
  out
}
