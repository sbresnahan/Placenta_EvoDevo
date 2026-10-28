validate_marker_table <- function(v,c) {
  require_columns(v,c("chrom","genomic_pos","ref","alt","focal_index","alternate_index"))
  if(!nrow(v))stop("No parental markers")
  if(anyNA(v)||any(!v$focal_index %in% 0:1)||any(!v$alternate_index %in% 0:1)||
    any(v$focal_index==v$alternate_index)||any(v$ref==v$alt)||any(nchar(v$ref)!=1|nchar(v$alt)!=1)||
    any(!v$ref %in% c("A","C","G","T"))||any(!v$alt %in% c("A","C","G","T")))
    stop("Markers must be called opposite homozygous biallelic SNPs")
  dup<-v[, .N,by=.(chrom,genomic_pos)][N>1]
  if(nrow(dup))v<-v[!dup,on=.(chrom,genomic_pos)]
  bases<-unlist(lapply(split(seq_len(nrow(v)),ceiling(seq_len(nrow(v))/20000)),function(i)
    as.character(Rsamtools::scanFa(Rsamtools::FaFile(c$paths$genome),
      GenomicRanges::GRanges(v$chrom[i],IRanges::IRanges(v$genomic_pos[i],width=1))))))
  match<-bases==v$ref
  if(mean(match)<c$markers$min_ref_fraction)stop("Marker/genome reference compatibility failed")
  v[match]
}
sample_fields <- function(s,format) {
  k<-strsplit(format,":",fixed=TRUE)[[1]];v<-strsplit(s,":",fixed=TRUE)[[1]]
  setNames(v[seq_along(k)],k)
}
parent_call <- function(fields,mode,min_dp=5,min_fraction=.9) {
  gt<-fields["GT"];allele<-if(!is.na(gt)&&gt %in% c("0/0","0|0"))0L else
    if(!is.na(gt)&&gt %in% c("1/1","1|1"))1L else NA_integer_
  if(is.na(allele))return(NA_integer_)
  if(mode=="mgp")return(if(!is.na(fields["FI"])&&fields["FI"]=="1")allele else NA_integer_)
  ad<-suppressWarnings(as.numeric(strsplit(fields["AD"],",",fixed=TRUE)[[1]]))
  if(length(ad)!=2||any(!is.finite(ad))||any(ad<0)||sum(ad)<min_dp||ad[allele+1]/sum(ad)<min_fraction)
    return(NA_integer_)
  allele
}
parent_calls <- function(values,formats,mode,min_dp=5,min_fraction=.9) {
  # Vectorized within FORMAT layouts: do not run an R function for every site/sample.
  result<-rep(NA_integer_,length(values))
  for(fmt in unique(formats)) {
    ii<-which(formats==fmt);keys<-strsplit(fmt,":",fixed=TRUE)[[1]]
    field<-tstrsplit(values[ii],":",fixed=TRUE)
    take<-function(k) {
      j<-match(k,keys)
      if(is.na(j)||j>length(field))rep(NA_character_,length(ii)) else field[[j]]
    }
    gt<-take("GT");allele<-rep(NA_integer_,length(ii))
    allele[gt %in% c("0/0","0|0")]<-0L
    allele[gt %in% c("1/1","1|1")]<-1L
    if(mode=="mgp") {
      allele[is.na(take("FI"))|take("FI")!="1"]<-NA_integer_
    } else {
      ad<-tstrsplit(take("AD"),",",fixed=TRUE)
      if(length(ad)!=2L)allele[]<-NA_integer_ else {
        ar<-suppressWarnings(as.numeric(ad[[1]]));aa<-suppressWarnings(as.numeric(ad[[2]]))
        depth<-ar+aa
        fraction<-ifelse(allele==0,ar/depth,aa/depth)
        good<-is.finite(ar)&is.finite(aa)&ar>=0&aa>=0&depth>=min_dp&fraction>=min_fraction
        good[is.na(good)]<-FALSE;allele[!good]<-NA_integer_
      }
    }
    result[ii]<-allele
  };result
}
vcf_markers <- function(path,c) {
  con<-gzfile(path,"rt");on.exit(close(con))
  repeat {
    h<-readLines(con,n=1,warn=FALSE)
    if(!length(h))stop("No VCF header")
    if(startsWith(h,"#CHROM"))break
  }
  header<-strsplit(h,"\t",fixed=TRUE)[[1]]
  parent_names<-c(c$markers$focal_sample,c$markers$alternate_sample)
  reference_proxy<-c$markers$mode=="mgp"&&isTRUE(c$markers$proxy)&&
    c$markers$focal_sample=="REFERENCE_C57BL_6J"
  if(any(!parent_names[if(reference_proxy)2 else 1:2] %in% header))
    stop("Missing configured parental VCF sample")
  dupnames<-unlist(c$markers$duplicate_samples)
  if(length(dupnames)&&any(!dupnames %in% header))stop("Missing duplicate-line samples")
  selected<-sort(unique(c(1:9,which(header %in% c(parent_names,dupnames)))))
  kept<-list()
  repeat {
    lines<-readLines(con,n=c$markers$chunk_rows,warn=FALSE);if(!length(lines))break
    if(any(lengths(strsplit(lines,"\t",fixed=TRUE))!=length(header)))stop("Malformed VCF chunk")
    d<-fread(text=paste(lines,collapse="\n"),header=FALSE,sep="\t",quote="",
      colClasses="character",select=selected)
    setnames(d,header[selected])
    d<-d[nchar(REF)==1 & nchar(ALT)==1 & REF!=ALT &
      REF %in% c("A","C","G","T") & ALT %in% c("A","C","G","T")]
    if(c$markers$mode=="panzea")d<-d[FILTER %in% c("PASS",".") & !grepl("(^|;)NI5(;|$)",INFO)]
    if(!nrow(d))next
    getcalls<-function(name)parent_calls(d[[name]],d$FORMAT,c$markers$mode,
        c$markers$min_depth,c$markers$min_allele_fraction)
    a<-if(reference_proxy)rep(0L,nrow(d)) else getcalls(parent_names[1])
    b<-getcalls(parent_names[2])
    good<-!is.na(a)&!is.na(b)&a!=b
    if(length(dupnames))for(k in seq_along(dupnames)) {
      dc<-getcalls(dupnames[k]);target<-if(k==1)a else b
      good<-good & (is.na(dc)|dc==target)
    }
    good[is.na(good)]<-FALSE
    if(any(good))kept[[length(kept)+1L]]<-data.table(chrom=d[["#CHROM"]][good],
      genomic_pos=as.integer(d$POS[good]),ref=d$REF[good],alt=d$ALT[good],
      focal_index=a[good],alternate_index=b[good])
  }
  v<-rbindlist(kept)
  if(!nrow(v))stop("No quality-passing markers in ",basename(path))
  v
}
prepare_markers <- function(c) {
  if(!isTRUE(c$markers$provenance_approved))stop("Parental marker provenance/build review required")
  s<-resolve_manifest_paths(fread(c$paths$marker_sources),"path",c)
  require_columns(s,c("source_id","path","build","sha256"));assert_unique(s,"source_id")
  if(!nrow(s)||any(s$build!=c$reference$build))stop("Marker/reference build mismatch")
  for(i in seq_len(nrow(s))) {
    if(!file.exists(s$path[i])||!identical(sha256(s$path[i]),s$sha256[i]))
      stop("Marker source missing or exact SHA256 not approved: ",s$source_id[i])
  }
  v<-if(c$markers$mode=="prepared")rbindlist(lapply(s$path,fread)) else
    rbindlist(lapply(s$path,vcf_markers,c=c))
  v<-validate_marker_table(v,c)
  setorder(v,chrom,genomic_pos);write_tsv(v,c$paths$markers);v
}
tx_variant_map <- function(v,ex,c,output_vcf,output_map,sample) {
  # v contains genomic SNPs and gt. A typed trio is NOT inferred from this map.
  require_columns(v,c("chrom","genomic_pos","ref","alt","gt"))
  assert_unique(v,c("chrom","genomic_pos"))
  e<-copy(ex);setorder(e,transcript_id,tx_exon_order)
  e[,offset:=c(0,head(cumsum(end-start+1L),-1L)),by=transcript_id]
  h<-GenomicRanges::findOverlaps(GenomicRanges::GRanges(v$chrom,IRanges::IRanges(v$genomic_pos,width=1)),
    GenomicRanges::GRanges(e$chrom,IRanges::IRanges(e$start,e$end)),ignore.strand=TRUE)
  d<-cbind(v[S4Vectors::queryHits(h)],e[S4Vectors::subjectHits(h),
    .(transcript_id,strand,start,end,offset)])
  d[,tx_pos:=as.integer(offset+ifelse(strand=="+",genomic_pos-start+1L,end-genomic_pos+1L))]
  complement<-function(x)chartr("ACGT","TGCA",x)
  d[,`:=`(tx_ref=ifelse(strand=="-",complement(ref),ref),
          tx_alt=ifelse(strand=="-",complement(alt),alt))]
  if(!nrow(d))stop("No exonic personalization SNPs")
  bases<-as.character(Rsamtools::scanFa(Rsamtools::FaFile(c$paths$transcripts),
    GenomicRanges::GRanges(d$transcript_id,IRanges::IRanges(d$tx_pos,width=1))))
  if(any(bases!=d$tx_ref))stop("Transcript REF round-trip failed")
  setorder(d,transcript_id,tx_pos);assert_unique(d,c("transcript_id","tx_pos"))
  header<-c("##fileformat=VCFv4.2",'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",sample,sep="\t"))
  dir.create(dirname(output_vcf),recursive=TRUE,showWarnings=FALSE)
  con<-file(output_vcf,"wt");on.exit(close(con));writeLines(header,con)
  writeLines(paste(d$transcript_id,d$tx_pos,".",d$tx_ref,d$tx_alt,".","PASS",".","GT",d$gt,sep="\t"),con)
  write_tsv(d,output_map);invisible(d)
}
build_reference <- function(c,native=TRUE,public=TRUE) {
  if(!isTRUE(c$reference$annotation_approved))stop("Exact annotation/build approval required")
  for(k in c("genome","annotation"))if(!file.exists(c$paths[[k]]))stop("Missing reference input: ",k)
  for(k in c("genome","annotation"))if(!identical(sha256(c$paths[[k]]),c$reference[[paste0(k,"_sha256")]]))
    stop("Reference checksum not approved: ",k)
  if(native)run_checked(c$tools$samtools,c("faidx",c$paths$genome))
  cs<-fread(c$paths$genome_fai,header=FALSE,select=1:2,col.names=c("chrom","length"))
  write_tsv(cs,c$paths$sizes)
  z<-read_gtf(c$paths$annotation)
  aliases<-c$reference$contig_aliases
  if(length(aliases)) {
    for(k in names(aliases))z[chrom==k,chrom:=as.character(aliases[[k]])]
  }
  # audit_gtf expects a headerless two-column sizes file.
  fwrite(cs,c$paths$sizes,sep="\t",col.names=FALSE)
  a<-audit_gtf(z,c$paths$sizes)
  z[type=="exon",exon_number:=a$exons$tx_exon_order[
    match(paste(transcript_id,start,end),paste(a$exons$transcript_id,a$exons$start,a$exons$end))]]
  z[type=="exon",attributes:=vapply(seq_len(.N),function(i)
    if(grepl('exon_number "',attributes[i],fixed=TRUE))
      sub('exon_number "[^"]*"',paste0('exon_number "',exon_number[i],'"'),attributes[i])
    else paste0(attributes[i],' exon_number "',exon_number[i],'";'),character(1))]
  fwrite(z[,.(chrom,source,type,start,end,score,strand,frame,attributes)],c$paths$gtf,
    sep="\t",quote=FALSE,col.names=FALSE)
  write_tsv(a$exons,c$paths$exons);write_tsv(a$transcripts,c$paths$tx_map)
  write_tsv(a$audit,c$paths$annotation_audit)
  fa<-Rsamtools::FaFile(c$paths$genome);open(fa);on.exit(close(fa))
  ids<-sort(a$transcripts$transcript_id);unlink(c$paths$transcripts)
  for(i in seq(1L,length(ids),by=1000L)) {
    e<-a$exons[transcript_id %in% ids[i:min(i+999L,length(ids))]]
    setorder(e,transcript_id,tx_exon_order)
    seqs<-Rsamtools::scanFa(fa,GenomicRanges::GRanges(e$chrom,IRanges::IRanges(e$start,e$end)))
    minus<-which(e$strand=="-")
    if(length(minus))seqs[minus]<-Biostrings::reverseComplement(seqs[minus])
    e[,sequence:=as.character(seqs)]
    t<-e[,.(sequence=paste0(sequence,collapse="")),by=transcript_id]
    x<-Biostrings::DNAStringSet(t$sequence);names(x)<-t$transcript_id
    if(any(Biostrings::width(x)!=a$transcripts$spliced_length[match(names(x),a$transcripts$transcript_id)]))
      stop("Transcript length mismatch")
    Biostrings::writeXStringSet(x,c$paths$transcripts,append=TRUE)
  }
  if(native)run_checked(c$tools$samtools,c("faidx",c$paths$transcripts))
  maps<-events<-list()
  for(w in c(25L,50L,100L)) {
    g<-group_promoters(a$transcripts,w)
    g[,`:=`(tss_unknown_repeat=FALSE,tss_te=FALSE,te_class=NA_character_,
      te_family=NA_character_,window_te=FALSE,first_exon_te=FALSE)]
    ev<-event_catalog(g,a$exons);ev[,`:=`(window=w,tss_te_label="not_assessed")]
    maps[[length(maps)+1L]]<-g[,.(transcript_id,gene_id,event_id,window=w)]
    events[[length(events)+1L]]<-ev
  }
  write_tsv(rbindlist(maps),c$paths$event_map);write_tsv(rbindlist(events),c$paths$events)
  # Human primary pair is structural, not a TE exposure comparison.
  ev<-rbindlist(events)[window==50 & distinguishing_bp>0]
  p<-rbindlist(lapply(split(ev,ev$gene_id),function(d) {
    if(nrow(d)<2)return(NULL)
    ab<-t(combn(sort(d$event_id),2))
    q<-data.table(gene_id=d$gene_id[1],gene_name=d$gene_name[1],event_a=ab[,1],event_b=ab[,2],
      pair_type="structural_promoter_pair",
      support=pmin(d$distinguishing_bp[match(ab[,1],d$event_id)],d$distinguishing_bp[match(ab[,2],d$event_id)]))
    setorder(q,-support,event_a,event_b);q[1]
  }),fill=TRUE)
  if(!ncol(p))p<-data.table(gene_id=character(),gene_name=character(),event_a=character(),event_b=character(),pair_type=character())
  write_tsv(p,c$paths$pairs)
  if(public&&c$profile=="public_cross") {
    v<-prepare_markers(c);v[,gt:="0/1"]
    tx_variant_map(v,a$exons,c,c$paths$tx_vcf,c$paths$tx_marker_map,c$cohort_id)
    make_index(c,c$paths$tx_vcf,c$cohort_id,c$paths$personal_fasta,c$paths$index)
  }
  write_tsv(data.table(endpoint=c("human_whole_gene_masking","TE_comparison"),status="deferred_unimplemented"),
    c$paths$deferred)
}
