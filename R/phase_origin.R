# Chromosome-wise origin of the first/second phased fetal haplotype.
# Operates on one trio's normalized, biallelic typed.query.tsv: no random
# SNP subsampling, no cohort-wide VCF kept in memory, and no new imputation.
phase_summary_path <- function(c,sample)paste0(cfg_path(c,"child_query",sample=sample),".phase_origin.tsv")

phase_origin_settings <- function(c) {
  x<-c$phase_origin
  if(is.null(x))x<-list()
  list(enabled=isTRUE(x$enabled),
    min_informative=if(is.null(x$min_informative))20L else as.integer(x$min_informative),
    min_agreement=if(is.null(x$min_agreement))0.95 else as.numeric(x$min_agreement),
    window_bp=if(is.null(x$window_bp))5000000L else as.integer(x$window_bp),
    min_window_votes=if(is.null(x$min_window_votes))5L else as.integer(x$min_window_votes),
    min_covered_windows=if(is.null(x$min_covered_windows))3L else as.integer(x$min_covered_windows))
}

split_diploid_gt <- function(gt) {
  z<-sub(":.*$","",as.character(gt))
  good<-grepl("^[01][/|][01]$",z)
  first<-second<-rep(NA_integer_,length(z))
  first[good]<-as.integer(substr(z[good],1L,1L))
  second[good]<-as.integer(substr(z[good],3L,3L))
  data.table::data.table(first=first,second=second,phased=good&grepl("|",z,fixed=TRUE))
}

infer_chromosome_origin <- function(anchors,settings) {
  require_columns(anchors,c("chrom","genomic_pos","child_gt","mother_gt","father_gt"))
  if(!nrow(anchors))stop("No phased trio variants supplied")
  ch<-split_diploid_gt(anchors$child_gt)
  mo<-split_diploid_gt(anchors$mother_gt)
  fa<-split_diploid_gt(anchors$father_gt)
  # A homozygous parent identifies the matching fetal haplotype at an
  # independently informative child-heterozygous SNP. Parental phase is
  # irrelevant; phased fetal haplotypes must remain consistent by chromosome.
  ok<-ch$phased & !is.na(ch$first) & !is.na(ch$second) & ch$first!=ch$second
  mom_hom<-!is.na(mo$first)&!is.na(mo$second)&mo$first==mo$second
  dad_hom<-!is.na(fa$first)&!is.na(fa$second)&fa$first==fa$second
  mv<-ifelse(ok&mom_hom,ifelse(ch$first==mo$first,1L,2L),NA_integer_)
  dv<-ifelse(ok&dad_hom,ifelse(ch$first==fa$first,2L,1L),NA_integer_)
  discordant_parent<-!is.na(mv)&!is.na(dv)&mv!=dv
  vote<-ifelse(discordant_parent,NA_integer_,ifelse(!is.na(mv),mv,dv))
  # Keep only score vectors, not a second whole-genome copy of the trio query.
  x<-data.table::data.table(chrom=anchors$chrom,genomic_pos=anchors$genomic_pos,
    vote=vote,discordant_parent=discordant_parent,phased_het=ok)
  x[,window:=floor((as.numeric(genomic_pos)-1)/settings$window_bp)]
  ans<-x[,{
    vv<-vote[!is.na(vote)]
    n1<-sum(vv==1L);n2<-sum(vv==2L)
    winner<-if(n1>=n2)1L else 2L
    agreement<-if(length(vv))max(n1,n2)/length(vv) else NA_real_
    blocks<-.SD[!is.na(vote),.(n=.N,n1=sum(vote==1L),n2=sum(vote==2L)),by=window]
    covered<-sum(blocks$n>=settings$min_window_votes)
    # Phase-set switches, including switching between unlabelled PS blocks,
    # invalidate chromosome-wide parental assignment instead of silently
    # reversing parent labels over part of the chromosome.
    conflicting_blocks<-if(nrow(blocks))sum(blocks$n>=settings$min_window_votes &
      if(winner==1L)blocks$n2/blocks$n>=0.75 else blocks$n1/blocks$n>=0.75) else 0L
    status<-if(length(vv)<settings$min_informative)"insufficient_informative_sites" else
      if(covered<settings$min_covered_windows)"insufficient_chromosome_coverage" else
      if(agreement<settings$min_agreement)"inconsistent_phase_or_trio" else
      if(conflicting_blocks>0L)"phase_switch_suspected" else "callable"
    list(maternal_haplotype=if(status=="callable")winner else NA_integer_,
      n_informative=length(vv),n_maternal_hap1_votes=n1,n_maternal_hap2_votes=n2,
      agreement=agreement,n_phased_het=sum(phased_het),
      n_parent_discordant=sum(discordant_parent),
      n_covered_windows=covered,n_conflicting_windows=conflicting_blocks,status=status)
  },by=chrom]
  data.table::setorder(ans,chrom)
  ans
}

orient_counts_chromosome <- function(mapped,anchors,phase) {
  # Vectorized Mendelian assignment is a conservative fallback for uncallable
  # chromosomes / unphased variants. This avoids an R loop over every imputed
  # trio SNP, which is prohibitively slow for a whole-genome GUSTO trio.
  a<-data.table::copy(anchors)
  require_columns(a,c("chrom","genomic_pos","ref","alt","child_gt"))
  data.table::setDT(a)
  keys<-c("chrom","genomic_pos","ref","alt")
  assert_unique(a,keys,"Phased trio anchors")
  child<-split_diploid_gt(a$child_gt)
  mother<-split_diploid_gt(a$mother_gt)
  father<-split_diploid_gt(a$father_gt)
  hetero<-!is.na(child$first)&!is.na(child$second)&child$first!=child$second
  mh<-!is.na(mother$first)&!is.na(mother$second)&mother$first==mother$second
  fh<-!is.na(father$first)&!is.na(father$second)&father$first==father$second
  mendelian<-hetero&(mh|fh)&!(mh&fh&mother$first==father$first)
  di<-ifelse(mh,mother$first,ifelse(fh,1L-father$first,NA_integer_))
  a[,`:=`(direct_index=ifelse(mendelian,di,NA_integer_),
    direct_status=ifelse(!hetero,"child_not_heterozygous",
      ifelse(mh&fh&mother$first==father$first,"mendelian_conflict",
        ifelse(mendelian,"typed_trio","ambiguous_transmission"))))]
  a<-merge(a,phase[,.(chrom,maternal_haplotype,status)],by="chrom",all.x=TRUE,sort=FALSE)
  gt<-split_diploid_gt(a$child_gt)
  a[,maternal_index:=ifelse(!is.na(maternal_haplotype)&gt$phased&
    !is.na(gt$first)&!is.na(gt$second)&gt$first!=gt$second,
    ifelse(maternal_haplotype==1L,gt$first,gt$second),NA_integer_)]
  a[,phase_status:=status]
  add<-merge(mapped,a[,.(chrom,genomic_pos,ref,alt,maternal_index,phase_status,
    direct_index,direct_status)],by=keys,all.x=TRUE,sort=FALSE)
  add[,`:=`(ref_maternal=ifelse(is.na(direct_index),NA,direct_index==0L),
    maternal_count=ifelse(is.na(direct_index),NA_real_,
      ifelse(direct_index==0L,ref_count,alt_count)),
    paternal_count=ifelse(is.na(direct_index),NA_real_,
      ifelse(direct_index==0L,alt_count,ref_count)),
    parent_status=ifelse(is.na(direct_index),direct_status,"typed_trio"))]
  add[is.na(direct_status),parent_status:="no_typed_anchor"]
  phased<-which(!is.na(add$maternal_index))
  if(length(phased)) {
    contradictory<-phased[add$parent_status[phased]=="typed_trio" &
      !is.na(add$ref_maternal[phased]) &
      add$ref_maternal[phased]!=(add$maternal_index[phased]==0L)]
    good<-setdiff(phased,contradictory)
    add[good,`:=`(ref_maternal=maternal_index==0L,
      maternal_count=ifelse(maternal_index==0L,ref_count,alt_count),
      paternal_count=ifelse(maternal_index==0L,alt_count,ref_count),
      parent_status="phased_chromosome_trio")]
    add[contradictory,`:=`(ref_maternal=NA,maternal_count=NA_real_,
      paternal_count=NA_real_,parent_status="phase_vs_mendelian_conflict")]
  }
  add[,c("maternal_index","phase_status","direct_index","direct_status"):=NULL]
  add[,locus_id:=paste(chrom,genomic_pos,ref,alt,sep=":")]
  add
}
