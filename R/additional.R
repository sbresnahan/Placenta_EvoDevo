# Additive, read-only reanalysis of completed primary POE counts.
# Never writes the frozen reference, pilot, quantification or model endpoints.
additional_dir <- function(cfg) file.path(dirname(cfg$paths$model_gene),"additional_v1")
additional_paths <- function(cfg) {
  d<-additional_dir(cfg)
  c(genes=file.path(d,"gene_POE.tsv"),
    structural=file.path(d,"separated_TSS_donor_POE.tsv"),
    annotation=file.path(d,"event_TE_annotation.tsv"),
    provenance=file.path(d,"analysis_provenance.tsv"))
}
additional_receipt <- function(cfg)file.path(cfg$paths$receipts,"additional_v1.rds")
additional_settings <- function(cfg) {
  opt<-cfg$additional_analysis
  if(is.null(opt))opt<-list()
  list(min_depth=if(is.null(opt$min_depth))10L else as.integer(opt$min_depth),
    min_per_direction=if(is.null(opt$min_per_direction))2L else as.integer(opt$min_per_direction),
    max_pairs_per_gene=if(is.null(opt$max_pairs_per_gene))10L else as.integer(opt$max_pairs_per_gene),
    te_rmsk=opt$te_rmsk)
}
additional_validate <- function(cfg,settings) {
  if(any(!is.finite(unlist(settings[c("min_depth","min_per_direction","max_pairs_per_gene")]))))
    stop("Nonfinite additional analysis settings")
  if(settings$min_depth<1||settings$min_per_direction<2||settings$max_pairs_per_gene<1)
    stop("Invalid additional analysis settings")
  if(!is.null(settings$te_rmsk)&&(!file.exists(settings$te_rmsk)||
      file.info(settings$te_rmsk)$size<=0))stop("TE annotation missing or empty")
}
additional_inference_labels <- function(x,scale="log2") {
  if(!"estimate"%in%names(x))x[,estimate:=NA_real_]
  if(!"se"%in%names(x))x[,se:=NA_real_]
  if(!"residual_df"%in%names(x))x[,residual_df:=NA_real_]
  if(!"p_value"%in%names(x))x[,p_value:=NA_real_]
  margin<-if(scale=="log2")log2(1.25) else log(1.25)
  x[,`:=`(nominal_p_lt_0_05=!is.na(p_value)&p_value<.05,
    ci90_low=ifelse(is.finite(estimate)&is.finite(se)&residual_df>0,
      estimate-qt(.95,residual_df)*se,NA_real_),
    ci90_high=ifelse(is.finite(estimate)&is.finite(se)&residual_df>0,
      estimate+qt(.95,residual_df)*se,NA_real_))]
  x[,equivalent_within_1_25_fold:=ifelse(is.finite(ci90_low)&is.finite(ci90_high),
    ci90_low> -margin&ci90_high<margin,NA)]
  x
}

# Freeze a SNP by maximum count of callable *independent biological units*,
# then minimum depth; never inspect the observed maternal imbalance to select it.
additional_locus <- function(d,group,units,settings) {
  require_columns(d,c(group,"locus_id","unit_id","maternal_count","paternal_count"))
  if(!nrow(d))return(d[0])
  z<-copy(d)
  z[,depth:=maternal_count+paternal_count]
  z<-z[is.finite(depth)&depth>=settings$min_depth]
  if(!nrow(z))return(z)
  assert_unique(z,c(group,"unit_id","locus_id"),"Additional analysis entity/locus/unit")
  scored<-z[,.(n_units=uniqueN(unit_id),min_depth=min(depth)),by=c(group,"locus_id")]
  setorderv(scored,c(group,"n_units","min_depth","locus_id"),c(1L,-1L,-1L,1L))
  selected<-scored[,.SD[1],by=group]
  merge(z,selected[,c(group,"locus_id"),with=FALSE],by=c(group,"locus_id"))
}
additional_public_genes <- function(counts,units,cfg,settings) {
  picked<-additional_locus(counts,"gene_id",units,settings)
  independent<-!anyNA(units$cluster_id)&&!anyDuplicated(units$cluster_id)&&
    all(nzchar(units$cluster_id))
  rows<-lapply(split(picked,picked$gene_id),function(d) {
    u<-units[match(d$unit_id,unit_id)]
    orient<-ifelse(u$maternal_parent==cfg$parents$focal,1L,-1L)
    n1<-sum(orient==1L);n2<-sum(orient==-1L)
    out<-data.table(gene_id=d$gene_id[1],locus_id=d$locus_id[1],
      n_focal_mother=n1,n_alternate_mother=n2,
      min_depth=min(d$depth),n_units=nrow(d),
      status="insufficient_reciprocal_units")
    if(n1<settings$min_per_direction||n2<settings$min_per_direction)return(out)
    y<-log2((d$maternal_count+.5)/(d$paternal_count+.5))-
      log2(cfg$dosage$maternal/cfg$dosage$paternal)
    fit<-public_fit(y,orient,independent)
    cbind(out[,status:=NULL],fit[,setdiff(names(fit),"n_units"),with=FALSE])
  })
  result<-rbindlist(rows,fill=TRUE)
  if(!ncol(result))result<-data.table(gene_id=character(),status=character(),p_value=numeric())
  if(!"p_value"%in%names(result))result[,p_value:=NA_real_]
  result[,q_value:=p.adjust(p_value,"BH")]
  additional_inference_labels(result,"log2")
  result[,`:=`(effect_type="gene_parent_of_origin",
    evidence="Ornaments_expected_counts_single_frozen_locus",
    hypothesis="maternal_vs_paternal_dosage_adjusted_reciprocal_intercept")]
  result
}

additional_event_pairs <- function(events,settings) {
  e<-copy(events[window==50])
  setorder(e,gene_id,event_id)
  rows<-lapply(split(e,e$gene_id),function(d) {
    if(nrow(d)<2L)return(NULL)
    comb<-t(combn(seq_len(nrow(d)),2))
    a<-d[comb[,1]];b<-d[comb[,2]]
    known_donor<-!is.na(a$first_donor)&!is.na(b$first_donor)
    equal_donor<-known_donor&a$first_donor==b$first_donor
    close_tss<-abs(a$tss_min-b$tss_min)<=50L
    # TSS contrast conditions on identical first donor; first-donor contrast
    # conditions on local TSS (~50bp). Neither is a pure causal perturbation.
    kind<-ifelse(equal_donor&!close_tss,"TSS_conditional_on_first_donor",
      ifelse(known_donor&!equal_donor&close_tss,"first_donor_conditional_on_TSS",NA_character_))
    out<-data.table(gene_id=a$gene_id,gene_name=a$gene_name,
      event_a=a$event_id,event_b=b$event_id,
      effect_type=kind,tss_distance_bp=abs(a$tss_min-b$tss_min),
      donor_distance_bp=abs(a$first_donor-b$first_donor),
      event_a_distinguishing_bp=a$distinguishing_bp,
      event_b_distinguishing_bp=b$distinguishing_bp)
    out[!is.na(effect_type)]
  })
  rbindlist(rows,fill=TRUE)
}

additional_public_pairs <- function(event_counts,events,units,cfg,settings) {
  pairs<-additional_event_pairs(events,settings)
  template<-data.table(gene_id=character(),gene_name=character(),event_a=character(),
    event_b=character(),effect_type=character(),status=character(),p_value=numeric())
  if(!nrow(pairs))return(template)
  # Bound high-isoform genes before expensive locus joins, on structural
  # coordinate support alone (never observed maternal/paternal imbalance).
  setorder(pairs,gene_id,effect_type,-event_a_distinguishing_bp,
    -event_b_distinguishing_bp,event_a,event_b)
  pairs<-pairs[,.SD[seq_len(min(.N,settings$max_pairs_per_gene*5L))],
    by=.(gene_id,effect_type)]
  raw<-copy(event_counts)
  raw[,depth:=maternal_count+paternal_count]
  raw<-raw[is.finite(depth)&depth>=settings$min_depth]
  assert_unique(raw,c("event_id","unit_id","locus_id"),"Event support per biological unit")
  independent<-!anyNA(units$cluster_id)&&!anyDuplicated(units$cluster_id)&&
    all(nzchar(units$cluster_id))
  by_event<-split(raw,raw$event_id)
  rows<-lapply(seq_len(nrow(pairs)),function(i) {
    p<-pairs[i]
    a<-by_event[[p$event_a]];b<-by_event[[p$event_b]]
    result<-copy(p)[,`:=`(status="no_shared_callable_locus",n_units=0L,
      locus_id=NA_character_,n_focal_mother=0L,n_alternate_mother=0L,
      direct_coordinate_both=FALSE)]
    if(is.null(a)||is.null(b))return(result)
    ab<-merge(a,b,by=c("unit_id","locus_id"),suffixes=c("_a","_b"),sort=FALSE)
    if(!nrow(ab))return(result)
    loci<-ab[,.(n_units=uniqueN(unit_id),min_depth=min(depth_a,depth_b)),by=locus_id]
    setorder(loci,-n_units,-min_depth,locus_id)
    d<-ab[locus_id==loci$locus_id[1]]
    u<-units[match(d$unit_id,unit_id)]
    orientation<-ifelse(u$maternal_parent==cfg$parents$focal,1L,-1L)
    n1<-sum(orientation==1L);n2<-sum(orientation==-1L)
    result[,`:=`(n_units=nrow(d),locus_id=loci$locus_id[1],
      n_focal_mother=n1,n_alternate_mother=n2,
      direct_coordinate_both=all(d$distinguishing_coordinate_a%in%TRUE&
        d$distinguishing_coordinate_b%in%TRUE),
      status="insufficient_reciprocal_units")]
    if(n1<settings$min_per_direction||n2<settings$min_per_direction)return(result)
    contrast<-log2((d$maternal_count_a+.5)/(d$paternal_count_a+.5))-
      log2((d$maternal_count_b+.5)/(d$paternal_count_b+.5))
    fit<-public_fit(contrast,orientation,independent)
    cbind(result[,status:=NULL],fit[,setdiff(names(fit),"n_units"),with=FALSE])
  })
  out<-rbindlist(rows,fill=TRUE)
  # Limit multiple predeclared comparisons per gene using callability only.
  setorder(out,gene_id,effect_type,-n_units,event_a,event_b)
  out<-out[,.SD[seq_len(min(.N,settings$max_pairs_per_gene))],
    by=.(gene_id,effect_type)]
  if(!"p_value"%in%names(out))out[,p_value:=NA_real_]
  out[,q_value:=p.adjust(p_value,"BH")]
  additional_inference_labels(out,"log2")
  out[,`:=`(evidence="Ornaments_expected_counts_same_locus_not_raw_molecules",
    hypothesis="reciprocal_intercept_difference_of_allelic_ratios")]
  out
}

additional_human_genes <- function(counts,cfg,settings) {
  # Exactly one representative per family and one orientation-independent
  # sentinel locus per gene. The baseline promoter model remains untouched.
  require_columns(counts,c("sample_id","family_id","gene_id","locus_id",
    "ref_maternal","maternal_count","paternal_count"))
  x<-copy(counts)
  x[,depth:=maternal_count+paternal_count]
  x<-x[is.finite(depth)&depth>=settings$min_depth&!is.na(ref_maternal)]
  setorder(x,gene_id,family_id,-depth,sample_id,locus_id)
  x<-x[,.SD[1],by=.(gene_id,family_id,locus_id)]
  scored<-x[,.(n_families=uniqueN(family_id),min_depth=min(depth)),by=.(gene_id,locus_id)]
  setorder(scored,gene_id,-n_families,-min_depth,locus_id)
  sel<-scored[,.SD[1],by=gene_id]
  x<-merge(x,sel[,.(gene_id,locus_id)],by=c("gene_id","locus_id"))
  setorder(x,gene_id,family_id,-depth,sample_id)
  x<-x[,.SD[1],by=.(gene_id,family_id)]
  cov<-fread(cfg$paths$covariates)
  assert_unique(cov,"sample_id","GUSTO covariates")
  x<-merge(x,cov,by="sample_id",all.x=TRUE)
  rows<-lapply(split(x,x$gene_id),function(d) {
    base<-data.table(gene_id=d$gene_id[1],locus_id=d$locus_id[1],
      n_families=nrow(d),n_ref_maternal=sum(d$ref_maternal),
      n_alt_maternal=sum(!d$ref_maternal),status="insufficient_families")
    if(nrow(d)<cfg$models$min_families)return(base)
    if(min(base$n_ref_maternal,base$n_alt_maternal)<cfg$models$min_orientation) {
      base$status<-"insufficient_cis_orientation";return(base)
    }
    names<-unlist(cfg$models$covariates)
    require_columns(d,names,"GUSTO gene covariates")
    z<-data.frame(cis=ifelse(d$ref_maternal,1,-1))
    if(length(names))for(nm in names)z[[nm]]<-d[[nm]]
    if(anyNA(z)){base$status<-"missing_covariates";return(base)}
    for(nm in names)if(is.numeric(z[[nm]]))z[[nm]]<-z[[nm]]-mean(z[[nm]])
    X<-tryCatch(model.matrix(~.,z),error=function(e)NULL)
    if(is.null(X)){base$status<-"invalid_covariates";return(base)}
    other<-setdiff(colnames(X),c("(Intercept)","cis"))
    if(length(other))X[,other]<-scale(X[,other,drop=FALSE],center=TRUE,scale=FALSE)
    y<-log((d$maternal_count+.5)/(d$paternal_count+.5))
    fit<-hc3_intercept(y,X)
    if(is.null(fit)){base$status<-"unidentifiable_or_leverage";return(base)}
    set.seed(cfg$models$seed)
    boots<-replicate(as.integer(cfg$models$bootstrap),{
      ii<-sample.int(nrow(d),nrow(d),replace=TRUE)
      if(qr(X[ii,,drop=FALSE])$rank!=ncol(X))NA_real_ else
        unname(lm.fit(X[ii,,drop=FALSE],y[ii])$coefficients[1])
    })
    fraction<-mean(is.finite(boots))
    interval<-if(fraction>=.9)quantile(boots,c(.025,.975),na.rm=TRUE) else c(NA_real_,NA_real_)
    base[,`:=`(estimate=fit$estimate,se=fit$se,p_value=fit$p,
      ci_low=unname(interval[1]),ci_high=unname(interval[2]),
      bootstrap_valid_fraction=fraction,
      residual_df=fit$df,status=if(fraction>=.9)"eligible_conditional_inference" else "unstable_bootstrap")]
    base
  })
  result<-rbindlist(rows,fill=TRUE)
  if(!ncol(result))result<-data.table(gene_id=character(),status=character(),p_value=numeric())
  if(!"p_value"%in%names(result))result[,p_value:=NA_real_]
  result[,`:=`(q_value=p.adjust(p_value,"BH"),effect_type="gene_parent_of_origin",
    evidence="expected_counts_per_trio_frozen_locus",
    hypothesis="log_maternal_paternal_ratio_adjusted_for_cis_and_covariates")]
  additional_inference_labels(result,"ln")
  result
}

# Human paired-event, same-locus additional analysis. All files stay private.
# Called by worker_additional for the human-trio profile.
additional_human_pairs <- function(event_counts,gene_counts,events,cfg,settings) {
  template<-data.table(gene_id=character(),gene_name=character(),event_a=character(),
    event_b=character(),effect_type=character(),status=character(),n_families=integer(),
    n_units=integer(),p_value=numeric(),q_value=numeric())
  pairs<-additional_event_pairs(events,settings)
  if(!nrow(pairs))return(template)
  require_columns(event_counts,c("sample_id","family_id","event_id","locus_id",
    "maternal_count","paternal_count"))
  require_columns(gene_counts,c("sample_id","gene_id","locus_id","ref_maternal"))
  evmap<-unique(events[,.(event_id,gene_id)])
  assert_unique(evmap,"event_id","Human structural event-to-gene map")
  e<-copy(event_counts)
  if(!"gene_id" %in% names(e))e<-merge(e,evmap,by="event_id",all.x=TRUE)
  e[,depth:=maternal_count+paternal_count]
  e<-e[!is.na(gene_id)&is.finite(depth)&depth>=settings$min_depth]
  assert_unique(e,c("sample_id","event_id","locus_id"),"Human event/sample/locus counts")
  # Orientation comes from the already completed, parent-assigned gene counts,
  # not from the observed event imbalance that is being tested.
  ori<-copy(gene_counts)[,.(ref_maternal={
    v<-unique(na.omit(ref_maternal))
    if(length(v)==1L)v[1] else NA
  }),by=.(sample_id,gene_id,locus_id)]
  assert_unique(ori,c("sample_id","gene_id","locus_id"),"Human gene cis orientation")
  e<-merge(e,ori,by=c("sample_id","gene_id","locus_id"),all.x=TRUE,
    suffixes=c("","_gene"))
  if("ref_maternal_gene"%in%names(e)) {
    e[,ref_maternal:=ref_maternal_gene]
    e[,ref_maternal_gene:=NULL]
  }
  cov<-fread(cfg$paths$covariates)
  assert_unique(cov,"sample_id","GUSTO covariates")
  if("family_id"%in%names(cov))cov[,family_id:=NULL]
  pairs<-pairs[order(gene_id,effect_type,-event_a_distinguishing_bp,
    -event_b_distinguishing_bp,event_a,event_b)]
  pairs<-pairs[,.SD[seq_len(min(.N,settings$max_pairs_per_gene*5L))],
    by=.(gene_id,effect_type)]
  by_event<-split(e,e$event_id)
  rows<-lapply(seq_len(nrow(pairs)),function(i) {
    p<-copy(pairs[i])
    result<-copy(p)[,`:=`(status="no_shared_callable_locus",n_families=0L,
      n_units=0L,locus_id=NA_character_,n_ref_maternal=0L,
      n_alt_maternal=0L,direct_coordinate_both=FALSE)]
    a<-by_event[[p$event_a]];b<-by_event[[p$event_b]]
    if(is.null(a)||is.null(b))return(result)
    d<-merge(a,b,by=c("sample_id","family_id","gene_id","locus_id"),
      suffixes=c("_a","_b"),sort=FALSE)
    if(!nrow(d))return(result)
    # Retain one locus selected only by family callability/depth, never effect.
    scored<-d[,.(n_families=uniqueN(family_id),
      min_depth=min(pmin(depth_a,depth_b))),by=locus_id]
    setorder(scored,-n_families,-min_depth,locus_id)
    d<-d[locus_id==scored$locus_id[1]]
    setorder(d,family_id,-depth_a,-depth_b,sample_id)
    d<-d[,.SD[1],by=family_id]  # no multiple children from one family
    n_ref<-sum(d$ref_maternal_a %in% TRUE)
    n_alt<-sum(d$ref_maternal_a %in% FALSE)
    coord<-if(all(c("distinguishing_coordinate_a","distinguishing_coordinate_b")%in%names(d)))
      all(d$distinguishing_coordinate_a%in%TRUE &
          d$distinguishing_coordinate_b%in%TRUE) else FALSE
    result[,`:=`(locus_id=scored$locus_id[1],n_families=nrow(d),
      n_units=nrow(d),n_ref_maternal=n_ref,n_alt_maternal=n_alt,
      direct_coordinate_both=coord,status="insufficient_families")]
    if(nrow(d)<cfg$models$min_families)return(result)
    if(anyNA(d$ref_maternal_a)||any(d$ref_maternal_a!=d$ref_maternal_b)) {
      result$status<-"missing_or_discordant_cis_orientation";return(result)
    }
    if(min(n_ref,n_alt)<cfg$models$min_orientation) {
      result$status<-"insufficient_cis_orientation";return(result)
    }
    d<-merge(d,cov,by="sample_id",all.x=TRUE)
    names<-unlist(cfg$models$covariates)
    require_columns(d,names,"GUSTO structural covariates")
    z<-data.frame(cis=ifelse(d$ref_maternal_a,1,-1))
    if(length(names))for(nm in names)z[[nm]]<-d[[nm]]
    if(anyNA(z)){result$status<-"missing_covariates";return(result)}
    for(nm in names)if(is.numeric(z[[nm]]))z[[nm]]<-z[[nm]]-mean(z[[nm]])
    X<-tryCatch(model.matrix(~.,z),error=function(e)NULL)
    if(is.null(X)||qr(X)$rank!=ncol(X)) {
      result$status<-"rank_deficient_covariates";return(result)
    }
    y<-log((d$maternal_count_a+.5)/(d$paternal_count_a+.5))-
      log((d$maternal_count_b+.5)/(d$paternal_count_b+.5))
    if(any(!is.finite(y))){result$status<-"invalid_ratio";return(result)}
    fit<-hc3_intercept(y,X)
    if(is.null(fit)){result$status<-"unidentifiable_or_leverage";return(result)}
    # Seed deterministic per pair, independent of parallel job ordering.
    seed_offset<-sum(utf8ToInt(paste(p$event_a,p$event_b,sep="/")))
    set.seed(as.integer((as.numeric(cfg$models$seed)+seed_offset)%%.Machine$integer.max))
    boots<-replicate(as.integer(cfg$models$bootstrap),{
      ii<-sample.int(nrow(d),nrow(d),replace=TRUE)
      if(qr(X[ii,,drop=FALSE])$rank!=ncol(X))NA_real_ else
        unname(lm.fit(X[ii,,drop=FALSE],y[ii])$coefficients[1])
    })
    fraction<-mean(is.finite(boots))
    interval<-if(fraction>=.9)quantile(boots,c(.025,.975),na.rm=TRUE) else c(NA_real_,NA_real_)
    result[,`:=`(estimate=fit$estimate,se=fit$se,p_value=fit$p,
      ci_low=unname(interval[1]),ci_high=unname(interval[2]),
      bootstrap_valid_fraction=fraction,residual_df=fit$df,
      status=if(fraction>=.9)"eligible_conditional_inference" else "unstable_bootstrap")]
    result
  })
  out<-rbindlist(rows,fill=TRUE)
  setorder(out,gene_id,effect_type,-n_families,event_a,event_b)
  out<-out[,.SD[seq_len(min(.N,settings$max_pairs_per_gene))],
    by=.(gene_id,effect_type)]
  if(!"p_value"%in%names(out))out[,p_value:=NA_real_]
  out[,q_value:=p.adjust(p_value,"BH")]
  additional_inference_labels(out,"ln")
  out[,`:=`(evidence="Ornaments_expected_counts_same_frozen_locus_one_family_per_observation",
    hypothesis="human_event_ratio_difference_adjusted_for_cis_and_covariates")]
  out
}

additional_te_annotation <- function(cfg,settings) {
  tx<-fread(cfg$paths$tx_map)
  m<-fread(cfg$paths$event_map)[window==50,.(transcript_id,event_id)]
  if(is.null(settings$te_rmsk))return(data.table(event_id=unique(m$event_id),
    tss_te_label="not_assessed",annotation_status="missing_build_matched_repeat_annotation"))
  # Genome-build match is an operator assertion; see runbook. RepeatMasker
  # coordinates must use this cohort's reference build before enabling.
  t<-annotate_repeats(tx,settings$te_rmsk)
  t<-merge(t,m,by="transcript_id")
  t[,.(tss_te_label=if(any(tss_unknown_repeat))"unknown_repeat" else
    if(all(tss_te))"TE" else if(!any(tss_te))"nonTE" else "mixed",
    n_transcripts=.N,te_class=paste(sort(unique(te_class[!is.na(te_class)])),collapse=";"),
    annotation_status="repeatmasker_reference_build_operator_verified"),by=event_id]
}


worker_additional <- function(cfg,root) {
  settings<-additional_settings(cfg);additional_validate(cfg,settings)
  reference_gate(cfg,root)
  if(!stage_valid(cfg$paths$model_receipt,
      unlist(cfg$paths[c("model_gene","model_promoter","model_pairs")])))
    stop("Complete frozen primary model outputs required before additional analysis")
  m<-load_manifests(cfg)
  if(cfg$profile=="public_cross"&&uniqueN(m$units$stage)>1L&&
     !fish_forced_five(cfg,m$units))
    stop("Mixed-stage fish cohort requires explicit five-F1 forced pooling policy")
  counts<-completed_counts(cfg,root,m)
  paths<-additional_paths(cfg)
  receipt<-additional_receipt(cfg)
  # Versioned output path: never replace or invalidate the primary receipts.
  if(stage_valid(receipt,paths))stop("additional_v1 already complete; no implicit overwrite")
  if(file.exists(receipt)||any(file.exists(paths)))
    stop("Partial/additional_v1 outputs exist: inspect before an explicit new analysis version")
  stage_begin(receipt)
  if(cfg$profile=="public_cross") {
    genes<-additional_public_genes(counts$gene,m$units,cfg,settings)
    pairs<-additional_public_pairs(counts$event,fread(cfg$paths$events),m$units,cfg,settings)
  } else {
    genes<-additional_human_genes(counts$gene,cfg,settings)
    pairs<-additional_human_pairs(counts$event,counts$gene,
      fread(cfg$paths$events),cfg,settings)
  }
  labels<-additional_te_annotation(cfg,settings)
  write_tsv(genes,paths[["genes"]]);write_tsv(pairs,paths[["structural"]])
  write_tsv(labels,paths[["annotation"]])
  write_tsv(data.table(cohort_id=cfg$cohort_id,species=cfg$species,
    stage=if(uniqueN(m$units$stage)>1L)"mixed_stage40_stage45_unadjusted" else m$units$stage[1],
    tissue=m$units$tissue[1],
    marker_proxy=isTRUE(cfg$markers$proxy),
    min_depth=settings$min_depth,min_per_direction=settings$min_per_direction,
    max_pairs_per_gene=settings$max_pairs_per_gene,
    source="frozen_completed_primary_counts",
    effect_scale=if(cfg$profile=="public_cross")"log2_ratio" else "natural_log_ratio",
    interpretation=if(fish_forced_five(cfg,m$units))
      "exploratory_stage_unadjusted_cluster_independence_unverified_small_N" else
      "exploratory_conditioned_on_quantifier_and_frozen_locus"),
    paths[["provenance"]])
  stage_finish(receipt,paths)
  cat("Additional read-only analysis complete: ",additional_dir(cfg),"\n",sep="")
}
