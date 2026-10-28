orient_public <- function(x,markers,unit,c) {
  require_columns(markers,c("chrom","genomic_pos","ref","alt","focal_index","alternate_index"))
  assert_unique(markers,c("chrom","genomic_pos"),"Marker coordinates")
  y<-merge(x,markers,by=c("chrom","genomic_pos","ref","alt"),all.x=TRUE,sort=FALSE)
  good<-!is.na(y$focal_index)&!is.na(y$alternate_index)&
    y$focal_index %in% 0:1 & y$alternate_index %in% 0:1 & y$focal_index!=y$alternate_index
  maternal<-if(unit$maternal_parent==c$parents$focal)y$focal_index else y$alternate_index
  y[,ref_maternal:=ifelse(good,maternal==0L,NA)]
  y[,maternal_count:=ifelse(good,ifelse(ref_maternal,ref_count,alt_count),NA_real_)]
  y[,paternal_count:=ifelse(good,ifelse(ref_maternal,alt_count,ref_count),NA_real_)]
  y[,parent_status:=ifelse(good,"external_parent_markers","missing_marker")]
  y
}
public_fit <- function(y,orientation,independence=TRUE) {
  if(any(!is.finite(y))||any(!orientation %in% c(-1,1)))stop("Invalid reciprocal model inputs")
  X<-cbind(1,orientation);rank<-qr(X)$rank;df<-length(y)-rank
  z<-data.table(n_units=length(y),rank=rank,residual_df=df,estimate=NA_real_,
    se=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p_value=NA_real_,status="descriptive_only")
  if(!independence) {z$status<-"independence_unresolved";return(z)}
  if(length(unique(orientation))<2L) {z$status<-"single_direction_cis_POE_confounded";return(z)}
  if(rank<2L||df<=0) {z$status<-"nonpositive_df_or_rank_deficient";return(z)}
  fit<-lm.fit(X,y);v<-sum(fit$residuals^2)/df*solve(crossprod(X));se_value<-sqrt(v[1,1])
  z[,`:=`(estimate=unname(fit$coefficients[1]),se=se_value,status="eligible_conditional_inference")]
  z[,`:=`(ci_low=estimate-qt(.975,df)*se,ci_high=estimate+qt(.975,df)*se,
    p_value=if(se>0)2*pt(-abs(estimate/se),df) else NA_real_)]
  if(se_value==0)z$status<-"degenerate_residual_variance"
  z
}
freeze_sentinels <- function(counts,id,units,min_depth=20) {
  require_columns(counts,c("unit_id",id,"locus_id","maternal_count","paternal_count"))
  assert_unique(counts,c("unit_id",id,"locus_id"),"Unique unit/entity/locus allocations")
  x<-copy(counts);x[,depth:=maternal_count+paternal_count]
  x<-x[unit_id %in% units & is.finite(depth) & depth>=min_depth]
  score<-x[,.(n_units=uniqueN(unit_id),min_depth=min(depth)),by=c(id,"locus_id")]
  score<-score[n_units==length(units)]
  setorderv(score,c(id,"min_depth","locus_id"),c(1,-1,1))
  sel<-score[, .SD[1],by=id]
  merge(x,sel[,c(id,"locus_id"),with=FALSE],by=c(id,"locus_id"))
}
public_endpoints <- function(gene,event,units,c) {
  # No dependence-based pseudo-replication: unresolved/duplicate clusters remain descriptive.
  independent<-!anyNA(units$cluster_id)&&all(nzchar(units$cluster_id))&&
    !anyDuplicated(units$cluster_id)
  fit_entities<-function(x,id,promoter=FALSE) {
    if(!nrow(x))return(data.table(entity_id=character(),status=character(),p_value=numeric(),q_value=numeric()))
    ans<-rbindlist(lapply(split(x,x[[id]]),function(d) {
      setorder(d,unit_id);u<-units[match(d$unit_id,unit_id)]
      y<-if(promoter)d$contrast else log2((d$maternal_count+.5)/(d$paternal_count+.5))-
        log2(c$dosage$maternal/c$dosage$paternal)
      z<-public_fit(y,ifelse(u$maternal_parent==c$parents$focal,1,-1),independent)
      z[,entity_id:=d[[id]][1]];z
    }),fill=TRUE)
    ans[,q_value:=p.adjust(p_value,"BH")];ans
  }
  g<-freeze_sentinels(gene,"gene_id",units$unit_id,c$models$min_depth)
  e<-freeze_sentinels(event,"event_id",units$unit_id,c$models$min_depth)
  pairs<-list()
  for(gid in sort(unique(e$gene_id))) {
    q<-e[gene_id==gid];ids<-sort(unique(q$event_id))
    if(length(ids)<2L)next
    candidates<-combn(ids,2,simplify=FALSE)
    scored<-rbindlist(lapply(candidates,function(p) {
      a<-q[event_id==p[1]];b<-q[event_id==p[2]]
      ab<-merge(a,b,by="unit_id",suffixes=c("_a","_b"))
      data.table(event_a=p[1],event_b=p[2],n_units=nrow(ab),min_depth=min(ab$depth_a,ab$depth_b))
    }))
    setorder(scored,-n_units,-min_depth,event_a,event_b);s<-scored[1]
    d<-merge(q[event_id==s$event_a],q[event_id==s$event_b],by="unit_id",suffixes=c("_a","_b"))
    d[,`:=`(gene_id=gid,event_a=s$event_a,event_b=s$event_b,
      contrast=log2((maternal_count_a+.5)/(paternal_count_a+.5))-
        log2((maternal_count_b+.5)/(paternal_count_b+.5)))]
    pairs[[length(pairs)+1L]]<-d
  }
  p<-rbindlist(pairs,fill=TRUE)
  if(!ncol(p))p<-data.table(gene_id=character(),unit_id=character())
  complete_callability<-function(result,universe,reason) {
    all<-data.table(entity_id=sort(unique(universe)))
    z<-merge(all,result,by="entity_id",all.x=TRUE)
    z[is.na(status),status:=reason]
    z
  }
  list(gene=complete_callability(fit_entities(g,"gene_id"),gene$gene_id,
      "no_frozen_sentinel_depth_every_unit"),
    promoter=complete_callability(fit_entities(p,"gene_id",TRUE),gene$gene_id,
      "fewer_than_two_fully_callable_groups"),paired=p)
}
