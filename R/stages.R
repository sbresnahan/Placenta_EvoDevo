reference_outputs <- function(c)c(unlist(c$paths[c("genome_fai","transcripts_fai","gtf","transcripts","exons","tx_map","event_map",
  "events","pairs","annotation_audit","deferred")]),if(c$profile=="public_cross")
    unlist(c$paths[c("markers","tx_vcf","tx_marker_map","personal_fasta","index")]))
reference_gate <- function(c,root) {
  if(!stage_valid(c$paths$reference_receipt,reference_outputs(c)))
    stop("Missing/incomplete reference receipt or outputs; run reference")
}
worker_reference <- function(c,root) {
  tool_gate(c);snapshot(c)
  stage_begin(c$paths$reference_receipt)
  for(p in reference_outputs(c))ensure_parent(p)
  build_reference(c);ensure_parent(c$paths$reference_receipt)
  stage_finish(c$paths$reference_receipt,reference_outputs(c))
}
worker_pilot <- function(c,root) {
  tool_gate(c);snapshot(c);reference_gate(c,root);m<-load_manifests(c)
  stage_begin(c$paths$pilot_receipt)
  if(c$profile=="human_trio")human_linkage(c,m)
  metrics<-list()
  for(sample in m$libraries$sample_id) {
    l<-m$libraries[sample_id==sample]
    diag<-rbindlist(lapply(c("fr","rf","unstranded"),function(strand)
      quant_sample(c,root,m,sample,"strand",strand,c$pilot$strand_fragments)))
    # Stranded choice must agree with the bounded diagnostic and local protocol review.
    selected<-diag[strand==l$strand_mode]$pseudoalignment_fraction
    if(!all(l$strand_verified %in% TRUE)||selected<c$pilot$min_pseudoalignment_fraction ||
       (l$strand_mode!="unstranded"&&
        selected<max(diag[strand!="unstranded"]$pseudoalignment_fraction)-c$pilot$strand_tolerance))
      stop("Strand/protocol/diagnostic approval failed: ",sample)
    q<-quant_sample(c,root,m,sample,"pilot",limit=c$pilot$fragments)
    if(q$fragments<c$pilot$fragments)stop("Pilot smaller than required throughput chunk")
    assign_sample(c,m,sample,"pilot")
    metrics[[length(metrics)+1L]]<-rbind(diag,q)
  }
  report<-rbindlist(metrics);nreads<-sum(m$sources$read_count,na.rm=TRUE)
  # ENA read_count counts reads, not paired fragments.
  denom<-ifelse(m$sources$layout=="PE",2,1)
  total<-sum(m$sources$read_count/denom)
  if(!is.finite(total)||total<=0)stop("Complete-cohort read counts required for projection")
  throughput<-report[mode=="pilot",sum(fragments)/sum(elapsed_seconds)]
  report[,`:=`(projected_quant_seconds=total/throughput,configured_total_fragments=total,
    prefix_sampling="nonrandom_prefix_not_whole_file_checksum_identity",
    actual_focal_parent=if(c$profile=="public_cross")c$parents$focal else "typed_trio",
    actual_alternate_parent=if(c$profile=="public_cross")c$parents$alternate else "typed_trio")]
  if(any(report$projected_quant_seconds>c$pilot$max_projected_seconds)||
      any(report$work_bytes_endpoint_lower_bound>c$pilot$max_work_bytes))
    stop("Measured projection exceeds configured local budget; no scale authorization")
  write_tsv(report,c$paths$pilot_report);ensure_parent(c$paths$pilot_receipt)
  stage_finish(c$paths$pilot_receipt,c$paths$pilot_report)
  cat("Technical pilot complete. Review the report before explicitly submitting scale.\n")
}
authorize_scale <- function(c,root) {
  reference_gate(c,root)
  if(!stage_valid(c$paths$pilot_receipt,c$paths$pilot_report))stop("Completed pilot required")
  if(!isTRUE(c$lsf$resource_limits_approved))stop("Local projected CPU/storage/budget acceptance required")
  ensure_parent(c$paths$scale_receipt)
  saveRDS(list(authorized=TRUE,approved_at=format(Sys.time(),tz="UTC",usetz=TRUE),
    approval="explicit_scale_command_and_local_resource_attestation"),c$paths$scale_receipt)
}
scale_gate <- function(c,root,m) {
  if(!file.exists(c$paths$scale_receipt))stop("Explicit scale authorization missing")
  r<-tryCatch(readRDS(c$paths$scale_receipt),error=function(e)NULL)
  if(!is.list(r)||!isTRUE(r$authorized))stop("Explicit scale authorization missing or invalid")
  if(!stage_valid(c$paths$pilot_receipt,c$paths$pilot_report))stop("Completed pilot required")
}
assignment_outputs <- function(c,sample,mode,l) {
  vapply(c("gene_counts","event_counts","assignment_qc"),function(k)
    cfg_path(c,k,sample=sample,mode=mode,strand=l$strand_mode),character(1))
}
worker_scale <- function(c,root) {
  tool_gate(c);snapshot(c);reference_gate(c,root);m<-load_manifests(c);scale_gate(c,root,m)
  allmetrics<-list()
  for(sample in m$libraries$sample_id) {
    l<-m$libraries[sample_id==sample]
    allmetrics[[length(allmetrics)+1L]]<-quant_sample(c,root,m,sample)
    qr<-cfg_path(c,"quant_receipt",sample=sample,mode="complete",strand=l$strand_mode)
    ar<-cfg_path(c,"assign_receipt",sample=sample,mode="complete",strand=l$strand_mode)
    outs<-assignment_outputs(c,sample,"complete",l)
    stage_begin(ar)
    assign_sample(c,m,sample);ensure_parent(ar);stage_finish(ar,outs)
  }
  write_tsv(rbindlist(allmetrics),c$paths$metrics)
  cat("Complete configured cohort processed. Models require a separate command.\n")
}
completed_counts <- function(c,root,m) {
  ans<-list(gene=list(),event=list())
  for(sample in m$libraries$sample_id) {
    if(identical(c$execution$backend,"seadragon")&&c$profile=="human_trio")child_index_gate(c,root,sample)
    l<-m$libraries[sample_id==sample]
    qr<-cfg_path(c,"quant_receipt",sample=sample,mode="complete",strand=l$strand_mode)
    ar<-cfg_path(c,"assign_receipt",sample=sample,mode="complete",strand=l$strand_mode)
    outs<-assignment_outputs(c,sample,"complete",l)
    if(!stage_valid(ar,outs))stop("Incomplete primary biological unit: ",sample)
    quant_outputs<-vapply(c("allele_counts","thetas","tpms","quant_log","resource_log"),
      function(k)cfg_path(c,k,sample=sample,mode="complete",strand=l$strand_mode),character(1))
    if(!stage_valid(qr,quant_outputs))stop("Quantification incomplete: ",sample)
    ans$gene[[sample]]<-fread(outs[1]);ans$event[[sample]]<-fread(outs[2])
  }
  lapply(ans,rbindlist,fill=TRUE)
}
worker_models <- function(c,root) {
  tool_gate(c);snapshot(c);reference_gate(c,root);m<-load_manifests(c);scale_gate(c,root,m)
  stage_begin(c$paths$model_receipt)
  x<-completed_counts(c,root,m)
  if(c$profile=="public_cross") {
    distinguishable<-fread(c$paths$events)[window==50 & distinguishing_bp>0,event_id]
    x$event<-x$event[event_id %in% distinguishable]
    out<-public_endpoints(x$gene,x$event,m$units,c)
    for(k in c("gene","promoter")) {
      out[[k]][,`:=`(cohort_id=c$cohort_id,stage=m$units$stage[1],species=c$species,
        marker_proxy=isTRUE(c$markers$proxy),proxy_reason=if(isTRUE(c$markers$proxy))c$markers$proxy_reason else "",
        focal_parent=c$parents$focal,alternate_parent=c$parents$alternate,
        fragment_assumptions=if(any(m$libraries$layout=="SE"))"configured_SE_assumptions_not_grid_tested" else "paired_end",
        fragment_mean_nt=if(any(m$libraries$layout=="SE"))c$fragments$fragment_mean_nt else NA_real_,
        fragment_sd_nt=if(any(m$libraries$layout=="SE"))c$fragments$fragment_sd_nt else NA_real_,
        fragment_parameter_origin=if(any(m$libraries$layout=="SE"))
          c$fragments$fragment_parameter_origin else "paired_end_quantifier",
        uncertainty="between_independent_units_conditional_on_Ornaments")]
      write_tsv(out[[k]],c$paths[[paste0("model_",k)]])
    }
    write_tsv(out$paired,c$paths$model_pairs)
  } else {
    if(!isTRUE(c$human$covariates_approved))stop("Covariate specification approval required")
    human_linkage(c,m)
    p<-choose_sentinels(x$event,fread(c$paths$pairs),c$models$min_depth)
    cov<-fread(c$paths$covariates);assert_unique(cov,"sample_id")
    if(any(!p$sample_id %in% cov$sample_id))stop("Missing exact covariate linkage")
    p<-merge(p,cov,by="sample_id",all.x=TRUE)
    out<-fit_all_pairs(p,unlist(c$models$covariates),c$models$bootstrap,c$models$min_families,
      c$models$min_orientation,c$models$seed)
    write_tsv(p,c$paths$model_pairs);write_tsv(out,c$paths$model_promoter)
    write_tsv(data.table(status="deferred_unimplemented"),c$paths$model_gene)
  }
  stage_finish(c$paths$model_receipt,
    unlist(c$paths[c("model_gene","model_promoter","model_pairs")]))
}
status_v2 <- function(c,root) {
  m<-load_manifests(c)
  actions<-list(
    reference=function()reference_gate(c,root),
    pilot=function() {
      if(!stage_valid(c$paths$pilot_receipt,c$paths$pilot_report))
        stop("Missing/incomplete pilot receipt or report")
    },
    scale=function()scale_gate(c,root,m),
    quantification=function(){scale_gate(c,root,m);completed_counts(c,root,m)},
    models=function() {
      if(!stage_valid(c$paths$model_receipt,
        unlist(c$paths[c("model_gene","model_promoter","model_pairs")])))stop("Missing/incomplete model receipt or outputs")
    })
  z<-rbindlist(lapply(names(actions),function(k)tryCatch({
    actions[[k]]();data.table(stage=k,state=if(k=="scale")"authorized" else "complete",reason="")
  },error=function(e)data.table(stage=k,state="blocked",reason=conditionMessage(e)))))
  write_tsv(z,c$paths$status);print(z)
}
export_v2 <- function(c) {
  if(c$profile!="human_trio")stop("Public outputs need no human-release bypass")
  e<-c$exports
  if(!isTRUE(e$enabled)||!isTRUE(e$policy_approved)||!isTRUE(e$complementary_suppression_reviewed)||
      is.null(e$review_token)||!nzchar(e$review_token))stop("Release disabled: exact manual policy review required")
  if(!file.exists(c$paths$candidate)||file.info(c$paths$candidate)$size<=0)
    stop("Reviewed candidate file missing or empty")
  x<-fread(c$paths$candidate);validate_aggregate(x,fread(c$paths$pairs),e$min_cell)
  write_tsv(x,c$paths$release)
}
review_v2 <- function(c) {
  if(c$profile!="human_trio"||!isTRUE(c$exports$policy_approved))stop("Locally approved human suppression policy required")
  positive(c$exports$min_cell,"approved min_cell")
  if(c$exports$min_cell!=floor(c$exports$min_cell))stop("Integer cell threshold required")
  x<-fread(c$paths$model_promoter)
  require_columns(x,export_columns)
  x<-x[status=="eligible_conditional_inference" & n_families>=c$exports$min_cell]
  x<-x[,..export_columns]
  validate_aggregate(x,fread(c$paths$pairs),c$exports$min_cell)
  write_tsv(x,c$paths$candidate)
  cat("Local review candidate only; inspect it before separately authorizing export. No transfer performed.\n")
}
