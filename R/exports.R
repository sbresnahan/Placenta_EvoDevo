export_columns <- c("gene_id","gene_name","event_a","event_b","pair_type","n_families",
 "status","estimate","se","p_value","q_value","ci_low","ci_high","residual_df",
 "bootstrap_valid_fraction","cis_basis","uncertainty","discordant")
validate_aggregate <- function(x, public_pairs, min_cell) {
  if (!is.numeric(min_cell) || length(min_cell)!=1 || is.na(min_cell) ||
      min_cell<1 || min_cell!=floor(min_cell)) stop("exports.min_cell must be a positive integer")
  if (any(!names(x) %in% export_columns)) stop("Non-allowlisted export columns")
  require_columns(x,export_columns,"Complete aggregate schema")
  assert_unique(public_pairs,"gene_id","Frozen public pairs")
  assert_unique(x,"gene_id")
  known <- merge(x,public_pairs[,.(gene_id,event_a,event_b,gene_name,pair_type)],
    by=c("gene_id","event_a","event_b"),all.x=TRUE,suffixes=c("","_public"))
  if (nrow(known)!=nrow(x) || anyNA(known$pair_type_public) ||
      any(known$gene_name!=known$gene_name_public) ||
      any(known$pair_type!=known$pair_type_public)) stop("Not public frozen annotation identities")
  if (anyNA(x$n_families) || any(x$n_families<min_cell) ||
      any(x$n_families!=floor(x$n_families))) stop("Cell suppression not satisfied")
  if (any(x$status!="eligible_conditional_inference")) stop("Export requires eligible_conditional_inference status")
  enums <- list(cis_basis=c("one_redundant_cis_term_removed","two_cis_terms"),
    uncertainty="between_families_conditional_on_Ornaments")
  for (nm in intersect(names(enums),names(x)))
    if (anyNA(x[[nm]]) || any(!x[[nm]] %in% enums[[nm]])) stop("Unexpected free text in export")
  numeric_columns<-intersect(names(x),c("n_families","estimate","se","p_value","q_value",
    "ci_low","ci_high","residual_df","bootstrap_valid_fraction"))
  for(nm in numeric_columns)
    if(!is.numeric(x[[nm]]) || any(!is.finite(x[[nm]]))) stop("Invalid numeric field: ",nm)
  for(nm in intersect(names(x),c("p_value","q_value","bootstrap_valid_fraction")))
    if(any(x[[nm]]<0 | x[[nm]]>1)) stop("Probability outside [0,1]")
  invisible(TRUE)
}
