#!/usr/bin/env Rscript
# Recover a partially completed GUSTO quantification run by retaining only
# biological units with complete, verified Ornaments/assignment results.
# No sample IDs are hard-coded; no model or completion receipt is fabricated.
args <- commandArgs(trailingOnly=TRUE)
usage <- paste(
  'Usage: bin/Rscript_sif bin/gusto_force_quantification.R CONFIG',
  '       [--preview | --force | --verify] [--max-excluded N] [--expect-excluded N]',
  sep='\n')
if (!length(args)) stop(usage)
config_file <- args[1L]
opts <- args[-1L]
mode <- '--preview'
operation_seen <- FALSE
max_excluded <- 10L
expected_excluded <- NA_integer_
i <- 1L
seen <- character()
while (i <= length(opts)) {
  flag <- opts[i]
  if (flag %in% seen) stop('Duplicate flag: ',flag)
  seen <- c(seen,flag)
  if (flag %in% c('--preview','--force','--verify')) {
    if (operation_seen) stop('Choose exactly one operation\n', usage)
    operation_seen <- TRUE
    mode <- flag
    i <- i+1L
  } else if (flag %in% c('--max-excluded','--expect-excluded')) {
    if (i==length(opts)) stop('Missing value for ',flag)
    value <- suppressWarnings(as.integer(opts[i+1L]))
    if (is.na(value) || !identical(as.character(value),opts[i+1L]) || value<0L)
      stop('Expected a nonnegative integer after ',flag)
    if (flag=='--max-excluded') max_excluded<-value else expected_excluded<-value
    i <- i+2L
  } else stop('Unknown argument: ',flag,'\n',usage)
}
if (max_excluded>100L) stop('max-excluded must be <=100; investigate mass failures')
if (!is.na(expected_excluded) && mode=='--verify' && expected_excluded!=0L)
  stop('--verify expects zero incomplete units in the active manifest')

filearg <- grep('^--file=',commandArgs(),value=TRUE)
if (!length(filearg)) stop('Run via Rscript, not source()')
root <- normalizePath(file.path(dirname(sub('^--file=','',filearg[1L])),'..'),mustWork=TRUE)
source(file.path(root,'R/library_paths.R'))
source(file.path(root,'R/config.R'))
cfg <- read_config_v2(config_file)
if (dir.exists(cfg$tools$r_library)) .libPaths(cfg$tools$r_library,include.site=FALSE)
suppressPackageStartupMessages(library(data.table))
for (module in c('common','workflow','parental','phase_origin','models','exports',
                 'acquisition','public','reference','pipeline','stages','additional',
                 'scheduler','seadragon','lsf')) {
  f <- file.path(root,'R',paste0(module,'.R'))
  if (file.exists(f)) source(f)
}
if (!identical(cfg$cohort_id,'GUSTO') || !identical(cfg$profile,'human_trio'))
  stop('This recovery implementation is restricted to the GUSTO human-trio profile')
keys <- c('sources','libraries','units','trios','covariates')
# Only manifest paths are direct cfg$paths entries. Quantification, assignment,
# and count paths are per-library templates resolved by cfg_path() below.
required <- keys
if (!all(required %in% names(cfg$paths)))
  stop('Missing config paths: ',paste(setdiff(required,names(cfg$paths)),collapse=', '))

libs <- fread(cfg$paths$libraries,colClasses='character')
require_columns(libs,c('sample_id','unit_id','child_id','strand_mode'))
if (nrow(libs)<10L || anyDuplicated(libs$sample_id) ||
    anyNA(libs[,.(sample_id,unit_id,child_id,strand_mode)]) ||
    any(!nzchar(libs$sample_id) | !nzchar(libs$unit_id) |
        !nzchar(libs$child_id) | !nzchar(libs$strand_mode)))
  stop('Invalid GUSTO library manifest')
if (libs[,any(uniqueN(child_id)!=1L),by=unit_id][,any(V1)])
  stop('One biological unit maps to multiple children')
if (libs[,any(uniqueN(unit_id)!=1L),by=child_id][,any(V1)])
  stop('One child maps to multiple biological units')

# Validate each sample's quantification receipt and assigned count files.
# Always use the existing stage_valid() rather than inferring success from LSF.
receipt_state <- function(i) {
  s <- libs$sample_id[i]; strand <- libs$strand_mode[i]
  path <- function(key) cfg_path(cfg,key,sample=s,mode='complete',strand=strand)
  q<-path('quant_receipt'); a<-path('assign_receipt')
  genes<-path('gene_counts'); events<-path('event_counts')
  if (!file.exists(q) || is.na(file.info(q)$size) || file.info(q)$size<=0)
    return('missing_quant_receipt')
  qr <- tryCatch(readRDS(q),error=function(e) NULL)
  if (!is.list(qr) || is.null(qr$metrics)) return('invalid_quant_receipt')
  if (!file.exists(a)) return('missing_assign_receipt')
  outputs<-c(genes,events)
  if (any(!file.exists(outputs)) || any(is.na(file.info(outputs)$size)) ||
      any(file.info(outputs)$size<=0)) return('missing_assigned_counts')
  ok <- tryCatch(stage_valid(a,outputs),error=function(e) FALSE)
  if (!isTRUE(ok)) return('invalid_assign_receipt')
  # Read a header/first row to reject truncated or unparsable counts cheaply.
  readable <- vapply(outputs,function(p) {
    d <- tryCatch(suppressWarnings(fread(p,nrows=1L)),error=function(e) NULL)
    !is.null(d) && ncol(d)>0L
  },logical(1))
  if (!all(readable)) return('unreadable_assigned_counts')
  'complete'
}
libs[,receipt_state:=vapply(seq_len(nrow(libs)),receipt_state,character(1))]
# A biological unit is retained only if every library comprising it is complete.
units <- libs[,.(n_libraries=.N,
                status=if (all(receipt_state=='complete')) 'complete' else 'incomplete',
                reasons=paste(sort(unique(receipt_state[receipt_state!='complete'])),collapse=';')),
              by=.(unit_id,child_id)]
setorder(units,unit_id)
excluded <- units[status!='complete']
keep_units <- units[status=='complete',unit_id]
kept <- libs[unit_id %in% keep_units]
cat(sprintf('GUSTO scale eligibility: %d/%d complete biological units, %d excluded; %d/%d libraries retained\n',
            length(keep_units),nrow(units),nrow(excluded),nrow(kept),nrow(libs)))
if (nrow(excluded)) print(excluded)
if (!is.na(expected_excluded) && nrow(excluded)!=expected_excluded)
  stop('Expected ',expected_excluded,' excluded biological units but found ',nrow(excluded))
if (mode=='--verify') {
  if (nrow(excluded)) stop('VERIFY FAIL: incomplete units remain in the active manifest')
  cat('VERIFY PASS: all active units have verified quantification and assignment outputs.\n')
  quit(save='no',status=0L)
}
if (nrow(excluded)==0L) {
  if (mode=='--force') stop('No units need exclusion; use --verify --submit to rerun collector')
  cat('PREVIEW: nothing to exclude. If collector failed, run --verify --submit.\n')
  quit(save='no',status=0L)
}
if (nrow(excluded)>max_excluded || length(keep_units)<ceiling(.90*nrow(units)))
  stop('Exclusion exceeds guardrails: --max-excluded or <90% of units retained')

# Assemble every manifest in memory before touching any file.
# Every modified manifest must be in the same active cohort-local directory.
manifest_dir <- normalizePath(dirname(cfg$paths$libraries),mustWork=TRUE)
if (grepl('/scripts(/|$)',manifest_dir)) stop('Refusing to edit repository manifests')
cohort_dir <- normalizePath(file.path(manifest_dir,'..'),mustWork=TRUE)
if (basename(cohort_dir)!='GUSTO' || basename(manifest_dir)!='manifests')
  stop('Active GUSTO manifests must reside in data/GUSTO/manifests')
filtered<-list(); originals<-list(); paths<-list()
for (key in keys) {
  p<-normalizePath(cfg$paths[[key]],mustWork=TRUE)
  if (dirname(p)!=manifest_dir) stop('Manifest outside cohort directory: ',p)
  x<-fread(p,colClasses='character')
  field<-switch(key,sources='sample_id',libraries='sample_id',units='unit_id',
                trios='child_id',covariates='sample_id')
  if (!field %in% names(x)) stop('Missing ',field,' in ',p)
  allowed<-switch(key,sources=kept$sample_id,libraries=kept$sample_id,
                  units=kept$unit_id,trios=kept$child_id,covariates=kept$sample_id)
  out<-x[get(field) %in% allowed]
  if (!nrow(out)) stop('Filtered manifest empty: ',key)
  if (key %in% c('libraries','units','trios','covariates') &&
      !setequal(unique(out[[field]]),unique(allowed)))
    stop('Filtered manifest lacks expected ',field,' keys: ',key)
  if (key=='sources' && !setequal(unique(out$sample_id),unique(allowed)))
    stop('Filtered sources lack expected samples')
  filtered[[key]]<-out; originals[[key]]<-x; paths[[key]]<-p
  cat(sprintf('  %-12s %d -> %d rows\n',key,nrow(x),nrow(out)))
}
if (mode=='--preview') {
  cat('PREVIEW ONLY. No data, manifests, or receipts changed.\n')
  quit(save='no',status=0L)
}

# Persist reproducible provenance before replacing cohort-local manifests.
archive<-file.path(cohort_dir,'private',paste0('force_quant_',format(Sys.time(),'%Y%m%dT%H%M%S'),
                                                   '_',Sys.getpid()))
if (dir.exists(archive)) stop('Backup already exists: ',archive)
dir.create(archive,recursive=TRUE,mode='0700')
if (!dir.exists(archive)) stop('Cannot create backup directory: ',archive)
archive_manifests<-file.path(archive,'manifests')
dir.create(archive_manifests,mode='0700')
for (key in keys) {
  dest<-file.path(archive_manifests,paste0(key,'.tsv'))
  if (!file.copy(paths[[key]],dest,overwrite=FALSE))
    stop('Backup failed for ',key,'; active manifests unchanged')
}
fwrite(libs[,.(sample_id,child_id,unit_id,receipt_state,
               excluded=!(unit_id %in% keep_units))],
       file.path(archive,'library_audit.tsv'),sep='\t')
fwrite(units[,.(unit_id,child_id,n_libraries,status,reasons)],
       file.path(archive,'unit_audit.tsv'),sep='\t')
writeLines(c(paste0('cohort=',cfg$cohort_id),
             paste0('excluded_units=',nrow(excluded)),
             paste0('retained_units=',length(keep_units)),
             paste0('timestamp_utc=',format(Sys.time(),tz='UTC',usetz=TRUE)),
             'original_quantification_receipts_unchanged=true',
             'original_assignment_receipts_unchanged=true'),
           file.path(archive,'provenance.txt'))

# Stage updates, replace them, and restore all originals if validation fails.
staged<-list()
for (key in keys) {
  tmp<-paste0(paths[[key]],'.force_quant_tmp_',Sys.getpid())
  if (file.exists(tmp)) stop('Temporary manifest exists: ',tmp)
  fwrite(filtered[[key]],tmp,sep='\t',quote=TRUE,na='')
  staged[[key]]<-tmp
}
changed<-character()
tryCatch({
  for (key in keys) {
    if (!file.rename(staged[[key]],paths[[key]])) stop('Cannot replace ',key)
    changed<-c(changed,key)
  }
  invisible(load_manifests(cfg))
},error=function(e) {
  for (key in changed) {
    if (!file.copy(file.path(archive_manifests,paste0(key,'.tsv')),
                   paths[[key]],overwrite=TRUE))
      warning('Manual manifest restoration required: ',key)
  }
  stop('Manifest update failed and rollback attempted: ',conditionMessage(e))
},finally={
  for (tmp in staged) if (file.exists(tmp)) unlink(tmp)
})
cat('SUCCESS: active manifests filtered; original manifests saved in ',archive,'\n',sep='')
cat('No output or receipt was marked complete. Next: run actual collector.\n')
