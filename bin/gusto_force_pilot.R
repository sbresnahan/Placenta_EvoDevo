#!/usr/bin/env Rscript
# Opt-in GUSTO pilot recovery. Retain only samples with a valid pilot receipt
# AND a parseable, matching pilot report. Never synthesize success receipts.
# The regular collect-pilot worker subsequently creates the cohort report/receipt.
args <- commandArgs(trailingOnly=TRUE)
usage <- paste(
  'Usage: bin/Rscript_sif bin/gusto_force_pilot.R CONFIG [--preview | --force | --verify] [--max-excluded N]',
  'Example: bin/Rscript_sif bin/gusto_force_pilot.R "$CONFIG" --force', sep='\n')
if (!length(args)) stop(usage)
config_file <- args[1]
opts <- args[-1]
if (length(opts) == 0L) opts <- '--preview'
if (anyDuplicated(opts[startsWith(opts, '--')])) stop('Duplicate flag\n', usage)
if (!all(opts[startsWith(opts, '--')] %in% c('--preview','--force','--verify','--max-excluded')))
  stop('Unknown option\n', usage)
remaining <- opts
if ('--max-excluded' %in% remaining) {
  j <- match('--max-excluded', remaining)
  if (j == length(remaining)) stop('Missing --max-excluded value')
  remaining <- remaining[-c(j,j+1L)]
}
if (!all(remaining %in% c('--preview','--force','--verify'))) stop('Unexpected positional argument')
if (sum(c('--preview','--force','--verify') %in% opts)>1L) stop('Choose exactly one operation')
pos <- match('--max-excluded', opts)
max_excluded <- if (is.na(pos)) 5L else {
  if (pos == length(opts)) stop('Missing --max-excluded value')
  as.integer(opts[pos+1L])
}
if (is.na(max_excluded) || max_excluded < 0L || max_excluded > 50L)
  stop('--max-excluded must be in 0..50')
if (!is.na(pos) && (pos+2L <= length(opts)) && !startsWith(opts[pos+2L],'--'))
  stop('Unexpected argument after --max-excluded')

script_arg <- grep('^--file=',commandArgs(),value=TRUE)
if (!length(script_arg)) stop('Run this file through Rscript')
root <- normalizePath(file.path(dirname(sub('^--file=','',script_arg[1])),'..'),mustWork=TRUE)
source(file.path(root,'R/library_paths.R'))
source(file.path(root,'R/config.R'))
cfg <- read_config_v2(config_file)
if (dir.exists(cfg$tools$r_library)) .libPaths(cfg$tools$r_library,include.site=FALSE)
suppressPackageStartupMessages(library(data.table))
for (mod in c('common','workflow','parental','phase_origin','models','exports',
              'acquisition','public','reference','pipeline','stages','additional',
              'scheduler','seadragon','lsf')) {
  f <- file.path(root,'R',paste0(mod,'.R'))
  if (file.exists(f)) source(f)
}
if (!identical(cfg$cohort_id,'GUSTO') || !identical(cfg$profile,'human_trio'))
  stop('This recovery command is restricted to the GUSTO human-trio cohort')
required <- c('sources','libraries','units','trios','covariates','pilot_report','pilot_receipt')
if (!all(required %in% names(cfg$paths)))
  stop('Missing config paths: ', paste(setdiff(required,names(cfg$paths)),collapse=', '))
libs <- fread(cfg$paths$libraries, colClasses='character')
if (!all(c('sample_id','child_id','unit_id') %in% names(libs)))
  stop('GUSTO libraries manifest must contain sample_id, child_id, unit_id')
if (nrow(libs) < 10L || anyDuplicated(libs$sample_id) || anyNA(libs$sample_id))
  stop('Expected distinct, nonmissing sample IDs for GUSTO')

validate_pilot <- function(sample) {
  rec <- cfg_path(cfg,'library_pilot_receipt',sample=sample)
  rep <- cfg_path(cfg,'library_pilot_report',sample=sample)
  if (!file.exists(rec)) return('missing_pilot_receipt')
  if (!file.exists(rep) || file.info(rep)$size <= 0) return('missing_pilot_report')
  ok <- tryCatch(stage_valid(rec,rep), error=function(e) FALSE)
  if (!isTRUE(ok)) return('invalid_pilot_receipt')
  d <- tryCatch(fread(rep),error=function(e) NULL)
  if (is.null(d) || !nrow(d) || !('sample_id' %in% names(d)))
    return('invalid_pilot_report')
  if (anyNA(d$sample_id) || !all(as.character(d$sample_id)==sample))
    return('invalid_pilot_report')
  if ('mode' %in% names(d) && !any(as.character(d$mode)=='pilot'))
    return('invalid_pilot_report')
  'retained'
}
state <- data.table(sample_id=libs$sample_id)
state[, status:=vapply(sample_id,validate_pilot,character(1))]
keep <- state[status=='retained',sample_id]
excluded <- state[status!='retained']
cat(sprintf('GUSTO pilot receipts: %d valid, %d excluded (out of %d)\n',
            length(keep),nrow(excluded),nrow(libs)))
if (nrow(excluded)) print(excluded)
if (!length(keep) || nrow(excluded)>max_excluded || length(keep)<ceiling(0.9*nrow(libs)))
  stop('Refusing partial pilot recovery: too many invalid/missing pilots')
if ('--verify' %in% opts) {
  if (nrow(excluded)) stop('Filtered cohort still contains invalid pilot samples')
  cat('VERIFY PASS: all active GUSTO libraries have valid pilot receipts.\n')
  quit(save='no',status=0L)
}
if (!nrow(excluded)) stop('Nothing to recover: no invalid pilot samples')

# Construct all replacement tables in memory before modifying any working file.
# No source manifests are modified in the scripts repo; only the active cohort
# manifest copies are changed, with byte-identical backups available.
filtered <- list()
source_paths <- list()
for (key in c('sources','libraries','units','trios','covariates')) {
  path <- cfg$paths[[key]]
  if (is.null(path) || !file.exists(path)) stop('Missing manifest ',key,': ',path)
  if (grepl('/scripts/',normalizePath(path,mustWork=TRUE),fixed=TRUE))
    stop('Refusing to edit repository source manifest: ',path)
  x <- fread(path,colClasses='character')
  if (!nrow(x)) stop('Cannot filter empty manifest: ',path)
  column <- switch(key, sources='sample_id',libraries='sample_id',
                   units='unit_id',trios='child_id',covariates='sample_id')
  if (!column %in% names(x)) stop('Missing ',column,' in ',path)
  allowed <- switch(key,
    sources=keep, libraries=keep, units=libs[sample_id %in% keep,unit_id],
    trios=libs[sample_id %in% keep,child_id], covariates=keep)
  filtered[[key]] <- x[get(column) %in% allowed]
  if (!nrow(filtered[[key]])) stop('Filtering empties ',key,'; stopping')
  source_paths[[key]] <- path
  cat(sprintf('  %-12s %d -> %d rows\n',key,nrow(x),nrow(filtered[[key]])))
}
if (!setequal(filtered$libraries$sample_id,keep) ||
    !setequal(filtered$sources$sample_id,keep) ||
    !setequal(filtered$covariates$sample_id,keep) ||
    !setequal(filtered$trios$child_id,libs[sample_id %in% keep,child_id]))
  stop('Inconsistent filtered manifest keys; no files changed')
if (!'--force' %in% opts) {
  cat('PREVIEW ONLY. Add --force to archive and update working GUSTO manifests.\n')
  quit(save='no',status=0L)
}
if (stage_valid(cfg$paths$pilot_receipt,cfg$paths$pilot_report))
  stop('Pilot is already complete. Refusing to change manifests.')

# Keep backups separate from analysis receipts, in the same cohort tree.
base_dir <- dirname(dirname(cfg$paths$pilot_receipt))
archive <- file.path(base_dir,'private',paste0('force_pilot_',format(Sys.time(),'%Y%m%dT%H%M%S')))
if (dir.exists(archive)) stop('Backup directory already exists: ',archive)
dir.create(archive,recursive=TRUE,mode='0700')
if (!dir.exists(archive)) stop('Cannot create backup directory: ',archive)
manifest_archive <- file.path(archive,'manifests')
dir.create(manifest_archive,mode='0700')
for (key in names(source_paths)) {
  p <- source_paths[[key]]
  if (!file.copy(p,file.path(manifest_archive,paste0(key,'.tsv')),overwrite=FALSE))
    stop('Cannot back up: ',p)
}
write.table(state,file=file.path(archive,'pilot_sample_audit.tsv'),sep='\t',
            row.names=FALSE,col.names=TRUE,quote=TRUE,na='')
writeLines(c(paste0('cohort=GUSTO'),paste0('retained=',length(keep)),
             paste0('excluded=',nrow(excluded)),
             paste0('utc=',format(Sys.time(),tz='UTC',usetz=TRUE)),
             'status=pilot_manifest_filter_applied; pilot receipt NOT modified'),
           con=file.path(archive,'audit.txt'))
# Stage each updated file, then replace atomically. Roll back on error.
staged <- list()
for (key in names(source_paths)) {
  p <- source_paths[[key]]
  tmp <- paste0(p,'.force_pilot_new_',Sys.getpid())
  if (file.exists(tmp)) stop('Staging file already exists: ',tmp)
  fwrite(filtered[[key]],tmp,sep='\t',quote=TRUE,na='')
  staged[[key]] <- tmp
}
changed <- character()
tryCatch({
  for (key in names(source_paths)) {
    p <- source_paths[[key]]
    if (!file.rename(staged[[key]],p)) stop('Cannot replace ',p)
    changed <- c(changed,key)
  }
  # Validate using the same manifest parser used by the production pipeline.
  invisible(load_manifests(cfg))
},error=function(e) {
  for (key in changed) {
    file.copy(file.path(manifest_archive,paste0(key,'.tsv')),
              source_paths[[key]],overwrite=TRUE)
  }
  stop('Rolled back active manifests: ',conditionMessage(e))
},finally={
  for (tmp in staged) if(file.exists(tmp))unlink(tmp)
})
cat('Filtered GUSTO manifests committed. Backup: ',archive,'\n',sep='')
cat('Pilot completion receipt has NOT been fabricated.\n')
cat('Next: run the genuine collect-pilot worker on these retained samples.\n')
