range_ok <- function(status,headers,offset,received,total) {
  if(offset==0 && status==200L)return(received==total)
  if(status!=206L)return(FALSE)
  h<-tolower(rawToChar(headers))
  m<-regmatches(h,regexec("content-range:[ \t]*bytes ([0-9]+)-([0-9]+)/([0-9]+)",h))[[1]]
  length(m)==4L && as.numeric(m[2])==offset &&
    as.numeric(m[3])-offset+1==received && as.numeric(m[4])==total
}
verify_file <- function(path,bytes,md5) {
  file.exists(path)&&is.finite(bytes)&&bytes>0&&file.info(path)$size==bytes&&
    is.character(md5)&&!is.na(md5)&&grepl("^[0-9a-f]{32}$",md5)&&
    identical(unname(tools::md5sum(path)),md5)
}
fetch_one <- function(url,dest,partial,response,bytes,md5,c,transport=curl::curl_fetch_disk) {
  bytes<-as.numeric(bytes) # curl must not receive bit64's underlying double representation.
  if(verify_file(dest,bytes,md5))return(invisible(dest))
  if(!is.finite(bytes)||bytes<=0||is.na(md5)||!grepl("^[0-9a-f]{32}$",md5))
    stop("Verified byte size and MD5 are required before acquisition: ",basename(dest))
  if(!grepl("^https://",url))stop("HTTPS resource required")
  for(p in c(dest,partial,response))dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
  for(attempt in seq_len(c$retrieval$attempts)) {
    offset<-if(file.exists(partial))file.info(partial)$size else 0
    if(offset>=bytes) {
      if(verify_file(partial,bytes,md5)) {
        if(!file.rename(partial,dest))stop("POSIX rename failed");return(invisible(dest))
      };unlink(partial);offset<-0
    }
    h<-curl::new_handle(timeout=c$retrieval$timeout_seconds,
      connecttimeout=c$retrieval$connect_timeout_seconds,followlocation=TRUE,
      maxfilesize=bytes-offset)
    if(offset>0)curl::handle_setopt(h,range=paste0(format(offset,scientific=FALSE),"-"))
    r<-tryCatch(transport(url,response,h),error=function(e)NULL)
    n<-if(file.exists(response))file.info(response)$size else 0
    if(!is.null(r)&&n>0&&range_ok(r$status_code,r$headers,offset,n,bytes)) {
      if(offset==0) {
        if(!file.rename(response,partial))stop("POSIX staging rename failed")
      } else {
        # Appending is restricted to a successfully validated 206 segment on POSIX scratch.
        inp<-file(response,"rb");out<-file(partial,"ab")
        tryCatch(repeat {
          chunk<-readBin(inp,"raw",n=1048576L);if(!length(chunk))break;writeBin(chunk,out)
        },finally={close(inp);close(out)})
        unlink(response)
      }
      if(verify_file(partial,bytes,md5)) {
        if(!file.rename(partial,dest))stop("POSIX final rename failed")
        return(invisible(dest))
      }
      if(file.info(partial)$size>=bytes)unlink(partial)
    } else unlink(response) # HTTP 200 is NEVER appended to an existing partial.
    if(attempt<c$retrieval$attempts)Sys.sleep(c$retrieval$retry_seconds)
  }
  stop("Acquisition retry/budget limit; partial is not validated: ",basename(dest))
}
fastq_id <- function(x)sub("(/[12])$","",sub("[ \t].*$","",substring(x,2L)))
validate_fastq <- function(f1,f2=NULL,limit=Inf,output1=NULL,output2=NULL) {
  a<-gzfile(f1,"rt");b<-if(!is.null(f2))gzfile(f2,"rt") else NULL
  oa<-if(!is.null(output1))gzfile(output1,"wt") else NULL
  ob<-if(!is.null(output2))gzfile(output2,"wt") else NULL
  on.exit({close(a);if(!is.null(b))close(b);if(!is.null(oa))close(oa);if(!is.null(ob))close(ob)})
  n<-0
  repeat {
    if(n>=limit)break
    k<-as.integer(min(10000,limit-n))*4L
    x<-readLines(a,n=k,warn=FALSE);y<-if(!is.null(b))readLines(b,n=k,warn=FALSE) else NULL
    if(!length(x)) {if(!is.null(b)&&length(y))stop("Mate count mismatch");break}
    if(length(x)%%4L||(!is.null(b)&&length(x)!=length(y)))stop("Truncated FASTQ/mate mismatch")
    ix<-seq(1L,length(x),4L)
    for(z in list(x,y))if(length(z)) {
      if(any(!startsWith(z[ix],"@"))||any(!startsWith(z[ix+2L],"+"))||
        any(nchar(z[ix+1L])!=nchar(z[ix+3L])))stop("Malformed FASTQ")
    }
    if(!is.null(y)&&any(fastq_id(x[ix])!=fastq_id(y[ix])))stop("Mate name mismatch")
    if(!is.null(oa))writeLines(x,oa)
    if(!is.null(ob))writeLines(y,ob)
    n<-n+length(ix)
  }
  if(n==0)stop("Empty FASTQ")
  n
}
retrieve_inputs <- function(c,native=TRUE) {
  m<-load_manifests(c)
  if(c$profile=="human_trio")stop("Private raw files are inventoried locally, never fetched by public retrieval")
  records<-list()
  for(i in seq_len(nrow(m$sources)))for(mate in seq_len(if(m$sources$layout[i]=="PE")2 else 1)) {
    s<-m$sources[i];id<-s$run_id
    p<-cfg_path(c,"fastq",run=id,mate=mate)
    fetch_one(s[[paste0("url_",mate)]],p,cfg_path(c,"partial",run=id,mate=mate),
      cfg_path(c,"response",run=id,mate=mate),s[[paste0("bytes_",mate)]],
      s[[paste0("md5_",mate)]],c)
    if(native)run_checked(c$tools$gzip,c("-t",p))
    records[[length(records)+1L]]<-data.table(id=id,mate=mate,path=p,bytes=file.info(p)$size)
  }
  for(i in seq_len(nrow(m$sources))) {
    s<-m$sources[i]
    validate_fastq(cfg_path(c,"fastq",run=s$run_id,mate=1),
      if(s$layout=="PE")cfg_path(c,"fastq",run=s$run_id,mate=2) else NULL)
  }
  r<-resolve_manifest_paths(fread(c$paths$resources),c("destination","partial","response"),c)
  require_columns(r,c("resource_id","url","destination","partial","response","bytes","md5","required","build"))
  assert_unique(r,"resource_id")
  for(i in which(r$required)) {
    fetch_one(r$url[i],r$destination[i],r$partial[i],r$response[i],r$bytes[i],r$md5[i],c)
    if(native&&endsWith(r$destination[i],".gz"))run_checked(c$tools$gzip,c("-t",r$destination[i]))
    if("decompress_to" %in% names(r)&&!is.na(r$decompress_to[i])&&nzchar(r$decompress_to[i])) {
      out<-resolve_value(r$decompress_to[i],dirname(attr(c,"file")),unlist(c$environment_allowlist))
      ensure_parent(out)
      if(native)run_checked(c$tools$gzip,c("-dc",r$destination[i]),stdout=out) else {
        inp<-gzfile(r$destination[i],"rb");dest<-file(out,"wb")
        tryCatch(withCallingHandlers(repeat {
          block<-readBin(inp,"raw",1048576L)
          if(!length(block))break
          writeBin(block,dest)
        },warning=function(w)stop(conditionMessage(w))),finally={close(inp);close(dest)})
      }
    }
    records[[length(records)+1L]]<-data.table(id=r$resource_id[i],mate=NA_integer_,
      path=r$destination[i],bytes=file.info(r$destination[i])$size)
  }
  write_tsv(rbindlist(records),c$paths$resource_receipt)
}
