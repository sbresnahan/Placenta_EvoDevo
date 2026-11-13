range_ok <- function(status,headers,offset,received,total) {
  if(offset==0 && status==200L)return(received==total)
  if(status!=206L)return(FALSE)
  h<-tolower(rawToChar(headers))
  m<-regmatches(h,regexec("content-range:[ \t]*bytes ([0-9]+)-([0-9]+)/([0-9]+)",h))[[1]]
  length(m)==4L && as.numeric(m[2])==offset &&
    as.numeric(m[3])-offset+1==received && as.numeric(m[4])==total
}
verify_file <- function(path,bytes) {
  file.exists(path)&&is.finite(bytes)&&bytes>0&&file.info(path)$size==bytes
}
retrieval_method <- function(c) if(is.null(c$retrieval$method)) "https" else c$retrieval$method
sra_receipt <- function(c,run) paste0(cfg_path(c,"fastq",run=run,mate=1),".sra.tsv")
verify_public_source <- function(c,s,mate) {
  expected_path<-cfg_path(c,"fastq",run=s$run_id,mate=mate)
  expected_mate<-as.character(mate)
  if(retrieval_method(c)=="sra"&&file.exists(sra_receipt(c,s$run_id))) {
    d<-tryCatch(fread(sra_receipt(c,s$run_id),colClasses="character"),error=function(e)NULL)
    if(!is.null(d)&&all(c("run_id","mate","path","bytes","transport") %in% names(d))) {
      expected_run<-s$run_id
      d<-d[run_id==expected_run & mate==expected_mate & path==expected_path & transport=="sra"]
      if(nrow(d)==1L&&grepl("^[1-9][0-9]*$",d$bytes)&&file.exists(expected_path)&&
         file.info(expected_path)$size==as.numeric(d$bytes))return(TRUE)
    }
  }
  # Original ENA downloads are accepted by their declared byte count.
  verify_file(expected_path,s[[paste0("bytes_",mate)]])
}
fetch_one <- function(url,dest,partial,response,bytes,c,transport=curl::curl_fetch_disk) {
  bytes<-as.numeric(bytes) # curl must not receive bit64's underlying double representation.
  if(verify_file(dest,bytes))return(invisible(dest))
  if(!is.finite(bytes)||bytes<=0)
    stop("Expected byte size required before acquisition: ",basename(dest))
  if(!grepl("^https://",url))stop("HTTPS resource required")
  for(p in c(dest,partial,response))dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
  for(attempt in seq_len(c$retrieval$attempts)) {
    offset<-if(file.exists(partial))file.info(partial)$size else 0
    if(offset>=bytes) {
      if(verify_file(partial,bytes)) {
        if(!file.rename(partial,dest))stop("POSIX rename failed");return(invisible(dest))
      };unlink(partial);offset<-0
    }
    h<-curl::new_handle(timeout=c$retrieval$timeout_seconds,
      connecttimeout=c$retrieval$connect_timeout_seconds,followlocation=TRUE,
      maxfilesize_large=bytes-offset)
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
      if(verify_file(partial,bytes)) {
        if(!file.rename(partial,dest))stop("POSIX final rename failed")
        return(invisible(dest))
      }
      if(file.info(partial)$size>=bytes)unlink(partial)
    } else unlink(response) # HTTP 200 is NEVER appended to an existing partial.
    if(attempt<c$retrieval$attempts)Sys.sleep(c$retrieval$retry_seconds)
  }
  stop("Acquisition retry/budget limit; partial is not validated: ",basename(dest))
}
fastq_id <- function(x) {
  id<-sub("(/[12])$","",sub("[ \t].*$","",substring(x,2L)))
  # fasterq-dump can append the mate after the SRA spot number: ERR123.4.1/.2.
  sub("^((SRR|ERR|DRR)[0-9]+\\.[0-9]+)\\.[12]$","\\1",id)
}
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
  if(retrieval_method(c)=="sra")stop("SRA retrieval requires the Seadragon host-only retrieval worker")
  m<-load_manifests(c)
  if(c$profile=="human_trio")stop("Private raw files are inventoried locally, never fetched by public retrieval")
  records<-list()
  for(i in seq_len(nrow(m$sources)))for(mate in seq_len(if(m$sources$layout[i]=="PE")2 else 1)) {
    s<-m$sources[i];id<-s$run_id
    p<-cfg_path(c,"fastq",run=id,mate=mate)
    fetch_one(s[[paste0("url_",mate)]],p,cfg_path(c,"partial",run=id,mate=mate),
      cfg_path(c,"response",run=id,mate=mate),s[[paste0("bytes_",mate)]],c)
    if(native)run_checked(c$tools$gzip,c("-t",p))
    records[[length(records)+1L]]<-data.table(id=id,mate=mate,path=p,bytes=file.info(p)$size)
  }
  for(i in seq_len(nrow(m$sources))) {
    s<-m$sources[i]
    validate_fastq(cfg_path(c,"fastq",run=s$run_id,mate=1),
      if(s$layout=="PE")cfg_path(c,"fastq",run=s$run_id,mate=2) else NULL)
  }
  r<-resolve_manifest_paths(fread(c$paths$resources),c("destination","partial","response"),c)
  require_columns(r,c("resource_id","url","destination","partial","response","bytes","required","build"))
  assert_unique(r,"resource_id")
  for(i in which(r$required)) {
    fetch_one(r$url[i],r$destination[i],r$partial[i],r$response[i],r$bytes[i],c)
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

# Seadragon: describe transfers in R, execute network requests only on the host.
host_retrieval_items <- function(c,resources_only=FALSE) {
  if(c$profile=="human_trio")stop("Private raw files must not be downloaded")
  m<-if(resources_only)NULL else load_manifests(c)
  method<-retrieval_method(c)
  if(!resources_only && method=="sra")require_columns(m$sources,"read_count")
  rows<-list()
  if(!resources_only)for(i in seq_len(nrow(m$sources)))for(mate in seq_len(if(m$sources$layout[i]=="PE")2 else 1)) {
    s<-m$sources[i];id<-s$run_id
    rows[[length(rows)+1L]]<-data.table(id=id,mate=mate,
      url=s[[paste0("url_",mate)]],destination=cfg_path(c,"fastq",run=id,mate=mate),
      partial=cfg_path(c,"partial",run=id,mate=mate),
      bytes=as.numeric(s[[paste0("bytes_",mate)]]),
      decompress_to=NA_character_,method=method,
      read_count=if(method=="sra")as.numeric(s$read_count) else NA_real_,
      sra_receipt=if(method=="sra")sra_receipt(c,id) else NA_character_)
  }
  r<-resolve_manifest_paths(fread(c$paths$resources),
    c("destination","partial","response"),c)
  require_columns(r,c("resource_id","url","destination","partial","response","bytes","required","build"))
  assert_unique(r,"resource_id")
  required_rows<-which(r$required %in% TRUE)
  if(resources_only && !length(required_rows))
    stop("No required reference resources in manifest: ",c$paths$resources,
      ". Populate it with bin/prepare_fish_refseq_resources.sh before retrieval.")
  for(i in required_rows) {
    resource_method<-if("method" %in% names(r)&&!is.na(r$method[i]))r$method[i] else "https"
    if(resource_method=="cyverse_md5")resource_method<-"cyverse" # legacy private manifests
    unzip<-if("decompress_to" %in% names(r)&&!is.na(r$decompress_to[i])&&nzchar(r$decompress_to[i]))
      resolve_value(r$decompress_to[i],dirname(attr(c,"file")),unlist(c$environment_allowlist)) else NA_character_
    rows[[length(rows)+1L]]<-data.table(id=r$resource_id[i],mate=NA_integer_,
      url=r$url[i],destination=r$destination[i],partial=r$partial[i],
      bytes=as.numeric(r$bytes[i]),decompress_to=unzip,
      method=resource_method,
      read_count=NA_real_,sra_receipt=NA_character_)
  }
  items<-rbindlist(rows,fill=TRUE)
  cyverse<-items$method=="cyverse"
  sized<-!cyverse
  if(anyNA(items[,.(id,url,destination,partial,method)])||
     any(!items$method %in% c("https","sra","cyverse"))||
     any(!grepl("^https://",items$url))||
     any(!is.finite(items$bytes[sized])|items$bytes[sized]<=0|
       items$bytes[sized]!=floor(items$bytes[sized])|items$bytes[sized]>2^53)||
     any(!is.na(items$bytes[cyverse])&
       (!is.finite(items$bytes[cyverse])|items$bytes[cyverse]<=0|
        items$bytes[cyverse]!=floor(items$bytes[cyverse])|items$bytes[cyverse]>2^53))||
     anyDuplicated(items$destination)||
     any(cyverse & (!grepl("^https://data[.]cyverse[.]org/dav-anon/",items$url)|
       !grepl("[.]vcf[.]gz$",items$url))))
    stop("Unsafe, duplicated or unsized public acquisition entry")
  sra<-items[method=="sra"]
  if(nrow(sra)&&(any(!grepl("^(SRR|ERR|DRR)[0-9]+$",sra$id))||
      any(!is.finite(sra$read_count)|sra$read_count<=0|sra$read_count!=floor(sra$read_count))||
      anyNA(sra$sra_receipt)))stop("SRA retrieval requires valid run IDs and read counts")
  items[,target:=ifelse(method=="sra",paste0("sra_",id),
    ifelse(!is.na(mate),paste0("fastq_",id,"_",mate),paste0("resource_",id)))]
  if(any(!grepl("^[A-Za-z0-9_.-]+$",items$target)))
    stop("Unsafe retrieval target identifier")
  items
}
prepare_host_retrieval <- function(c,root,job,target="shared") {
  resources_only<-target=="resources_finalize" || startsWith(target,"resource_")
  items<-host_retrieval_items(c,resources_only=resources_only)
  if(!target %in% c("shared","finalize","resources_finalize")) {
    # Resolve the function argument before entering data.table's column scope:
    # `items[target == target]` would otherwise select every resource.
    selection<-which(items[["target"]]==target)
    if(!length(selection))stop("Unknown retrieval target: ",target)
    items<-items[selection]
  }
  write_tsv(items,file.path(job,"retrieve_items.tsv"))
  lines<-character()
  if(target %in% c("finalize","resources_finalize"))return(write_plan(lines,job,"pre"))
  for(i in seq_len(nrow(items))) {
    item<-items[i]
    if(item$method=="sra") {
      if(item$mate!=1L)next # One SRA archive produces both mates.
      mates<-items[id==item$id & method=="sra"]
      second<-mates[mate==2L]
      layout<-if(nrow(second))"PE" else "SE"
      lines<-c(lines,plan_line(file.path(root,"bin/fetch_sra.sh"),c(
        item$id,layout,item$destination,if(nrow(second))second$destination else "-",
        c$paths$scratch,c$retrieval$sra_threads,
        format(item$read_count,scientific=FALSE,trim=TRUE),item$sra_receipt,
        format(item$bytes,scientific=FALSE,trim=TRUE),
        if(nrow(second))format(second$bytes,scientific=FALSE,trim=TRUE) else "0",
        format(c$retrieval$attempts,scientific=FALSE,trim=TRUE),
        format(c$retrieval$retry_seconds,scientific=FALSE,trim=TRUE))))
    } else if(item$method=="cyverse") {
      lines<-c(lines,plan_line(file.path(root,"bin/fetch_cyverse.sh"),
        c(c$tools$curl,item$url,item$destination,item$partial,
          c$retrieval$attempts,c$retrieval$timeout_seconds,
          c$retrieval$connect_timeout_seconds,c$retrieval$retry_seconds)))
    } else {
      lines<-c(lines,plan_line(file.path(root,"bin/fetch_verified.sh"),
        c(c$tools$curl,item$url,item$destination,item$partial,
          format(item$bytes,scientific=FALSE,trim=TRUE),
          c$retrieval$attempts,c$retrieval$timeout_seconds,
          c$retrieval$connect_timeout_seconds,c$retrieval$retry_seconds)))
    }
    if(!is.na(item$decompress_to)) {
      ensure_parent(item$decompress_to)
      lines<-c(lines,plan_line(c$tools$gzip,c("-dc",item$destination),
        stdout=item$decompress_to))
    }
  }
  write_plan(lines,job,"pre")
}
finish_host_retrieval <- function(c,job,target="shared") {
  items<-fread(file.path(job,"retrieve_items.tsv"),colClasses=c(bytes="numeric"))
  m<-if(any(items$method=="sra"))load_manifests(c) else NULL
  for(i in seq_len(nrow(items))) {
    item<-items[i]
    valid<-if(item$method=="sra") {
      s<-m$sources[run_id==item$id]
      nrow(s)==1L&&verify_public_source(c,s,item$mate)
    } else if(item$method=="cyverse") {
      receipt<-paste0(item$destination,".cyverse.tsv")
      d<-tryCatch(fread(receipt,colClasses="character"),error=function(e)NULL)
      if(is.null(d)||!all(c("path","bytes","url","method") %in% names(d))||nrow(d)!=1L)
        FALSE
      else if(d$path!=item$destination||d$url!=item$url||
        !d$method %in% c("cyverse","cyverse_md5")||!grepl("^[1-9][0-9]*$",d$bytes))FALSE
      else if(!file.exists(item$destination)||file.info(item$destination)$size!=as.numeric(d$bytes))FALSE
      else {items$bytes[i]<-as.numeric(d$bytes);TRUE}
    } else verify_file(item$destination,item$bytes)
    if(!valid)
      stop("Host retrieval missing/invalid: ",basename(item$destination))
    if(!is.na(item$decompress_to)&&
       (!file.exists(item$decompress_to)||file.info(item$decompress_to)$size<=0))
      stop("Host decompression missing/empty: ",basename(item$decompress_to))
  }
  if(!target %in% c("shared","finalize","resources_finalize"))return(invisible(TRUE))
  # Transfer tools already completed. Do not re-read entire FASTQs at finalization.
  items[,bytes:=vapply(destination,function(p)file.info(p)$size,numeric(1))]
  write_tsv(items[,.(id,mate,path=destination,bytes,method)],c$paths$resource_receipt)
}
