# Deterministic software/coordinate fixtures only. No biological or benchmark data.
make_host_fixture <- function(root,dir) {
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  put<-function(name,lines) {
    p<-file.path(dir,name);writeLines(lines,p);Sys.chmod(p,"0755");p
  }
  module<-put("modules.sh",c("# deterministic module initialization fixture",
    paste0("export PATH=",shQuote(dir),':"$PATH"')))
  conda<-put("conda",c("#!/bin/bash",
    '[[ "$*" == "shell.bash hook" ]] || exit 9',
    "printf '%s\\n' 'conda() { [[ \"$1\" == activate ]] || return 8; export CONDA_PREFIX=\"/fixture/$2\"; export POE_FIXTURE_ACTIVE=\"$2\"; }'"))
  real_samtools<-Sys.which("samtools")
  if(!nzchar(real_samtools))stop("Tiny coordinate integration requires local samtools")
  samtools<-put("samtools",c("#!/bin/bash","set -euo pipefail",
    '[[ "${POE_FIXTURE_ACTIVE:-}" == samtools-1.16.1 ]] || exit 71',
    paste("exec",shQuote(real_samtools),'"$@"')))
  bcftools<-put("bcftools",c("#!/bin/bash","set -euo pipefail",
    '[[ "${POE_FIXTURE_ACTIVE:-}" == bcftools-1.16 ]] || exit 72',
    '[[ "${1:-}" != --version ]] || { printf "MOCK_bcftools_not_real_version\\n"; exit 0; }',
    'op="$1"; shift; dest=""; selected=""; input="${!#}"',
    'while [[ $# -gt 0 ]]; do case "$1" in -o) dest="$2"; shift 2 ;; -s) selected="$2"; shift 2 ;; *) shift ;; esac; done',
    'case "$op" in',
    'view|norm) cp "$input" "$dest" ;;',
    "query) awk -v selected=\"$selected\" 'BEGIN{FS=OFS=\"\\t\"; n=split(selected,ids,\",\")} /^#CHROM/{for(j=1;j<=n;j++)for(k=10;k<=NF;k++)if($k==ids[j])idx[j]=k;next} !/^#/{printf \"%s\\t%s\\t%s\\t%s\",$1,$2,$4,$5;for(j=1;j<=n;j++){if(!idx[j])exit 81; printf \"\\t%s\",$idx[j]}printf \"\\n\"}' \"$input\" > \"$dest\" ;;",
    '*) exit 73 ;; esac'))
  ornaments<-put("ornaments",c("#!/bin/bash","set -euo pipefail",
    '[[ "${POE_FIXTURE_ACTIVE:-}" == ornaments ]] || exit 74',
    'op="$1"; shift',
    '[[ "$op" != quant || "${POE_FIXTURE_FAIL:-no}" != yes ]] || exit 83',
    '[[ "$op" != version ]] || { printf "MOCK_ORNAMENTS_NOT_A_REAL_BUILD\\n"; exit 0; }',
    'dest=""; index=""; vcf=""; input="${!#}"',
    'while [[ $# -gt 0 ]]; do case "$1" in -o) dest="$2"; shift 2 ;; -i) index="$2"; shift 2 ;; --vcf) vcf="$2"; shift 2 ;; *) shift ;; esac; done',
    'case "$op" in',
    'index) printf "DETERMINISTIC_MOCK_INDEX\\n" > "$index" ;;',
    'quant) mkdir -p "$dest"',
    "awk 'BEGIN{FS=OFS=\"\\t\"}!/^#/{print $1,$2,$4,$5,\".\",12,8,0}' \"$vcf\" > \"$dest/allele_counts.txt\"",
    'printf "MOCK_THETAS\\n" > "$dest/thetas.txt"; printf "MOCK_TPMS\\n" > "$dest/tpms.txt"',
    'printf "[quant] processed 3 reads, 3 reads pseudoaligned\\n" ;;',
    '*) exit 75 ;; esac'))
  personalizer<-put("create_personalized_transcriptome",c("#!/bin/bash","set -euo pipefail",
    '[[ "${POE_FIXTURE_ACTIVE:-}" == ornaments ]] || exit 76',
    'input=""; dest=""; while [[ $# -gt 0 ]]; do case "$1" in -f) input="$2"; shift 2 ;; -o) dest="$2"; shift 2 ;; *) shift ;; esac; done',
    'cp "$input" "$dest"'))
  python<-put("python",c("#!/bin/bash",'printf "MOCK_PYTHON_NOT_EXECUTED\\n"'))
  activate<-put("activate",c(paste0("export PATH=",shQuote(dir),':"$PATH"'),
    "export POE_FIXTURE_ACTIVE=ornaments"))
  # Known fixture resource values, never used as biological performance estimates.
  timer<-put("time",c("#!/bin/bash","set -euo pipefail",
    '[[ "$1" == -v && "$2" == -o ]] || exit 77; out="$3"; shift 3',
    '"$@"',
    "printf '\\tMaximum resident set size (kbytes): 128\\n\\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:01.00\\n' > \"$out\""))
  singularity<-put("singularity",c("#!/bin/bash","set -euo pipefail",
    '[[ -z "${LD_LIBRARY_PATH:-}" && -z "${PYTHONPATH:-}" && -z "${PYTHONHOME:-}" ]] || exit 78',
    '[[ "$1" == exec ]] || exit 79; shift',
    'while [[ $# -gt 0 ]]; do case "$1" in --cleanenv) shift ;; --bind) shift 2 ;; --env) export "$2"; shift 2 ;; *) shift; break ;; esac; done',
    '[[ "$1" == Rscript ]] || exit 80; shift',
    paste("exec",shQuote(file.path(R.home("bin"),"Rscript")),'"$@"')))
  sif<-put("fixture.sif","NOT_A_REAL_CONTAINER")
  site<-put("site.sh",c(
    paste0("export POE_MODULE_INIT=",shQuote(module)),
    paste0("export POE_CONDA=",shQuote(conda)),
    "export POE_SAMTOOLS_ENV=samtools-1.16.1","export POE_BCFTOOLS_ENV=bcftools-1.16",
    paste0("export POE_ORNAMENTS_ACTIVATE=",shQuote(activate)),
    "export POE_ORNAMENTS_SOURCE=''",
    paste0("export POE_R_SIF=",shQuote(sif)),
    paste0("export POE_R_LIBRARY=",shQuote(.libPaths()[1])),
    paste0("export POE_SINGULARITY=",shQuote(singularity))))
  list(dir=dir,site=site,sif=sif,ornaments=ornaments,personalizer=personalizer,time=timer,
    bcftools=bcftools,samtools=samtools,activate=activate)
}
make_profile_fixture <- function(root,dir,host,human=FALSE) {
  dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  original<-file.path(root,"config",if(human)"GUSTO.example.yaml" else "maize_B73_W22.example.yaml")
  x<-yaml::read_yaml(original);x$cohort_id<-if(human)"human_fixture" else "cross_fixture"
  for(k in names(x$paths))x$paths[[k]]<-file.path(dir,"products",k)
  x$paths$genome<-file.path(dir,"genome.fa");x$paths$annotation<-file.path(dir,"annotation.gtf")
  x$paths$genome_fai<-paste0(x$paths$genome,".fai")
  x$paths$transcripts_fai<-paste0(x$paths$transcripts,".fai")
  x$paths$job_registry<-file.path(dir,"jobs.tsv")
  for(k in names(x$templates))x$templates[[k]]<-file.path(dir,"templates",paste0(k,"_{sample}_{mode}_{strand}"))
  # Preserve output contracts and injectivity across library/mode/strand/run.
  qp<-file.path(dir,"quant/{sample}/{mode}/{strand}")
  x$templates$quant_dir<-qp
  for(k in c("allele_counts","thetas","tpms"))x$templates[[k]]<-file.path(qp,paste0(k,".txt"))
  for(k in c("event_counts","gene_counts","assignment_qc","quant_log","resource_log"))
    x$templates[[k]]<-file.path(qp,paste0(k,".tsv"))
  for(k in c("quant_receipt","assign_receipt"))x$templates[[k]]<-file.path(qp,paste0(k,".rds"))
  for(k in c("child_receipt","human_receipt","library_pilot_receipt","library_pilot_report","library_metrics"))
    x$templates[[k]]<-file.path(dir,"per_library/{sample}",k)
  for(k in c("child_query","personal_query","child_vcf","personal_vcf","child_tx_vcf","child_tx_map","child_fasta","child_index"))
    x$templates[[k]]<-file.path(dir,"per_child/{sample}",k)
  x$templates$prefix<-file.path(dir,"prefix/{sample}/{mode}/{strand}/{run}_{mate}.gz")
  for(k in c("fastq","partial","response"))x$templates[[k]]<-file.path(dir,"raw/{run}_{mate}.gz")
  x$execution$site_file<-host$site;x$execution$posix_work_approved<-TRUE
  x$execution$image<-host$sif
  x$reference$genome_source<-NULL;x$reference$annotation_approved<-TRUE
  x$reference$contig_aliases<-list()
  x$markers$mode<-"prepared";x$markers$proxy<-FALSE;x$markers$provenance_approved<-TRUE
  x$ornaments$binary_sha256<-sha256(host$ornaments)
  x$ornaments$personalizer_sha256<-sha256(host$personalizer)
  # Resolve repository tools against the actual checkout, not fixture config location.
  for(k in names(x$tools))if(startsWith(x$tools[[k]],"../bin/"))
    x$tools[[k]]<-file.path(root,sub("^\\.\\./","",x$tools[[k]]))
  x$tools$r_library<-.libPaths()[1];x$tools$time<-host$time
  x$pilot$fragments<-3;x$pilot$strand_fragments<-3
  x$pilot$max_projected_seconds<-100;x$pilot$max_work_bytes<-1e8
  x$human$genotype_qc_approved<-TRUE;x$human$covariates_approved<-TRUE
  x$lsf$resource_limits_approved<-TRUE
  for(k in c("memory_units_verified","memory_per_job_verified","execution_environment_verified"))x$lsf[[k]]<-TRUE
  dir.create(dirname(x$paths$work),recursive=TRUE,showWarnings=FALSE)
  genome<-paste(rep("A",200),collapse="")
  writeLines(c(">chrT",genome),x$paths$genome)
  gtf<-character()
  for(t in c("txA","txB")) {
    starts<-c(if(t=="txA")1L else 61L,101L);ends<-starts+7L
    attr<-sprintf('gene_id "GENE"; gene_name "GENE"; transcript_id "%s";',t)
    gtf<-c(gtf,paste("chrT","fixture","transcript",min(starts),max(ends),".","+",".",attr,sep="\t"),
      paste("chrT","fixture","exon",starts,ends,".","+",".",paste0(attr,' exon_number "',1:2,'";'),sep="\t"))
  }
  writeLines(gtf,x$paths$annotation)
  x$reference$genome_sha256<-sha256(x$paths$genome);x$reference$annotation_sha256<-sha256(x$paths$annotation)
  samples<-c("fixtureA","fixtureB")
  libs<-data.table(sample_id=samples,unit_id=samples,layout="PE",strand_mode="fr",strand_verified=TRUE)
  if(human)libs[,`:=`(child_id=c("childA","childB"),identity_verified=TRUE)]
  units<-data.table(unit_id=samples,cluster_id=samples,independence_verified=TRUE,
    maternal_parent=if(human)"typed_mother" else c(x$parents$focal,x$parents$alternate),
    paternal_parent=if(human)"typed_father" else c(x$parents$alternate,x$parents$focal),
    stage="FIXTURE",tissue="FIXTURE")
  source<-data.table(run_id=samples,sample_id=samples,layout="PE",url_1="",url_2="",md5_1="",md5_2="",
    bytes_1=0,bytes_2=0,read_count=6)
  for(i in 1:2)for(k in 1:2) {
    p<-file.path(dir,"raw",paste0(samples[i],"_",k,".gz"));dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
    con<-gzfile(p,"wt");writeLines(unlist(lapply(1:3,function(j)c(paste0("@fixture",j,"/",k),"AAAA","+","!!!!"))),con);close(con)
    source[i,(paste0("md5_",k)):=unname(tools::md5sum(p))]
    source[i,(paste0("bytes_",k)):=file.info(p)$size]
    source[i,(paste0("local_",k)):=p]
  }
  write_tsv(source,x$paths$sources);write_tsv(libs,x$paths$libraries);write_tsv(units,x$paths$units)
  resource<-data.table(resource_id=character(),url=character(),destination=character(),partial=character(),
    response=character(),bytes=numeric(),md5=character(),required=logical(),build=character())
  write_tsv(resource,x$paths$resources)
  markers<-data.table(chrom="chrT",genomic_pos=c(2L,62L),ref="A",alt="G",focal_index=0L,alternate_index=1L)
  mp<-file.path(dir,"fixture_markers.tsv");write_tsv(markers,mp)
  write_tsv(data.table(source_id="coordinate_fixture",path=mp,build=x$reference$build,sha256=sha256(mp)),x$paths$marker_sources)
  if(human) {
    write_tsv(data.table(child_id=c("childA","childB"),mother_id=c("motherA","motherB"),
      father_id=c("fatherA","fatherB"),family_id=c("familyA","familyB")),x$paths$trios)
    write_tsv(data.table(sample_id=samples),x$paths$covariates)
    header<-c("##fileformat=VCFv4.2",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchildA\tmotherA\tfatherA\tchildB\tmotherB\tfatherB")
    records<-paste("chrT",c(2L,62L),".","A","G",".","PASS",".","GT",
      "0/1","0/0","1/1","0/1","1/1","0/0",sep="\t")
    for(p in c(x$paths$typed_vcf,x$paths$personalization_vcf))writeLines(c(header,records),p)
  }
  cfg<-file.path(dir,"profile.yaml");yaml::write_yaml(x,cfg)
  list(config=cfg,c=read_config_v2(cfg),samples=samples)
}
