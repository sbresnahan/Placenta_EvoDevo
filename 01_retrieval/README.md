# 01 — Public retrieval

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

Downloads public FASTQs and resources marked `required=TRUE` from the source
and resource tables. Checks MD5/bytes, gzip integrity and FASTQ structure/mates.
Writes raw/reference files and `resource_receipt.tsv`.

### NCBI SRA alternative on Seadragon

For a public cohort, set `retrieval.method: sra` and `retrieval.sra_threads: 4`
in its **private** YAML; allocate at least four CPUs to `lsf.stages.retrieve.cpus`.
The `maize_B73_Ki11.example.yaml` demonstrates these settings. Existing private
YAMLs are not automatically updated by changing an example.

The LSF transfer-node Bash worker loads `sratoolkit`, calls `prefetch` with the
**original** `SRR`/`ERR`/`DRR` run ID, validates the SRA archive, then calls
`fasterq-dump --split-files --threads 4 --temp <cohort-scratch>/sra_tmp/...`.
Its output is gzip-compressed into the original `{run}_{mate}.fastq.gz` paths.
Reference resources still use verified HTTPS downloads. No R process retrieves
data or submits LSF jobs.

SRA-derived compressed files will not have ENA's original gzip MD5/size.
The worker records the generated files' own size/MD5 in a per-run `.sra.tsv`
receipt, checks `gzip -t`, and R rechecks the receipt, full FASTQ structure,
mate pairing and expected read count. Downstream quantification verifies the
same local checksum receipt. A missing SRA accession, partial/missing mate,
unpaired reads or unverified pre-existing destination fails explicitly.

NCBI's `prefetch` and `fasterq-dump` may need substantial space for the `.sra`
cache and **uncompressed** FASTQ/temp files; reserve several times the FASTQ
size under `paths.scratch`, not `$HOME` or a default `/tmp`. Completed `.sra`
archives remain in `scratch/sra_cache` for restartability; delete only after
verified completion if space is needed. Run time can exceed the example's
two-hour LSF transfer wall limit: confirm your site's transfer-queue limits
and submit sufficiently sized batches. The real availability of each ERR run
at NCBI and transfer-node access to NCBI must be checked by the first live job.

Next: [reference preparation](../02_reference/README.md).
Human inputs use local file paths in their source table.
