# 01 — Public retrieval

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

Prepares a host-only `HOST_PLAN` script; run the printed script with `bash` to
submit one transfer job per SRA run (both mates), or per HTTPS file, followed
by a final receipt job dependent on all acquisitions. The default concurrency
cap is `retrieval.max_parallel_jobs: 2` to limit simultaneous scratch use;
raise it only after checking available storage. HTTPS files are checked against
expected byte sizes during transfer. SRA conversion requires successful commands,
nonempty expected mate files and a size-only local receipt, without a second
pass through the FASTQs. The final job checks file presence and sizes and
writes `resource_receipt.tsv` without rereading entire FASTQs.

When the retrieval jobs finish, confirm the cohort's inputs are ready:

```bash
tests/test_inputs.sh "$CONFIG"
```

### NCBI SRA alternative on Seadragon

For a public cohort, set `retrieval.method: sra` and `retrieval.sra_threads: 4`
in its **private** YAML; allocate at least four CPUs to `lsf.stages.retrieve.cpus`.
The `maize_B73_Ki11.example.yaml` demonstrates these settings. Existing private
YAMLs are not automatically updated by changing an example.

The LSF transfer-node Bash worker loads `sratoolkit`, calls `prefetch` with the
**original** `SRR`/`ERR`/`DRR` run ID, then calls
`fasterq-dump --split-files --threads 4 --temp <cohort-scratch>/sra_tmp/...`.
Failed `prefetch` calls retry up to `retrieval.attempts` times with increasing
delays based on `retrieval.retry_seconds`; cached SRA data is preserved between
attempts and subsequent submissions.
Its output is gzip-compressed into the original `{run}_{mate}.fastq.gz` paths.
Reference resources use HTTPS downloads with byte-size checks, except CyVerse
VCFs, which use gzip-integrity checks and can be staged in the cohort's
`reference/cyverse_downloads/` directory under their original filenames.
No R process retrieves data or submits LSF jobs.

SRA-derived compressed files will not have ENA's original compressed byte size.
The worker records generated FASTQ sizes in a per-run `.sra.tsv` receipt and
does not perform MD5, gzip-integrity or full FASTQ-content scans afterward.
Existing legacy per-run receipts are accepted by file size; this trades away
silent-corruption detection for speed. A missing SRA accession, missing mate,
unexpected unpaired reads, or unreceipted existing destination still fails.
When `pigz` is installed, compression uses the configured SRA thread count;
otherwise it falls back to single-threaded `gzip`.

NCBI's `prefetch` and `fasterq-dump` may need substantial space for the `.sra`
cache and **uncompressed** FASTQ/temp files; reserve several times the FASTQ
size under `paths.scratch`, not `$HOME` or a default `/tmp`. Completed `.sra`
archives remain in `scratch/sra_cache` for restartability; remove only after
all required outputs have been produced if space is needed. Individual run
time can exceed the example's two-hour LSF transfer wall limit: confirm your
site's transfer-queue limits. The real availability of each ERR run
at NCBI and transfer-node access to NCBI must be checked by the first live job.

Next: [reference preparation](../02_reference/README.md).
Human inputs use local file paths in their source table.
