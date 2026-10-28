# 01 — Public retrieval

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

Downloads public FASTQs and resources marked `required=TRUE` from the source
and resource tables. Checks MD5/bytes, gzip integrity and FASTQ structure/mates.
Writes raw/reference files and `resource_receipt.tsv`.

Next: [reference preparation](../02_reference/README.md).
Human inputs use local file paths in their source table.
