# 01 — Public retrieval

**Purpose:** acquire only the configured public FASTQs/resources with verified
size/checksum and integrity contracts. There is no human-data retrieval route.

**Prerequisites:** successful preflight, approved internet-enabled execution
host/queue, writable cohort storage and complete required asset contracts.

**Inputs:** operational sources/resources TSVs, URLs, expected sizes/checksums
and absolute destinations under the project data root.

From the script root, inspect then submit:

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
```

Wait for the job and inspect logs. R handles bounded HTTPS/range resume;
host gzip verifies integrity. No automatic trimming, SRA conversion or count
table fallback is supplied. Alternatively, run
`01_retrieval/01_retrieve_public.lsf` directly on the approved host with the
same explicit environment, not in addition to a submitted job.

**Outputs:** validated cohort raw/reference files and `paths.resource_receipt`.
Missing required checksum/size contracts block acquisition.

**Next:** review the actual genome/annotation and compatible marker-source
paths/builds, then [reference preparation](../02_reference/README.md).
