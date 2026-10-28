# 04 — Ornaments indices and quantification

**Purpose:** build public shared or human child-specific indices, validate
bounded pilots and quantify complete primary libraries.

**Prerequisites:** validated transcript/reference/SNP inputs, native identity,
actual raw reads and approved protocol/resources. Scale also needs the exact
accepted pilot-report hash.

**Execution:** native installed `create_personalized_transcriptome` and
`ornaments` run on the host after activation of the configured environment.
R prepares arguments and validates outputs through `bin/Rscript_sif`.

The reference graph already submits `01_build_index.lsf`. From the script root,
after that graph finishes:

```bash
bin/submit_pipeline.sh "$CONFIG" pilot
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

Wait for all library jobs and `10_collect_pilot.lsf`; inspect the report and
measured resources. Pilot completion never submits scale. After review:

```bash
read -r -p "Accepted pilot-report SHA256: " PILOT_SHA
[[ "$PILOT_SHA" =~ ^[a-f0-9]{64}$ ]] || { printf 'Invalid SHA256\n' >&2; exit 1; }
bin/submit_pipeline.sh "$CONFIG" scale --accept "$PILOT_SHA"
bin/submit_pipeline.sh "$CONFIG" scale --submit --accept "$PILOT_SHA"
```

**Outputs:** public `paths.personal_fasta`/`paths.index` or private child index
templates; per-library pilot reports; all-library `paths.pilot_report`;
quantifier `allele_counts.txt`, `thetas.txt`, `tpms.txt`; logs, metrics and
validated receipts. `11_collect_quantification.lsf` requires quantification
and parental assignment for every configured library.

`loop_ornaments.sh` is an alternative host submission wrapper for pilot/scale,
not an extra processing step. Use one interface, not both.

**Next:** [parental assignment](../05_parental_expression/README.md) runs
automatically per library in scale. After its collector completes, run models.
Expected allocations are not observed joint parent/event-informative molecules.
