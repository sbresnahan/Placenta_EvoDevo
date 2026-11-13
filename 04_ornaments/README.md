# 04 — Ornaments

The reference graph builds shared public or per-child human indices.
Native `create_personalized_transcriptome` and `ornaments` run on the host;
R prepares arguments and validates results.

## Pilot

```bash
bin/submit_pipeline.sh "$CONFIG" pilot --submit
```

Each library runs strand diagnostics and a bounded throughput prefix.
The collector writes `pilot_report.tsv` with pseudoalignment, time, RSS,
index size, endpoint work bytes, and informative runtime projections.

## Full quantification

After the pilot completes:

```bash
bin/submit_pipeline.sh "$CONFIG" scale --submit
```

Outputs: `allele_counts.txt`, `thetas.txt`, `tpms.txt`, logs and metrics.
The graph also submits dependent [parental assignment](../05_parental_expression/README.md)
and an all-library collector.
