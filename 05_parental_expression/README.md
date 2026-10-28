# 05 — Parental expression

**Purpose:** orient expected allelic allocations as maternal/paternal expression
at the same genomic locus.

**Prerequisites:** validated complete-library quantification, correct coordinate
maps and public parental markers or independently QC'd typed Mendelian anchors.

**Inputs:** count/quantification receipt, exon/transcript/event maps, parent
orientation and exact primary-library identity.

The scale graph submits `01_assign_parents.lsf` after its library's
quantification. No additional manual command is needed. The graph is submitted
once using the accepted report hash:

```bash
bin/submit_pipeline.sh "$CONFIG" scale --submit --accept "$PILOT_SHA"
```

`PILOT_SHA` must already contain the exact reviewed report hash as in the
[walkthrough](../docs/WALKTHROUGH.md). Do not resubmit scale to run assignment
again without first reconciling jobs.

**Outputs:** `templates.gene_counts`, `templates.event_counts`,
`templates.assignment_qc` and per-library assignment receipts. Missing parental
evidence stays NA, not zero; neither repeated SNPs nor technical runs inflate n.

**Next:** wait for all assignments and the scale collector, then
[models](../06_promoter_models/README.md).
