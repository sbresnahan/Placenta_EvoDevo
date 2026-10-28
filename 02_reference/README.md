# 02 — Reference catalog

**Purpose:** create indexed genome/transcript references, normalized annotation
and promoter/event maps.

**Prerequisites:** successful preflight; actual compatible genome/annotation;
annotation review; complete downstream public
marker contracts or private typed-child inputs.

**Inputs:** `paths.genome`, `paths.annotation`, build/contig settings and
reference approvals. Human `reference.genome_source` is copied to
`paths.genome` before indexing so the shared source is not rewritten.

From the script root:

```bash
bin/submit_pipeline.sh "$CONFIG" reference
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

This submits a **whole reference graph**, not just this module:
catalog → public markers → shared index, or catalog → private child preparation
and index per library. Do not submit those jobs again separately.

**Outputs:** adjacent FASTA indices, `paths.gtf`, `paths.transcripts`,
`paths.exons`, `paths.tx_map`, `paths.events`, `paths.event_map`, `paths.pairs`,
`paths.annotation_audit`, deferred-endpoint table and catalog receipt.
Catalog completion is not public index completion.

**Next:** wait for [genotypes](../03_genotypes/README.md) and
[index jobs](../04_ornaments/README.md), then run the pilot.
