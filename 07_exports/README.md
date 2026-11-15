# 07 — Human aggregate export

Set `exports.min_cell` to a positive integer, then prepare the candidate:

```bash
bin/Rscript_sif 07_exports/prepare_candidate.R "$CONFIG"
less -S "$COHORT_DIR/candidate.tsv"
```

Set `exports.enabled: true` and write the aggregate:

```bash
bin/Rscript_sif 07_exports/release_aggregates.R "$CONFIG"
```

These actions validate annotation identities, allowed columns, eligible model
status and cell/family thresholds. Outputs are `paths.candidate` and
`paths.release`.

See the [GUSTO recipe](../docs/Placenta_EvoDevo_step_by_step.md#part-6-human-gusto--private-trio-linked-placenta-analysis).
