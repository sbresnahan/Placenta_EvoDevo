# 06 — Models and callability

**Purpose:** assess eligible gene/promoter comparisons and fit public reciprocal
or private paired-family models. Failed gates produce statuses, not forced
p-values.

**Prerequisites:** all configured primary units and valid reference, pilot,
scale, quantification and assignment receipts; reviewed independence and,
for human, approved linked covariates.

From the script root, after scale completes:

```bash
bin/submit_pipeline.sh "$CONFIG" models
bin/submit_pipeline.sh "$CONFIG" models --submit
```

Wait for the model job, then:

```bash
bin/Rscript_sif poe.R "$CONFIG" status
```

**Outputs:** `paths.model_gene`, `paths.model_promoter`, `paths.model_pairs`
and model receipt. Human gene analysis is deferred; inspect promoter results.
Public results require reciprocal direction/rank/independence/coverage checks.
Human inference also requires the specified family/orientation, residual-df
and bootstrap-validity floors.

**Next:** interpret eligible statuses and qualifications. GUSTO aggregates
may proceed to separate [local disclosure review](../07_exports/README.md);
there is no automatic export.
