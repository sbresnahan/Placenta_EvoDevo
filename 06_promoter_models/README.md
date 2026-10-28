# 06 — Models

After complete-library quantification and assignment:

```bash
bin/submit_pipeline.sh "$CONFIG" models --submit
```

Public models use reciprocal-cross orientation and dosage-adjusted gene/promoter
contrasts. Cluster IDs determine independent-unit eligibility.
Human promoter models use paired-family contrasts, cis orientation and configured
covariates.

Outputs: `model_gene.tsv`, `model_promoter.tsv`, `model_pairs.tsv`, and the
model completion record. Each result records eligibility status.
Human gene modeling records `deferred_unimplemented`.

Optional next action: [human aggregate export](../07_exports/README.md).
