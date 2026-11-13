# 06 — Models

After complete-library quantification and assignment:

```bash
bin/submit_pipeline.sh "$CONFIG" models --submit
```

Public models use reciprocal-cross orientation and dosage-adjusted gene/promoter
contrasts. Under the exploratory callability policy, the coverage-selected
sentinel requires ≥10 total assigned counts in ≥3 units, including both cross
directions for inference. Cluster IDs are retained as provenance but do not
block model fitting; dependence between sibling embryos or shared ears remains
an uncorrected standard-error limitation. Only the models stage needs rerunning
after a threshold/independence policy change if quantification is complete.
Numerically degenerate residual variance is recorded as
`degenerate_residual_variance`: its effect estimate is retained, but the SE,
confidence interval and p/q-values are unavailable. This prevents spurious
small p-values caused by tied ratios in three-unit fits. For the recommended
primary sensitivity, set `models.min_callable_units: 4` while retaining the
three-unit results as exploratory. Separate mouse X/mitochondrial effects
from autosomal imprinting candidates when interpreting results.
Human promoter models use paired-family contrasts, cis orientation and configured
covariates.

Outputs: `model_gene.tsv`, `model_promoter.tsv`, `model_pairs.tsv`, and the
model completion record. Each result records eligibility status.
Human gene modeling records `deferred_unimplemented`.

Optional next action: [human aggregate export](../07_exports/README.md).
