# GUSTO sample-agnostic quantification recovery

Use this only after all full-scale GUSTO quantification/assignment jobs have
finished or failed. It is for *partial successful* scale executions; do not use
it to conceal a widespread processing problem.

## Behavior

- Discovers failures from *actual per-library receipts and count files*, without
  hard-coded sample IDs or assumptions about the number of failed samples.
- Excludes the entire biological unit if any library of that unit is incomplete.
- Retains only units whose `quant_receipt` is readable and contains metrics,
  and whose `assign_receipt` validates the gene and event count outputs.
- Defaults to at most 10 excluded units and at least 90% of biological units
  retained. Override the maximum explicitly with `--max-excluded N`, up to 100.
- Writes a dated audit and original manifest backups under
  `data/GUSTO/private/force_quant_*/`. Only the five *cohort-local* active
  manifests (sources, libraries, units, trios, covariates) are filtered.
- Does not delete the original data, rerun quantification, fabricate completion
  receipts, or change the source manifests committed under `scripts/`.
- Launches the normal `collect-quantification` worker without the now-invalid
  dependency on the original failed jobs. The worker must produce the actual
  quantification completion receipt.

## Commands (Seadragon)

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export SCRIPTS_DIR="$PROJECT/scripts"
export COHORT_DIR="$PROJECT/data/GUSTO"
export CONFIG="$SCRIPTS_DIR/config/GUSTO.private.yaml"
cd "$SCRIPTS_DIR"

# Review discovery without modifying files.
bin/gusto_force_quantification.sh "$CONFIG" --preview

# ONLY if all GUSTO sample jobs are terminal, cancel obsolete invalid dependencies:
bjobs -p    # inspect pending dependencies; cancel obsolete jobs by verified IDs
# bkill <stale_assignment_job_ids> <stale_collector_job_id>

# Excludes precisely the incomplete units and runs the collector.
bin/gusto_force_quantification.sh "$CONFIG" --force --submit
# Optional safeguard: add --expect-excluded N for the number in the preview

# Once collector DONE:
bin/Rscript_sif poe.R "$CONFIG" status

# If the first --force successfully filtered manifests but collector submission
# failed, DO NOT re-force; resume without changing manifests:
bin/gusto_force_quantification.sh "$CONFIG" --verify --submit

# Proceed ONLY after quantification reports complete.
bin/submit_pipeline.sh "$CONFIG" models --submit
```

If the discovery output is not exactly the expected failures, investigate. In
particular, missing `assign_receipt` files for successful quantification jobs may
indicate still-running assignments rather than unusable samples.

The filtered cohort sample size must be reported when interpreting GUSTO results.
