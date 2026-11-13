# GUSTO `--force` pilot recovery

Run only when all pilot child jobs have left the queue and a small subset failed. This is
an explicitly opt-in **sample exclusion**, not a waiver of QC for included samples.

From the Seadragon login host:

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export SCRIPTS_DIR="$PROJECT/scripts"
export COHORT_DIR="$PROJECT/data/GUSTO"
export CONFIG="$SCRIPTS_DIR/config/GUSTO.private.yaml"
cd "$SCRIPTS_DIR"

# Confirm all pilot jobs have reached a terminal status.
bjobs -a -w -u "$USER" | grep 'GUSTO_pilot_' | awk '{n[$3]++}END{for (k in n) print k,n[k]}'

# Check each included child has a valid .rds pilot receipt and per-child QC TSV.
bin/gusto_force_pilot.sh "$CONFIG" --preview

# Back up cohort-local manifests, exclude only invalid pilot samples,
# and submit a new collector with no dependency on the two failed jobs.
bin/gusto_force_pilot.sh "$CONFIG" --force --submit
```

If submission fails after filtering, resume the collector with
`bin/gusto_force_pilot.sh "$CONFIG" --verify --submit`. This verifies all
remaining child receipts before submitting a new collector, without filtering
the manifests again. Do not resubmit if the collector is already active.

The script requires at least 90% of samples to have valid completed pilot receipts,
with at most 5 excluded by default. Increase the explicit cap only if justified:
`--max-excluded N` (max 50). This modifies ONLY the **active cohort-local**
`{sources,libraries,units,trios,covariates}.tsv` manifests; original repository
manifests and all raw/reference/index/pilot results are unchanged. Backups and
an included/excluded sample audit are retained in `data/GUSTO/private/force_pilot_*/`.
The new collector runs the existing `runtime.R collect-pilot finish` action and is
responsible for writing a genuine `pilot_report.tsv` and `pilot_receipt.rds`.
The script explicitly refuses to fabricate a completion receipt.

Inspect the new collector job and report:

```bash
bjobs -a -u "$USER" | grep GUSTO_force_collect_pilot
bin/Rscript_sif poe.R "$CONFIG" status
column -t -s $'\t' "$COHORT_DIR/pilot_report.tsv" | less -S
```

Once `pilot complete` and QC is acceptable, submit the original scale pipeline
using the **filtered 195-child working manifests**:

```bash
bin/submit_pipeline.sh "$CONFIG" scale --submit
```

The two excluded children will not be scheduled for scale. Do not restore the full
manifests while that scale/model workflow is running. If later reintroducing them,
restore all five backup manifests together and reassess pilot, scale, and model
receipts; a 195-sample completed run does not establish completion for 197 samples.

Limitations: This patch adds a new opt-in command and does not change the generic
scheduler or the treatment of future scale-stage failures. If scale assignment
fails for additional samples, those failures require a new auditable exclusion
and collection step before models.
