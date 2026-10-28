# Running and monitoring

## Submit a stage

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve --submit
bin/submit_pipeline.sh "$CONFIG" reference --submit
bin/submit_pipeline.sh "$CONFIG" pilot --submit
bin/submit_pipeline.sh "$CONFIG" scale --submit
bin/submit_pipeline.sh "$CONFIG" models --submit
```

Run each command when the preceding stage has finished. Human analyses begin
with `reference`. Omit `--submit` to print a plan.

## Inspect jobs and results

```bash
bjobs -a
ls -lt "$COHORT_DIR/logs/lsf/"
ls -lt "$COHORT_DIR/private/"jobs_*.submitted.tsv
bin/Rscript_sif poe.R "$CONFIG" status
```

The terminal prints submitted job IDs. Submission records contain graph nodes,
dependencies, and scheduler IDs. Logs are named by cohort, action, and job ID.

```bash
read -r -p "Job ID: " JOB_ID
bjobs -l "$JOB_ID"
read -r -p "Log filename: " JOB_LOG
less -S "$JOB_LOG"
```

| Status | Meaning |
|---|---|
| `PEND` | Queued |
| `RUN` | Running |
| `DONE` | Scheduler success |
| `EXIT` | Job failed; inspect its logs |
| `complete` | Stage receipt and required outputs are present |
| `requested` | Scale was explicitly requested |
| `blocked` | A prerequisite or required output is missing/incomplete |

## Rerun after an error

Read the failed job's logs, correct the reported input or setting, and submit
the affected stage again. Check already submitted jobs before repeating a
partially submitted graph. For deliberate reruns, move aside existing outputs
when a tool's existing-file protection requires it, after the jobs using them
have finished.

Download retrieval resumes validated partial files. Locks record their owner
under `receipts/locks/`; clear an abandoned lock after its owning process/job
has ended. Input or setting changes require rerunning affected stages and
their dependents.

Pilot reports retain resource usage and projections as informative outputs.
Queue, cores, memory and wall time remain explicit job settings.
