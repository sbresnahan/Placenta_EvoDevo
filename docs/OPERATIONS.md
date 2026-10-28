# Operations and troubleshooting

Use [the step-by-step workflow](Placenta_EvoDevo_step_by_step.md) for the normal run order. All examples
assume the script root is your current directory and `CONFIG` names one
absolute private YAML.

## Inspect jobs and completion

```bash
bjobs
bjobs -a
bin/Rscript_sif poe.R "$CONFIG" status
```

For one scheduler job:

```bash
read -r -p "Numeric LSF job ID: " JOB_ID
[[ "$JOB_ID" =~ ^[0-9]+$ ]] || { printf 'Numeric job ID required\n' >&2; exit 1; }
bjobs -l "$JOB_ID"
```

Use site accounting/log tools if the job is no longer in recent LSF history.
Do not infer success from its absence. `status` validates reference, pilot,
scale authorization, complete quantification and models; it does not query LSF.
It writes `paths.status` with complete/authorized/blocked and reason fields.
Human reference status covers the shared reference; inspect every child/index
job and its logs as well.

The driver writes cohort logs under `paths.logs/lsf/` with numeric `%J` IDs.
Tool/resource logs use `templates.quant_log` and `templates.resource_log`.
Planning registries derive from `paths.job_registry` with unique suffixes:

- `.tsv`: dry-run/host-planned graph. `DRYRUN_` and `HOST_PLANNED_` are not jobs.
- `.host.sh`: generated single-use host submission plan.
- `.submitted.tsv`: actual accepted numeric IDs recorded as submission proceeds.
- `.host.sh.started`: replay guard, not a completion receipt.

Keep these files private for GUSTO. A `.started` marker only says execution
began; check scheduler records, submitted IDs and output receipts.

## Host submission interface

Use `bin/submit_pipeline.sh` for normal stages. R prepares quoted job commands;
the shell executes host `bsub`. Calling `run_pipeline.R` or `poe.R` with
`--submit` alone prepares a host plan but does not submit it.

The driver uses ordinary per-library jobs, not arrays. Reference/marker/index
dependencies and assignment/collection dependencies use numeric `done(jobID)`.
No long-running monitor or automatic polling loop is required.

An external numeric dependency is optional:

```bash
read -r -p "Prerequisite numeric LSF job ID: " AFTER_JOB
[[ "$AFTER_JOB" =~ ^[0-9]+$ ]] || { printf 'Numeric job ID required\n' >&2; exit 1; }
bin/submit_pipeline.sh "$CONFIG" reference --after "$AFTER_JOB"
```

This is a dry-run; `--submit` is still required for actual submission. A
dependency does not waive input, receipt or approval checks. Do not chain pilot
completion to scale; inspect the report, then explicitly submit scale.

## Manual submission

Prefer the driver so dependencies, resources and logs remain consistent.
For a deliberately manual single-library pilot, after successful reference
preparation and site review:

```bash
export SCRIPTS_DIR="$PWD"
export POE_SITE_FILE="$PWD/config/seadragon.sh"
mkdir -p /rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/data/logs
read -r -p "Exact primary-library sample_id from the manifest: " POE_SAMPLE
export POE_SAMPLE
[[ "$CONFIG" == /* && -f "$CONFIG" ]] || { printf 'Absolute CONFIG required\n' >&2; exit 1; }
bsub -env "all,CONFIG=\"$CONFIG\",SCRIPTS_DIR=\"$SCRIPTS_DIR\",POE_SITE_FILE=\"$POE_SITE_FILE\",POE_SAMPLE=\"$POE_SAMPLE\"" \
  < 04_ornaments/02_run_pilot.lsf
```

This submits a real job. Do not also submit the same work through the driver.
It does not run the all-library collector or authorize scale. Only use exact
manifest library IDs. Shared workers do not need a library selection.

Checked-in headers default to the project's `data/logs/` directory for manual
jobs. The driver overrides them with cohort-specific logs. Create the header
directory first or pass existing `-o`/`-e` destinations. If needed, override
`-q`, `-n`, `-M`, `-R` and `-W` with approved values; headers are templates,
not automatic interpretations of edited YAML resources.

Set `CONFIG`, `SCRIPTS_DIR` and `POE_SITE_FILE` explicitly for stdin submission;
LSF spool paths do not locate the project. Driver `-env` values cannot contain
commas, quotes or newlines; spaces are supported.

`00_run_pipeline.sh` is an optional short driver job controlled by `STAGE`,
`SUBMIT=yes`, optional `AFTER`. It exits after graph submission.
The ordinary interactive shell interface is simpler.

## Restart and failure handling

| Problem | Action |
|---|---|
| Missing native tool | Inspect the configured environment and command paths; complete [setup review](SETUP.md#installed-tool-identity) |
| Resource missing checksum/size | Complete the intended asset contract; do not replace it with an arbitrary file |
| Missing reference/marker approval | Verify build, exact bytes and parent/proxy provenance first |
| Failed download | Inspect response/error; retry through the supported retrieval stage after resolving the cause |
| Native tool failure or empty allele output | Inspect private tool/LSF logs; a success receipt must not be fabricated |
| One missing pilot library | Resolve it and rerun needed work; collector requires all libraries |
| Resource budget too small | Revise reviewed settings; rerun affected prerequisites/pilot before accepting a new report |
| Missing/incomplete receipt | Resolve the failure or missing output, then rerun the stage; never fabricate a success record |
| Lock collision | Inspect recorded owner, live jobs and processes; remove only after verifying the lock is stale |
| Partial scheduler submission | Reconcile actual IDs in `.submitted.tsv` and LSF before any retry |
| Technical success but ineligible model | Report design/callability limitation; do not lower inference gates |

Submission can fail after earlier jobs have already been accepted. Missing or
malformed job IDs must be reconciled with LSF, not assumed to mean no submission.
Do not replay a `.host.sh` plan or erase its guard to force replay. Generate a
fresh plan only after accounting for outstanding jobs and output ownership.

Locks protect shared reference and per-library mutations, plus collectors/
models. No automatic stale-lock deletion is provided. Do not rebuild a reference
while downstream jobs are reading it.

Completion records indicate that work finished and required nonempty outputs exist; they do not detect changed inputs, settings, software, or output bytes. An explicit stage rerun performs the requested work rather than automatically reusing prior quantification. A stage records that a rerun has started, so failure cannot leave that stage's earlier success in effect.

Before a deliberate rerun, reconcile active jobs and move aside affected outputs/records if existing-file protections require it. Rerun affected downstream stages yourself; do not delete active locks or overwrite files another job is reading.

Download resume accepts only verified matching HTTP range responses; an
HTTP 200 response is not appended to a partial. Retries/timeouts are bounded.
The supported transport is ENA HTTPS; no automatic SRA conversion or trimming
fallback is supplied.

## Keep human analysis and release separate

All participant tables, command plans, paths, logs, counts, snapshots, receipts
and indices stay private; entry points use `umask 077`. Aggregate review and
release are explicit local actions, with no automatic transfer.

Use the walkthrough's separate release configuration after recording completed
analysis/result integrity. Do not restart analysis with the release configuration. Aggregate validation does not
replace model-receipt validation or local disclosure-policy review.
