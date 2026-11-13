# Fish force exclusion of one failed pilot

This is an explicit, fish-only operator override for the low-pseudoalignment
`poeProA6xB3_pla8_stage45` sample. It never fabricates pilot receipts.

The active `sources`, `libraries`, and `units` files are reduced from six to five;
backups and exclusion provenance are stored under `manifests/fish_force_exclusion_*`.
The private YAML policy is changed from `quantification_only` to
`force_pooled_unadjusted`. Original stage40/stage45 labels remain in the units
file. The existing reference and other five pilot files are not modified.

The resulting design has A-mother n=2 (one stage40, one stage45) and B-mother
n=3 (all stage40). The two factors are confounded in this tiny cohort. A pooled,
stage-unadjusted model is explicitly exploratory. Litter/cluster identities remain
unknown; nominal p-values should not be treated as confirmatory.

## Before applying

Confirm no relevant pilot/collector job is running. Cancel the old pending
collector (LSF job 81176532 from the original six-sample graph), if still PEND:

    bkill 81176532

## On Seadragon

    cd "$SCRIPTS_DIR"

    python3 bin/force_drop_fish_failed_pilot.py \
      --cohort-dir "$COHORT_DIR" --config "$CONFIG" \
      --exclude poeProA6xB3_pla8_stage45 --preview

Review the preview and, only if it identifies the correct one failed sample:

    python3 bin/force_drop_fish_failed_pilot.py \
      --cohort-dir "$COHORT_DIR" --config "$CONFIG" \
      --exclude poeProA6xB3_pla8_stage45 --apply

Reconcile new active manifest and receipts, then schedule only incomplete
pilot libraries and collector:

    bin/submit_pipeline.sh "$CONFIG" pilot --resume
    bin/submit_pipeline.sh "$CONFIG" pilot --resume --submit

The resume graph schedules no repeated successful pilot workers. After the
collector completes, inspect pilot diagnostics, then submit scale normally.

The forced models retain a `mixed_stage40_stage45_unadjusted` label and an
explicit `exploratory_stage_unadjusted_cluster_independence_unverified_small_N`
uncertainty label. The existing pooled fit is not stage-adjusted or replicated
at litter level. Do not assume that force makes its P-values valid.
