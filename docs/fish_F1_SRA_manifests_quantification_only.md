# Fish placental SRA manifests: six-F1 quantification-only activation

## Scope

The 53,734 reconstructed Dryad SNP markers have 100% REF agreement against
`GCF_027474105.1`, and 11 GEO aliases were matched exactly to 11 SRA runs.
This procedure stages six reciprocal F1 placentas for reference construction,
pilot, and quantification. Five within-line controls are retained in the
versioned audit, not analyzed as reciprocal F1s.

**Do not run the existing pooled public-cross POE models on these samples.**
The six F1s contain one A-mother stage40, two A-mother stage45, and three
B-mother stage40 placentas. The old orientation-only model is not adjusted
for stage. The new `fish_stage_policy: quantification_only` option allows
reference/pilot/scale, but blocks primary and additional POE models.

The NCBI SRA provenance URLs in the source manifest are accession landing
pages, **not** ENA FASTQ links. In SRA mode the pipeline acquires or reuses
reads by run accession; the manifest's mate byte sizes are obtained from the
existing `.sra.tsv` receipts.

Biological cluster/litter IDs are unknown and are recorded as missing rather
than inferred to be independent. `strand_mode=unstranded` is provisional and
must be checked by the pipeline's pilot strand diagnostics.

## On Seadragon after pulling the patched repo

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export SCRIPTS_DIR="$PROJECT/scripts"
export COHORT_DIR="$PROJECT/data/fish_Pprolifica_AB"
export CONFIG="$SCRIPTS_DIR/config/fish_Pprolifica_AB.private.yaml"
export POE_SITE_FILE="$SCRIPTS_DIR/config/seadragon.sh"
export CROSS="$COHORT_DIR/manifests/fish_geo_crosswalk_v1/fish_placenta_geo_dryad_crosswalk.tsv"
export STAGED="$COHORT_DIR/manifests/fish_sra_f1_v1"
cd "$SCRIPTS_DIR"
umask 077

python3 bin/prepare_fish_sra_manifests.py \
  --crosswalk "$CROSS" \
  --runinfo manifests/fish_Pprolifica_AB_SraRunTable.csv \
  --raw-dir "$COHORT_DIR/raw" \
  --strand-mode unstranded \
  --out-dir "$STAGED"

cat "$STAGED/fish_sra_manifest_report.json"
column -t -s $'\t' "$STAGED/fish_Pprolifica_AB_units.tsv"
```

The script reads all 11 SRA receipts; it does not download anything. It
aborts if either existing FASTQ does not match its receipt's path/size.
It refuses to overwrite versioned outputs.

## Activate reviewed manifests

Before activation verify the six-fish stage structure (A mother: one stage40,
two stage45; B mother: three stage40), the `run_id` to sample mapping, and
that the default active fish source/library/unit TSVs have not been filled
since the last check. This procedure backs up all three originals:

```bash
for key in sources libraries units; do
  active="$COHORT_DIR/manifests/fish_Pprolifica_AB_${key}.tsv"
  test -f "$active" || { echo "Missing $active" >&2; exit 1; }
  test "$(wc -l < "$active")" -eq 1 || {
    echo "Refusing to overwrite populated $active" >&2; exit 1;
  }
done

for key in sources libraries units; do
  active="$COHORT_DIR/manifests/fish_Pprolifica_AB_${key}.tsv"
  staged="$STAGED/fish_Pprolifica_AB_${key}.tsv"
  cp -p -- "$active" "${active}.pre_fish_F1_v1.bak"
  cp -- "$staged" "${active}.pending"
  mv -- "${active}.pending" "$active"
done

# This setting affects only the fish cohort and authorizes quantification,
# NOT a pooled POE model. Leave all original stage40/stage45 labels intact.
cp -p "$CONFIG" "${CONFIG}.pre_fish_stage_policy.bak"
if ! grep -q '^fish_stage_policy:' "$CONFIG"; then
  printf '\nfish_stage_policy: quantification_only\n' >> "$CONFIG"
fi

bin/Rscript_sif poe.R "$CONFIG" status
bin/submit_pipeline.sh "$CONFIG" reference
```

If the reference preview succeeds and `tests/test_inputs.sh "$CONFIG"` passes:

```bash
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

Check LSF receipts and reference outputs before proceeding to pilot:

```bash
bjobs -u "$USER" -w
bin/Rscript_sif poe.R "$CONFIG" status
```

Do **not** submit `models` or `additional`; the guard should reject them
until a separate, reviewed stage-aware analysis is implemented. Stage40 has
only one A-mother F1, so even a stage-adjusted model has very limited power.
