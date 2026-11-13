# Synchronize regulatory configuration for all six cohorts

The primary pipeline's `config/<cohort>.private.yaml` files predate the opt-in
regulatory analysis. `use_cohort` deliberately copies an example only when its
private file is missing, so old private configs do **not** inherit new blocks.

From the actual Seadragon checkout:

```bash
export SCRIPTS_DIR=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo/scripts
cd "$SCRIPTS_DIR"

# Show changes first. The tool does not rewrite the original settings.
python3 bin/sync_regulatory_configs.py --dry-run --create-missing

# Insert only missing regulatory keys in existing private configs; create
# missing private configs by copying their complete cohort-specific examples.
python3 bin/sync_regulatory_configs.py --create-missing
```

You can synchronize only one cohort:

```bash
python3 bin/sync_regulatory_configs.py fish_Pprolifica_AB
```

The updater is standard-library Python, idempotent, and checks the cohort ID,
species and schema. Existing values, including site paths, runtime settings,
customized resource paths, and credentials, are left unchanged. When an existing
private config needs changes, it is first backed up as
`config/<cohort>.private.yaml.bak.<timestamp>` and updated atomically.
Re-running the command does not overwrite the previous values. Default fish CPU
count is four; an explicit larger initial CPU value can be set with
`--fish-cpus 8`. The option only fills missing keys; manually customized counts
are never overwritten. The requested slots must be available in the LSF
allocation (`LSB_DJOB_NUMPROC >= regulatory.fish.cpus`).

The configuration sections are species-aware:

- All six cohorts: `bootstrap`, `seed`, `min_depth`, `max_pairs_per_gene`.
- `fish_Pprolifica_AB`: the documented RepeatModeler/RepeatMasker container
  settings, explicit `regulatory.fish.cpus: 4`, and a matching `lsf.stages.fish_repeats` (`long`,
  4 CPUs, 64 GB, 120 hours) starting allocations. Existing LSF values are
  retained. These are editable scheduling defaults, not measured full-genome
  resource requirements; check site queue policies and observed use.
- `GUSTO`: `human_context.require_qtl: true`. Supply the real private atlas
  QTL tables and phenotype crosswalk before running human-context; nothing is
  fabricated or downloaded by this configuration updater.

This changes no cohort data, receipts or primary pipeline output. It does **not**
submit jobs or trigger reference/retrieval/model reruns.

For submission with `bin/submit_pipeline.sh`, the fish LSF stage settings are
now populated automatically. Review these starting allocations against local
queue policy before launching full-genome discovery. The generator does not
replace any previously customized values. The human regulatory-context
analysis additionally requires reviewed atlas data and phenotype mapping.
