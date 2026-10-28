# 00 — Preflight

**Purpose:** check the selected cohort settings/manifests and the actual
native/R environments without installing software, acquiring data or submitting
jobs.

**Prerequisites:** complete [shared setup](../docs/Placenta_EvoDevo_step_by_step.md#part-1-shared-setup--do-this-once),
select a profile, copy its operational manifests into the data directory and
edit its private YAML. Set `SCRIPTS_DIR`, `POE_SITE_FILE` and absolute `CONFIG`.

**Inputs:** private YAML, trusted site file, manifest tables, installed tools,
R image/library and existing native environments.

From the script root:

```bash
00_setup/preflight.sh "$CONFIG"
00_setup/preflight.sh "$CONFIG" --execute
```

**Outputs:** metadata validation; accepted `paths.native_inventory` on execution
success, or a printed pending inventory path to inspect if environment checks fail.
Acceptance requires available native commands and compatible R packages.
See [installed-tool identity](../docs/SETUP.md#installed-tool-identity).

**Next:** public [retrieval](../01_retrieval/README.md), or private GUSTO input/
reference preparation. Preflight alone does not certify all files, site policy
or biological callability.
