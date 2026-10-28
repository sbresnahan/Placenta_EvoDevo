# 00 — Setup

Set `SCRIPTS_DIR`, `POE_SITE_FILE`, and the absolute `CONFIG` filename using
the [walkthrough](../docs/Placenta_EvoDevo_step_by_step.md#part-1-shared-setup).

```bash
00_setup/preflight.sh "$CONFIG"
```

Checks native commands, scheduler availability, the R image/library, and
required R packages. `--execute` performs the same check.

Expected output: `Required tools and R packages are accessible.`
The inventory records command paths and versions.

Next: [public retrieval](../01_retrieval/README.md), or GUSTO input preparation.
