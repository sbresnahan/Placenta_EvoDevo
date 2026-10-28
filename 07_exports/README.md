# 07 — Human aggregate review

**Purpose:** prepare a private aggregate candidate, then write only the exact
locally approved aggregate table. These are explicit R actions, not LSF or
transfer jobs.

**Prerequisites:** completed human analysis/result-integrity review, approved
local disclosure policy and integer cell threshold. Keep analysis settings
fixed; use the separate release configuration described in the
[GUSTO recipe](../docs/Placenta_EvoDevo_step_by_step.md#part-6-human-gusto--private-trio-linked-placenta-analysis).

**Inputs:** eligible `paths.model_promoter`, public annotation identities,
`exports.*` approvals and the actual manually reviewed candidate.

From the script root, with `RELEASE_CONFIG` set to that private release copy:

```bash
bin/Rscript_sif 07_exports/prepare_candidate.R "$RELEASE_CONFIG"
```

Review the actual candidate and manual/complementary suppression.
Only after completing the separate release approvals:

```bash
bin/Rscript_sif 07_exports/release_aggregates.R "$RELEASE_CONFIG"
```

**Outputs:** `paths.candidate`, then the approved `paths.release`.
Unknown columns, identities, free text or insufficient cell sizes block release.
The operator must review the actual candidate before export; earlier reviewed
bytes are not compared automatically. These actions validate candidate/schema/policy,
not the entire analysis receipt chain.

**Next:** no automatic transfer or publication. Keep participant tables,
logs, allocations, indices, receipts and candidates on approved private storage.
