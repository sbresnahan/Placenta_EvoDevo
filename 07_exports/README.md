# 07 — Human aggregate review

**Purpose:** prepare a private aggregate candidate, then write only the exact
locally approved aggregate table. These are explicit R actions, not LSF or
transfer jobs.

**Prerequisites:** completed human analysis/result-integrity review, approved
local disclosure policy and integer cell threshold. Keep analysis settings
fixed; use the separate release configuration described in the
[GUSTO walkthrough](../docs/WALKTHROUGH.md#human-gusto).

**Inputs:** eligible `paths.model_promoter`, public annotation identities,
`exports.*` approvals and the exact reviewed candidate hash.

From the script root, with `RELEASE_CONFIG` set to that private release copy:

```bash
bin/Rscript_sif 07_exports/prepare_candidate.R "$RELEASE_CONFIG"
```

Review the candidate, its exact SHA256 and manual/complementary suppression.
Only after completing the separate release approvals:

```bash
bin/Rscript_sif 07_exports/release_aggregates.R "$RELEASE_CONFIG"
```

**Outputs:** `paths.candidate`, then the approved `paths.release`.
Unknown columns, identities, free text, insufficient cell sizes or changed
candidate bytes block release. These actions validate candidate/schema/policy,
not the entire analysis receipt chain.

**Next:** no automatic transfer or publication. Keep participant tables,
logs, allocations, indices, receipts and candidates on approved private storage.
