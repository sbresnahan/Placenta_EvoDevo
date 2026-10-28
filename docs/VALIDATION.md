# Validation and site checks

## Local coverage

The supplied software checks cover Mendelian orientation, exon/strand mapping,
missing evidence, metadata, acquisition, model gates, completion records,
host/R separation, scheduler interfaces, and privacy allowlists.
The completion-only update is limited to syntax and static consistency checks;
no regression suite or cohort analysis was run for that update.
Earlier assertion counts do not certify the changed interfaces.

**Ornaments, bcftools, the Python personalizer, Singularity and LSF are mocked
in host-integration tests.** Counts and resource values are deterministic
software fixtures, not biological findings or measured site performance.
The small human fixture deliberately fails the inference family-count floor.
Separate algebra tests cover eligible and blocked human models.

Scheduler tests use a fake `bsub` to check real numeric-ID substitution,
dependencies, malformed responses, failures and replay prevention. R does not
call host `bsub` directly.

Completion status records successful stages and required nonempty outputs.
It does not determine whether earlier results are still current. The operator
must rerun affected stages and their dependents after changes.

## Verify on Seadragon before full runs

- Actual modules/conda tools, command availability and expected output schema.
- R image execution, library/package compatibility and filesystem binds.
- Scheduler permissions, queue policy, memory limit/reservation semantics and
  environment propagation.
- New script/data root availability, POSIX behavior, write access and total
  storage capacity.
- Intended shared human genome and annotation availability/builds at their
  configured source locations.
- Approved internet host and actual public reference/marker/read contracts.
- Independent protocol, biological units, strain/proxy qualification and raw
  pilot callability.
- Complete private GUSTO inventory/linkage/genotype/covariate review.
- Measured index/quantification memory/time and total cohort storage budget.
- Local human disclosure/suppression policy and exact aggregate review.

Local testing does not establish a successful real-cohort Ornaments run,
actual Seadragon compatibility or any biological result. The 8-GB-class example
requests are not measured GUSTO/alternative-index requirements. Process RSS and
elapsed time do not measure total workflow cost; endpoint work size does not
measure disk high-water usage.

Human whole-gene analysis/masking and TE comparisons are outside implemented
scope. Technical success and structural promoter support are distinct from
inferential eligibility and evidence for imprinting/conflict.
