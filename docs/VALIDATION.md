# Validation and site checks

## Local coverage

The local suite contains **276 assertions**:

| Suite | Assertions | Scope |
|---|---:|---|
| `test_core.R` | 47 | Mendelian orientation, exon/strand mapping, conservation, missing evidence and keys |
| `test_workflow.R` | 19 | Receipts, linkage, protocol/identity, HC3 algebra and privacy allowlists |
| `test_v2.R` | 93 | Configurations, HTTP resume, checksums, MGP calls, promoter/sentinel selection, design/model and authorization gates |
| `test_cli.R` | 8 | Relocation, quoting, stage blockers and CLI behavior |
| `test_seadragon.R` | 109 | Host/R separation, activations, native identity, fan-out, collectors, locks and mocked integration |

Run in a disposable local development checkout with the required R packages
and system utilities. First stage the supplied manifest templates in isolated
cohort `manifests/` directories and point each example YAML's manifest entries
there. Keep supplied GUSTO schemas empty for the negative tests; never point
these tests at populated private cohort tables. Leave other deployment paths
as configured so the CLI suite can test its own output relocation.

```bash
Rscript tests/run_all.R
```

This is a software test, not the command to analyze a cohort. Package acceptance
uses this isolated manifest-path setup without changing tests or runtime code.
A tiny coordinate
fixture uses real local samtools when installed and real R FASTA/GTF processing.
The fixture is coordinate truth, not an experimental observation.

**Ornaments, bcftools, the Python personalizer, Singularity and LSF are mocked
in host-integration tests.** Counts and resource values are deterministic
software fixtures, not biological findings or measured site performance.
The small human fixture deliberately fails the inference family-count floor.
Separate algebra tests cover eligible and blocked human models.

Scheduler tests use a fake `bsub` to check real numeric-ID substitution,
dependencies, malformed responses, failures and replay prevention. R does not
call host `bsub` directly.

Documentation/package checks cover Markdown links, shell-block syntax, path
contracts, unchanged non-path scientific settings and archive file identity/
permissions after extraction into a path containing spaces. Operational
manifest copies are checked in a temporary local deployment layout.

## Verify on Seadragon before full runs

- Actual modules/conda tools, installed source/build correspondence and output
  schema.
- R image execution, library/package compatibility and filesystem binds.
- Scheduler permissions, queue policy, memory limit/reservation semantics and
  environment propagation.
- New script/data root availability, POSIX behavior, write access and total
  storage capacity.
- Exact shared human genome and annotation availability/hashes at their
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
measure disk high-water usage. Repeated image hashing adds unbenchmarked site
I/O overhead.

Human whole-gene analysis/masking and TE comparisons are outside implemented
scope. Technical success and structural promoter support are distinct from
inferential eligibility and evidence for imprinting/conflict.
