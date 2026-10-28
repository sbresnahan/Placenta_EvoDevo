# Software checks

This update uses R/shell syntax checks, static configuration/schema checks,
documentation-link checks, and sequential edit applicability checks.
Regression suites and cohort analyses were not run for this update.

The supplied tests cover parental orientation, coordinate mapping, acquisition,
completion records, model eligibility, host/R separation, scheduler interfaces,
and aggregate validation. Host-integration fixtures use mocked Ornaments,
bcftools, personalizer, Singularity, and LSF commands.

Runtime validation includes actual file access, MD5/byte and FASTQ integrity,
reference/genotype compatibility, count schemas, pilot strand diagnostics,
model eligibility, and aggregate-column/cell thresholds.

The [walkthrough](Placenta_EvoDevo_step_by_step.md) provides the execution order.
