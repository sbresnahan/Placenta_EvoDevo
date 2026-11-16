# H6/H7 implementation and aggregate archive

Start with `docs/Seadragon_analysis_runbook.md`, `docs/H6_H7_MODELS.md`
and `docs/rerun_exposure_nb.md`. Output definitions are in
`docs/exposure_output_schema.md`.

- H6: continuous matching-family TE-RNA × allele NB2 models.
- H7: GUSTO maternal PFAS → TE-RNA → signed POE mediation; modification optional.
- Targets use the active native promoter window with insertion-level provenance.
- Explicit private inputs, reviewed clusters and measured pilots are required.
- H1 and H5 significance do not select mediation jobs.

Required inputs: private count batches, sample crosswalk, TE locus/count/QC
tables, GUSTO PFAS dictionary, frozen functional memberships and H5 context.
Templates contain schemas and require actual reviewed inputs before execution.

Export with `poe_exposure.R CONFIG export`, then
`bin/export_results.R --exposure-nb --out FILE.tar.gz CONFIG...` through the
site R wrapper. Inventories record filenames and sizes. Maintain immutable
inputs and fresh output directories for each preparation.
