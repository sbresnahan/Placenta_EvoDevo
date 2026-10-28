# 03 — Genotype preparation

**Purpose:** prepare public parent-discriminating SNPs or private typed-child
and personalization inputs with consistent reference coordinates.

**Prerequisites:** reference catalog and approved matching sources. Public
marker-source paths/build/SHA256 must be complete. Human inputs need exact
trio linkage and independent RNA/DNA/genotype QC.

**Inputs and outputs:**

| Branch | Inputs | Outputs |
|---|---|---|
| Public | Panzea, joint MGP or prepared opposite-homozygote source | `paths.markers`, `paths.tx_vcf`, `paths.tx_marker_map`, marker receipt |
| Human | Typed trio VCF and separate child-personalization VCF | Per-library typed/personal query/VCF files, transcript mapping, child receipt |

These workers are automatically submitted by the reference graph:

```bash
bin/submit_pipeline.sh "$CONFIG" reference
```

This command is a dry-run. Use the [walkthrough](../docs/WALKTHROUGH.md) to
submit once; do not launch a second graph just to run this module.
`01_prepare_markers.lsf` handles public markers.
`02_prepare_child.lsf` uses host bcftools and R transcript/linkage checks.

Typed anchors orient human parents; imputed personalization variants do not
become typed anchors. Missing genotypes remain missing. A prepared asset must
actually exist; this stage does not perform WGS variant calling or imputation.

**Next:** the automatically dependent [Ornaments index](../04_ornaments/README.md).
