# 03 — Genotype preparation

Included in the `reference` graph.

| Profile | Inputs | Outputs |
|---|---|---|
| Public | Panzea, joint MGP, or prepared opposite-homozygote markers | Genomic markers, transcript VCF, coordinate map |
| Human | Typed trio VCF, child-personalization VCF, linked pedigree | Per-child queries, VCFs, transcript mappings |

Typed anchors supply maternal/paternal orientation. Missing genotypes remain
missing. Source/build, genotype, and REF checks run during preparation.

Next: the dependent [Ornaments index](../04_ornaments/README.md).
