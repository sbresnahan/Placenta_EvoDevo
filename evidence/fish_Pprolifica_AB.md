# Poeciliopsis prolifica reciprocal-cross placenta: resource audit

Scientific source: Zhang et al. (2026), *Science Advances*,
"Independent origin of genomic imprinting in placental fish,"
https://doi.org/10.1126/sciadv.aeh5487.

## Verified dataset identifiers

| Asset | Source | Readiness |
|---|---|---|
| Placental, embryonic and maternal follicle RNA-seq | GSE318355 / SRP674561; https://www.ncbi.nlm.nih.gov/Traces/study?acc=SRP674561 | Public study; exact run-to-unit/direction selection requires GEO/SRA metadata review |
| Paper scripts and processed SNP/allelic count tables | https://doi.org/10.5061/dryad.ns1rn8q8r | Published 2026-09-30; DataS1 is placenta, S2 embryo, S3 maternal follicle |
| Line A reference FASTA and RefSeq annotation | GCF_027474105.1; https://pmc.ncbi.nlm.nih.gov/articles/PMC10655195/ | Reference assembly/annotation published; coordinate compatibility to be verified |
| Line A PacBio Revio HiFi whole-genome DNA | PRJNA1509661 / SRP725350 / SRX34719858 / SRR40078298 / SAMN62272763; https://www.ncbi.nlm.nih.gov/sra/SRX34719858 | One raw WGS run, ~12.3 Gb; no precalled VCF or aligned BAM identified in cited source listings |
| Line B independent genomic DNA | No matching genotype dataset verified | Missing; RNA-based fixed differences from the paper may supply conditional SNP sites |

The paper used Trimmomatic/STAR/TopHat alignment followed by GATK 3.6 RNA
SNP calling to identify fixed Line A versus Line B exonic variants. The new
PacBio Revio run was extracted from an **adult Line A male** (per NCBI
BioSample SAMN62272763); it cannot independently infer Line B genotypes.
Line A HiFi genotypes can be used to cross-check the published SNP list,
and kinetic methylation tags permit a separate epigenetic investigation.

The published study finds 84 imprinted protein-coding genes. Do not import that
number as an estimate from this pipeline. Our model is an independent,
exploratory reanalysis conditional on original RNA-based variant confidence,
sample/parent assignment, Ornaments read assignment, depth, stage and tissue.

## Conditional implementation

- Profile: `config/fish_Pprolifica_AB.example.yaml` (`public_cross`, `fish`).
- RNA manifest curation: `bin/prepare_fish_manifests.py`, output from official
  ENA `read_run` metadata and manually verified GEO cross direction.
- Target only the fetal **pericardial sac**, not mixed embryonic tissue or
  maternal follicle; separate stages into separate cohorts if necessary.
- Markers: `prepared` opposite-homozygous fixed SNPs with exact matching
  GCF_027474105.1 genome coordinate and allele orientation. No reference-only
  assumption for Line B. Supplementary RNA genotypes should remain labeled
  as RNA-derived even after independent Line A HiFi confirmation.
- Run the ordinary retrieval/reference/pilot/scale/models pipeline only after
  these checks. HiFi mapping/calling is presently an upstream input-preparation
  workflow, not a separately implemented Seadragon stage.
