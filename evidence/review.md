# Nutritive-tissue POE resources

The evidence catalog contains 26 study/cohort-stage entries: 14 maize and
12 mouse. These are not 26 independent datasets; rows can represent stages or
components of a project. The catalog is a bounded review, not an exhaustive
resource census. The four configured public profiles are the execution
selections; other catalog entries are not automatically included.

## Readiness definitions

Resource evaluation considers original raw-read accessions, RNA protocol,
tissue/stage, physical runs versus independent units, parental genomic evidence,
actual strain identity, compatible genome/SNP/annotation coordinates and
promoter identifiability. Primary methods, official ENA/GEO/BioProject records
and genotype resources support resolved entries. Unresolved candidates remain
explicit. Accessions, source URLs, caveats and evidence citations are recorded
in [dataset_evidence.csv](dataset_evidence.csv).

Read availability, marker compatibility, and promoter callability are evaluated
separately. Reference-strain genotypes describe the named reference stock.
Biological-unit metadata determine the independent-unit denominator.

## Maize

**B73/Ki11, 20 DAP:** promoter-focused configured profile with four biological
pools, each containing two technical runs. The deposited long-read-informed
annotation and B73v4 resources support reference preparation. [9] Structural
promoter groups become callable endpoints after raw-read and model eligibility
checks.

**B73/W22, 14 DAP, Anderson 2021:** balanced reciprocal option with three
ear-pool replicates per direction and PE125 RNA. PRJNA623806 contains 18 runs
across parental pairs; the configured six B73/W22 runs total 18.265 GB compressed
ENA FASTQs. Compatible SNP/promoter assets remain conditional; a concatenated
genome is not a ready Ornaments VCF. [49]

**Higgins 2025:** developmental extension at 11, 14, 17 and 21 DAP, with 23
libraries and one missing W22-maternal 21-DAP replicate. PRJNA1110585 describes
the time course; PRJNA1123772 describes the multigenotype component. B73v5
cannot directly reuse AGPv4 coordinates. Select genotypes explicitly and
deduplicate reused Waters/Anderson data. [545, 557, 558]

**Dong 2023:** 11-DAP B73/Mo17/CAU5 reciprocal design with three reported
replicates per direction. PRJNA765150 mixes assays/tissues. Public WGS
identities include SRR12415217 Mo17, SRR12415218 B73 and SRR3124079 CAU5.
A validated ready VCF and stock-compatibility review are still needed;
the project does not perform new WGS calling. [509, 521]

Waters 2011 pools hybrid replicates before sequencing; Waters 2013 needs a
run-to-ear audit and includes reused B73/Mo17 data. Xin 2013 early whole kernels
are not isolated endosperm. [511, 518, 550] Mu and ROS1ab experiments offer
perturbation contexts with sibling/introgression/EMS-background caveats.
[516, 551] The 2026 spatial/single-nucleus resource offers cell-type context,
but cells/spots are not independent ears and promoter/raw-resource readiness
is unresolved. [530]

## Mouse

**Babr/CAST, E7.5 ExE:** five configured WT SE libraries. Mean 250 nt and SD
50 nt are one assumed fragment setting, not a sensitivity grid. A compatible
Babr genome-wide variant source is unresolved; Babraham C57BL/6J ancestry
supports only an explicit B6J-reference proxy assumption, not an original-parent
genotype measurement. [506]

**FVB/NJ × CAST/EiJ, E12.5 placenta:** GSE75957 includes four female placentas,
two per cross direction, with PE100 stranded total RNA. Four E16.5 libraries
form a separate developmental option. Selected run metadata total 22.014 GB
at E12.5 and 46.555 GB across both stages. Litter IDs determine independent-unit
eligibility. The configured cohort contains only E12.5 WT reciprocal samples.
Maternal decidua/blood can contribute maternal expression. [481]

Joint MGP v5 GRCm38 supports FVB_NJ and CAST_EiJ calls with per-parent `FI=1`.
Another strain can cause a joint site-filter failure without invalidating two
high-confidence parental calls. Conversely, absence from a strain-specific
variant-only VCF is not a reference genotype. T2T assemblies are not direct
substitutes for GRCm38 calls. [546, 547, 555, 556]

The 2026 layer-resolved study PRJNA1020824 contains 146 RNA and 16 bisulfite
PE runs across crosses/tissues. Exact sample/animal/strain mapping requires a
supplemental crosswalk; layers/lanes from one placenta are dependent. Its
preprint is not another cohort. [476] Wang 2011 has one AKR/PWD placenta per
direction; 15 lanes cannot resolve zero residual df with n=2/rank=2, and
PWD is not PWK. [544] BDF1/CAST ExE requires B6/DBA-aware markers.
TSC cultures, Airn deletions and broader atlases are separate conditional
contexts. [29, 529, 503, 494]

## Human extensions

**ELGAN:** collaboration/access opportunity, not evidence of absent parental
data. The project record describes parent-offspring genotypes and placental
multiomics for ELGAN/GUSTO. Exact usable overlap and access need investigator
confirmation; project n=86/n=102 are not callable pair-family counts. [149]

**Queensland Family Cohort:** EGA dataset EGAD50000001850, study
EGAS50000001301, DAC EGAC50000000755; eight complete families with placental
PE150 RNA and parental WGS. This can support technical work but not the
project's ten-family inference floor. ONT is genomic DNA, not long-read
RNA. [136, 383]

**Mali:** dbGaP phs001782. Distinguish 111 whole-transcriptome biopsies from
59 placentas from the targeted-RNA component with 262 biopsies/251 placentas.
Targeted parental-DNA panels limit novel-promoter anchors; repeated biopsies
are not independent placentas. [144] Pilvar 2019 is conditional on verified
joint raw-RNA/parental-genotype access. [36, 384]

See [human_cohort_recommendations.csv](human_cohort_recommendations.csv) for
the structured recommendations. These extensions are not supplied GUSTO inputs.

## Interpretation

The workflow tests conditional parent-of-origin/promoter discordance, not a
conflict mechanism by itself. Maternal admixture, stock mismatch, cis effects,
developmental composition and quantifier allocation are alternative explanations.

Raw expression, genotype provenance and model eligibility are separate gates.
Public reciprocal models use the configured dosage and rank-based df;
one cross direction cannot separate strain from parent of origin.
Human inference requires typed anchors and the specified family/cis/HC3/
bootstrap conditions. Coordinate support is not direct joint-molecule evidence.
Human whole-gene analysis/masking and TE comparisons are outside scope.
