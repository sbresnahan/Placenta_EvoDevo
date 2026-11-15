# P. prolifica Line A/B marker reconstruction close to the published placental panel

**Reference:** Zhang et al., *Independent origin of genomic imprinting in placental fish*, *Science Advances* (2026), DOI `10.1126/sciadv.aeh5487`; Dryad DOI `10.5061/dryad.ns1rn8q8r`.

## Target and limitations

The authors report **55,224 fixed-difference placental SNPs** among **193,524 placental SNPs**. Their analysis used joint RNA-seq variant calling and within-line crosses, with a low-depth cutoff of six RNA reads per tissue. The deposited Dryad CSV tables contain allele *read counts*, not the exact original VCF genotypes or the inclusion flag for the published 55,224 SNPs. Accordingly, exact site-wise reproduction cannot be established from the tables alone.

We reconstruct **53,734 RNA-supported candidate SNPs**, a numerical difference of **1,490 sites (2.70%)** from the reported benchmark. This is **count similarity, not demonstrated SNP-set overlap**. We did not tune the depth cutoff to force the published count.

### What the reconstruction actually tests

- Start with `DataS1.csv` (placental SNP coordinates; all consequence classes).
- Pool REF and ALT RNA counts **within the same named biological fish** across placenta (`DataS1`), embryo (`DataS2`), and maternal follicle (`DataS3`), matching `SNP_ID`, scaffold, exact-position string, REF and ALT. Do not count the same fish in two tissues as two independent genotype replicates.
- Controls: Line A individuals `A6xA6_1`, `A6xA6_2`, `A6xA6_7`; Line B individuals `B3xB3_4`, `B3xB3_11`. Only these within-line controls inform line alleles; no hybrid F1 expression is used.
- A fish is callable at a site with `REF + ALT >= 6` pooled reads. An observed Line A control must have **only REF reads**; an observed Line B control must have **only ALT reads**. Require at least one callable individual from each line and no contradictory callable controls from either line.
- The supplied publication specifies A=REF, B=ALT, so candidates have `focal_index=0` (A), `alternate_index=1` (B). Sites supporting reversed orientation are **excluded**, never auto-flipped.
- Exclude invalid/negative parental allele counts, cross-tissue ID/coordinate/allele disagreements, ambiguous or imprecisely printed genomic positions and invalid REF/ALT alleles.

### Results on deposited tables

| Class | SNPs |
|---|---:|
| Replicated in ≥2 distinct Line A and ≥2 distinct Line B individuals | 27,739 |
| Replicated in ≥2 individuals of only one line | 19,024 |
| Evidence from exactly one individual in each line | 6,971 |
| **Total reconstructed** | **53,734** |
| Article benchmark (not site-wise verified) | 55,224 |

These are RNA-derived **proxy genotypes**. Monoallelic RNA can result from expression or mapping effects, and pooled same-fish tissues do not constitute independent germline genotype validation. Retain evidence tiers for sensitivity analyses.

## Reproduce on Seadragon

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export SCRIPTS_DIR="$PROJECT/scripts"
export COHORT_DIR="$PROJECT/data/fish_Pprolifica_AB"
DRYAD="$COHORT_DIR/reference/doi_10_5061_dryad_ns1rn8q8r__v20260930"
OUT="$COHORT_DIR/reference/dryad_lineAB_published_like_v1"
cd "$SCRIPTS_DIR"
python3 bin/prepare_fish_published_like_markers.py \
  --datas1 "$DRYAD/DataS1.csv" \
  --datas2 "$DRYAD/DataS2.csv" \
  --datas3 "$DRYAD/DataS3.csv" \
  --out-dir "$OUT"
cat "$OUT/lineA_lineB.published_like_report.json"
```

This first run creates **candidate** `lineA_lineB.published_like_candidates.tsv`, annotated with evidence tiers; it does **not** produce the six-column active pipeline marker file.

## Validate the exact reference assembly before promotion

The Dryad tables use `SCAFFOLD...` labels, which are not automatically identical to the `GCF_027474105.1` FASTA contigs. Determine a verified alias map or liftover first. Simple name ordering is not evidence of coordinate identity.

```bash
GENOME="$COHORT_DIR/reference/genome.fa"   # verify actual configured assembly/filename
samtools faidx "$GENOME"
VALID="$COHORT_DIR/reference/dryad_lineAB_published_like_genome_checked_v1"
python3 bin/prepare_fish_published_like_markers.py \
  --datas1 "$DRYAD/DataS1.csv" \
  --datas2 "$DRYAD/DataS2.csv" \
  --datas3 "$DRYAD/DataS3.csv" \
  --out-dir "$VALID" --genome "$GENOME"
# If and only if a verified two-column alias map is needed, add:
# --contig-map "$COHORT_DIR/reference/verified_contig_aliases.tsv"
```

By default, ≥99% of candidates must resolve to the FASTA with matching reference bases; otherwise the script refuses to write `lineA_lineB.published_like.markers.tsv`. The genomic audit details each mismatch.

Before setting the resulting six-column marker file as the active cohort source, complete the GEO-to-SRA biological sample crosswalk, verify same-coordinate genome/GTF compatibility and assess marker callability in the six reciprocal F1 placentas. The SRA sample labels are not reconstructed from SNP imbalance.
