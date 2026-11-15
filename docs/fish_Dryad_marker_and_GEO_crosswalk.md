# Fish Line A/B marker reconstruction and GEO crosswalk

**Source:** Zhang et al., *Independent origin of genomic imprinting in placental fish*, Science Advances (2026), DOI `10.1126/sciadv.aeh5487`. Dryad: `10.5061/dryad.ns1rn8q8r`, `DataS1.csv` and `DataS2.csv`. RNA runs: GEO `GSE318355`, SRA `SRP674561`. Genome specified by this repository: `GCF_027474105.1`.

## Scientific restrictions

- DataS1 has 11 placental samples: 3 A6×A6 controls, 2 B3×B3 controls, 3 A6-mother×B3-father F1, and 3 B3-mother×A6-father reciprocal F1. Only the *five parental-control samples* inform the marker allele assignment. The F1 expression fractions are **not** used to define parental genotypes.
- The inferred marker alleles are **RNA-supported proxies for line differences**, not sequenced germline genotypes. The embryo table (DataS2) provides an additional tissue concordance check using the same named controls; embryos are often paired with their placentas, so it is **not** independent biological genotype replication. Allelic bias, mapping bias, residual segregating variation, and differences between line individuals remain caveats.
- Dryad genomic positions use `SCAFFOLD...` contig names. Never assume these are identical to the FASTA identifiers of `GCF_027474105.1`; use a verified assembly/contig mapping. Do not convert imprecise printed positions to integers.
- Even after `--genome` reference checking, the final `fixed_exonic.markers.tsv` is a conservative *coding-effect-annotated RNA-derived proxy*. The downstream transcript GTF and marker/genome checks must also succeed. Do not claim verified germline origin from it.
- `cluster_id` and the biological-independence of multiple embryos from the same named cross are **not established**. The crosswalk is always `include=NO` until reviewed.

## Prepare Dryad markers (run on host, Python 3.6+)

```bash
export PROJECT=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_EvoDevo
export SCRIPTS_DIR="$PROJECT/scripts"
export COHORT_DIR="$PROJECT/data/fish_Pprolifica_AB"
cd "$SCRIPTS_DIR"

DRYAD="$COHORT_DIR/reference/doi_10_5061_dryad_ns1rn8q8r__v20260930"
OUT="$COHORT_DIR/reference/dryad_lineAB_candidates_v1"
python3 bin/prepare_fish_dryad_markers.py \
  --datas1 "$DRYAD/DataS1.csv" --datas2 "$DRYAD/DataS2.csv" \
  --out-dir "$OUT"
cat "$OUT/lineA_lineB.marker_report.json"
```

Outputs (no production marker is written without FASTA):

- `lineA_lineB.RNA_controls_candidates.tsv`: opposite near-homozygous A/B alleles in parental-control RNA, any SNP annotation.
- `lineA_lineB.RNA_controls_coding_candidates.tsv`: coding-consequence-annotated subset only.
- `lineA_lineB.RNA_controls_embryo_corroborated.tsv`: candidates with supporting A- and B-line embryo control RNA.
- `lineA_lineB.RNA_controls_coding_embryo_corroborated.tsv`: coding subset with embryo corroboration.
- `lineA_lineB.marker_report.json`: thresholds and input exclusions.

Defaults: per-control total ref+alt depth >=20, ALT fraction <=0.10 or >=0.90, >=2 callable A controls and both B controls, **no conflicting callable control**. All F1 samples are ignored for line-genotype inference. A site is embryo corroborated if at least one A and one B embryo control has >=20 reads, both agree on the alleles, and no other callable embryo control contradicts them. This is a separate-tissue QC tier, not an additional independent embryo sample. Invalid negative counts, non-SNP alleles, and imprecise/noninteger SNP positions are excluded. Duplicate genomic positions fail closed.

### Production-marker promotion requires the exact reference FASTA

Verify `genome.fa` is the correct assembly and index it on the host:

```bash
GENOME="$COHORT_DIR/reference/genome.fa"
samtools faidx "$GENOME"
head -5 "${GENOME}.fai"
head -5 "$OUT/lineA_lineB.RNA_controls_coding_candidates.tsv"
```

If the Dryad `SCAFFOLD...` identifiers **do not** match FASTA contig identifiers, create a two-column TSV of verified equivalent genomic coordinates with header:

```
dryad_scaffold\tfasta_contig
SCAFFOLD001\t<verified_reference_contig_ID>
```

A name resemblance or name ordering is **not** sufficient to infer equivalence. If the assemblies differ in coordinates, a simple alias map is invalid and a verified liftover/sequence alignment is necessary.

Run validation in a *new* output directory; if an alias map is required include `--contig-map PATH`:

```bash
python3 bin/prepare_fish_dryad_markers.py \
  --datas1 "$DRYAD/DataS1.csv" \
  --datas2 "$DRYAD/DataS2.csv" \
  --out-dir "$COHORT_DIR/reference/dryad_lineAB_genome_checked_v1" \
  --genome "$GENOME"
```

The script requires at least 99% of embryo-corroborated coding candidates to map to the FASTA with the correct REF nucleotide. It writes `lineA_lineB.genome_validation.tsv` for inspection. It creates `lineA_lineB.fixed_exonic.markers.tsv` **only on passing**. An existing production marker file is never overwritten. Change output directory for new analyses. Only *after* validation, update the cohort marker-source manifest or copy this versioned file to the configured `reference/lineA_lineB.fixed_exonic.markers.tsv`. Document the RNA-derived proxy explicitly.

## GEO-to-SRA crosswalk

The packaged NCBI RunInfo establishes `SRR` ↔ `GSM` and tissue, but **does not establish** the 11 `GSM` ↔ Dryad sample labels. Do not guess correspondence from numerical ordering. Obtain official GEO family SOFT on a machine/transfer node with network access:

```bash
mkdir -p "$COHORT_DIR/reference/geo"
curl --fail --location --retry 3 \
  'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE318nnn/GSE318355/soft/GSE318355_family.soft.gz' \
  -o "$COHORT_DIR/reference/geo/GSE318355_family.soft.gz"
```

Then:

```bash
python3 bin/build_fish_geo_crosswalk.py \
  --runinfo manifests/fish_Pprolifica_AB_SraRunTable.csv \
  --datas1 "$DRYAD/DataS1.csv" \
  --geo-soft "$COHORT_DIR/reference/geo/GSE318355_family.soft.gz" \
  --out "$COHORT_DIR/manifests/fish_Pprolifica_AB_geo_crosswalk_review.tsv"

column -t -s $'\t' "$COHORT_DIR/manifests/fish_Pprolifica_AB_geo_crosswalk_review.tsv" | less -S
```

The script accepts only a literal, exact Dryad sample alias in the **corresponding** GEO sample's title, source, characteristics, or description. No exact alias means `NO_EXACT_ALIAS_IN_GEO` (manual source review); more than one is flagged as ambiguous. It never sets `include=YES` or invents `cluster_id`.

Even when all 11 aliases are matched, the final reviewed cohort manifests must distinguish the 6 reciprocal F1s from 5 parental controls; resolve replicate biological-independence and strand mode; and ensure marker reference compatibility. The normal `prepare_fish_manifests.py` expects ENA FASTQ metadata and is not yet an SRA-only manifest adapter. Do **not** submit `reference`, `pilot`, or `scale` until this limitation is resolved.
