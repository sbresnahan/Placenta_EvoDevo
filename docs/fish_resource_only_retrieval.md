# Fish RefSeq resource-only retrieval

`fish_Pprolifica_AB` initially has header-only sample linkage manifests. The normal
`retrieve` command deliberately requires verified sample/library/unit links, and
must continue to do so before the fish reference and quantification stages. The
already-staged 11 SRA runs are not a substitute for that crosswalk.

Use the explicit `retrieve --resources-only` mode ONLY to stage the exact RefSeq
FASTA, annotation and assembly report through the existing LSF retrieval worker,
content-length verification and `resource_receipt.tsv`. This mode skips only
sample linkage validation; it does not relax downstream model gates.

## Network-enabled manifest preparation

Set `PROJECT`, `COHORT_DIR`, `SCRIPTS_DIR`, `CONFIG`, and `POE_SITE_FILE` as usual.
The manifest helper rejects existing rows and validates the final response's
Content-Length for each URL before writing the cohort-local manifest. Execute on
a transfer host because login hosts may be network-restricted:

```bash
mkdir -p "$COHORT_DIR/logs/lsf"
bsub -q transfer -n 1 -M 2 -R 'rusage[mem=2]' -W 00:30 \
  -J fish_refseq_manifest \
  -o "$COHORT_DIR/logs/lsf/fish_refseq_manifest.%J.out" \
  -e "$COHORT_DIR/logs/lsf/fish_refseq_manifest.%J.err" \
  -env "all,COHORT_DIR=$COHORT_DIR" \
  bash "$SCRIPTS_DIR/bin/prepare_fish_refseq_resources.sh"
```

Only after this job succeeds, inspect
`$COHORT_DIR/manifests/fish_Pprolifica_AB_resources.tsv`. If NCBI does not serve
one of the three named resources, the helper stops without writing a partial
manifest. The code does not invent file sizes or silently substitute assemblies.

## Submit within the existing retrieval pipeline

```bash
bin/submit_pipeline.sh "$CONFIG" retrieve --resources-only
bin/submit_pipeline.sh "$CONFIG" retrieve --resources-only --submit
```

The default `retrieve` invocation still requires valid sample manifests; it
should not be used until sample crosswalk curation is complete.

After the transfer finishes, confirm the genome and annotation are nonempty and
check the receipt. Index `genome.fa` for the contig audit:

```bash
ls -lh "$COHORT_DIR/reference/genome.fa" "$COHORT_DIR/reference/source_annotation.gtf"
ls -lh "$COHORT_DIR/reference/"*assembly_report.txt
ls -lh "$COHORT_DIR/resource_receipt.tsv"
bin/native/samtools faidx "$COHORT_DIR/reference/genome.fa"
```

CAUTION: The Dryad `SCAFFOLD...` coordinates must still be reconciled to this
exact FASTA, with allele-by-allele validation; and the 11 GSM/SRR records must
still be matched to the original Line A/B and developmental-stage sample labels.
The normal reference / pilot / scale / models validation remains in force.
