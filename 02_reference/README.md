# 02 — Reference preparation

Inputs: compatible genome, annotation, build/contig settings, and parental
marker sources or typed-child inputs.

```bash
bin/submit_pipeline.sh "$CONFIG" reference --submit
```

The graph builds the transcript/promoter catalog, prepares genotypes, and
builds Ornaments indices. Public cohorts use three dependent jobs.
Human cohorts use one catalog plus preparation/index jobs per child.

Outputs include normalized annotation, transcripts, coordinate maps,
promoter events/pairs, annotation audit, markers and indices.
The human shared genome is copied into the cohort before indexing.

Next: [pilot quantification](../04_ornaments/README.md).
