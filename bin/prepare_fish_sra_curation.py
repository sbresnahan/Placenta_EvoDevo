#!/usr/bin/env python3
"""Scaffold fish cross assignments from public NCBI SRA RunInfo (no guessing)."""
import argparse
import csv
from pathlib import Path

FIELDS = ("run_accession", "geo_accession", "biosample", "experiment_accession",
          "source_tissue", "include", "unit_id", "maternal_parent", "paternal_parent",
          "stage", "tissue", "cluster_id", "strand_mode", "notes")


def convert(input_path):
    with input_path.open(newline="", encoding="utf-8-sig") as handle:
        runinfo = list(csv.DictReader(handle))
    assert len(runinfo) == 29, f"Expected 29 RNA runs, found {len(runinfo)}"
    assert len({r['Run'] for r in runinfo}) == len(runinfo), "Duplicate SRR"
    out = []
    for r in sorted(runinfo, key=lambda d: d["Run"]):
        assert r["SRA Study"] == "SRP674561"
        assert r["LibraryLayout"] == "PAIRED"
        assert r["Organism"] == "Poeciliopsis prolifica"
        tissue = r["tissue"].strip().lower()
        assert tissue in {"placenta", "embryo", "follicle"}, (r["Run"], tissue)
        out.append(dict(run_accession=r["Run"], geo_accession=r["Sample Name"],
                        biosample=r["BioSample"], experiment_accession=r["Experiment"],
                        source_tissue=tissue, include="NO", unit_id="",
                        maternal_parent="", paternal_parent="", stage="",
                        tissue="pericardial_sac" if tissue == "placenta" else tissue,
                        cluster_id="", strand_mode="",
                        notes=("Cross direction and biological unit unverified"
                               if tissue == "placenta" else
                               "Nonplacental; exclude from primary POE")))
    assert sum(x["source_tissue"] == "placenta" for x in out) == 11
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runinfo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Refusing to overwrite existing curation: {args.output}")
    rows = convert(args.runinfo)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} SRA runs, including 11 candidate placenta runs: {args.output}")
    print("All include=NO: manually verify crosses and biological units using GEO before analysis")


if __name__ == '__main__':
    main()
