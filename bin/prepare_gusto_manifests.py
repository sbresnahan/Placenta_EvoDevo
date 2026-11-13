#!/usr/bin/env python3
"""Create PRIVATE GUSTO input manifests from the verified RNA/VCF linkage.

This does not read sequence data, infer parentage or invent FASTQ read counts.
Run on Seadragon with access to all private source paths. Outputs must remain
outside git; existing manifests are never overwritten without --overwrite.
"""
import argparse
import csv
import gzip
from pathlib import Path


def table(path, header=False):
    with open(path, newline="") as stream:
        return list(csv.reader(stream, delimiter="\t"))


def vcf_samples(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        for line in stream:
            if line.startswith("#CHROM\t"):
                return line.rstrip("\n\r").split("\t")[9:]
    raise ValueError("VCF has no #CHROM sample header")


def build(mapping, matrix, samples, fastq_dir, strand, stage, check_fastq=True):
    by_child = {}
    by_run = set()
    for row in mapping:
        if len(row) != 2:
            raise ValueError("FASTQ ID mapping must have two columns without a header")
        infant, fastq = row
        child = "B" + infant
        run = "J" + fastq
        if child in by_child or run in by_run or not fastq.isdigit():
            raise ValueError("Duplicate or invalid infant/FASTQ identity")
        by_child[child] = run
        by_run.add(run)
    if not matrix or matrix[0][0] != "chr":
        raise ValueError("Invalid chromosome origin matrix header")
    matrix_ids = {col[:-2] for col in matrix[0][1:] if col.endswith("_L")}
    if any((child + "_R") not in matrix[0] for child in matrix_ids):
        raise ValueError("Incomplete Left/Right origin columns")
    seen, selected = set(), []
    for vcf_id in samples:
        if "_" not in vcf_id:
            raise ValueError("Unrecognized GUSTO VCF sample name: " + vcf_id)
        _, child = vcf_id.rsplit("_", 1)
        if child not in by_child or child not in matrix_ids:
            raise ValueError("VCF/matrix/FASTQ sample linkage conflict: " + vcf_id)
        if child in seen:
            raise ValueError("Duplicate VCF child: " + child)
        seen.add(child)
        selected.append((by_child[child], child))
    selected.sort()
    sources, libraries, units, trios, covariates = [], [], [], [], []
    for run, child in selected:
        mates = [fastq_dir / f"{run}_{mate}.fq.gz" for mate in (1, 2)]
        if check_fastq and any(not p.is_file() or p.stat().st_size == 0 for p in mates):
            raise FileNotFoundError(f"Missing or empty FASTQ for {run}: {mates}")
        sizes = [p.stat().st_size if p.is_file() else 0 for p in mates]
        sources.append(dict(run_id=run, sample_id=run, layout="PE", url_1="", url_2="",
                            bytes_1=sizes[0], bytes_2=sizes[1], read_count="",
                            local_1=str(mates[0]), local_2=str(mates[1])))
        libraries.append(dict(sample_id=run, unit_id=run, child_id=child,
                              layout="PE", strand_mode=strand))
        units.append(dict(unit_id=run, cluster_id=child,
                          maternal_parent="precomputed_maternal", paternal_parent="precomputed_paternal",
                          stage=stage, tissue="placenta"))
        # One infant/family is provisionally assumed; audit relatedness before models.
        trios.append(dict(child_id=child, family_id=child))
        covariates.append(dict(sample_id=run))
    absent = sorted(set(by_child) - seen)
    return (dict(sources=sources, libraries=libraries, units=units,
                 trios=trios, covariates=covariates), absent)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--origin-matrix", required=True, type=Path)
    parser.add_argument("--id-map", required=True, type=Path)
    parser.add_argument("--fastq-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--strand-mode", required=True, choices=("fr", "rf", "unstranded"))
    parser.add_argument("--stage", default="unspecified")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-fastq-check", action="store_true", help="Only for preparing offline fixtures")
    args = parser.parse_args()
    records, absent = build(table(args.id_map), table(args.origin_matrix),
                            vcf_samples(args.vcf), args.fastq_dir, args.strand_mode,
                            args.stage, not args.skip_fastq_check)
    destinations = {key: args.out_dir / f"GUSTO_{key}.tsv" for key in records}
    if not args.overwrite and any(p.exists() for p in destinations.values()):
        raise FileExistsError("Existing GUSTO manifests; inspect then explicitly use --overwrite")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    columns = {
        "sources": ("run_id", "sample_id", "layout", "url_1", "url_2", "bytes_1",
                    "bytes_2", "read_count", "local_1", "local_2"),
        "libraries": ("sample_id", "unit_id", "child_id", "layout", "strand_mode"),
        "units": ("unit_id", "cluster_id", "maternal_parent", "paternal_parent", "stage", "tissue"),
        "trios": ("child_id", "family_id"),
        "covariates": ("sample_id",),
    }
    for key, rows in records.items():
        with destinations[key].open("w", newline="") as out:
            writer = csv.DictWriter(out, delimiter="\t", fieldnames=columns[key], lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    print(f"Wrote {len(records['libraries'])} private GUSTO libraries to {args.out_dir}")
    print(f"Absent from VCF ({len(absent)}): {', '.join(absent)}")
    print("FASTQ read_count intentionally blank; no read totals were invented.")
    print("Confirm strand mode, independent family IDs, cohort stage and FASTQ paths before submission.")


if __name__ == "__main__":
    main()
