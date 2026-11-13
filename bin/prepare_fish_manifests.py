#!/usr/bin/env python3
"""Build reviewed P. prolifica reciprocal-cross manifests from ENA RunInfo.

Only rows explicitly marked YES in a manually reviewed curation table are used.
Parent-of-origin and placental tissue are never guessed from SRA/GEO labels.
This script does not download reads, infer parental genotypes, or process HiFi data.
"""

import argparse
import csv
import pathlib
import re
import sys

NAME = "fish_Pprolifica_AB"
PARENTS = {"Line_A", "Line_B"}
ID = re.compile(r"^[A-Za-z0-9_.-]+$")

SOURCES = ("run_id", "sample_id", "study_accession", "sample_accession",
           "experiment_accession", "library_name", "layout", "read_count",
           "base_count", "url_1", "md5_1", "bytes_1", "url_2", "md5_2", "bytes_2")
LIBRARIES = ("sample_id", "unit_id", "layout", "strand_mode")
UNITS = ("unit_id", "cluster_id", "maternal_parent", "paternal_parent", "stage", "tissue")
CURATION = ("run_accession", "include", "unit_id", "maternal_parent",
            "paternal_parent", "stage", "tissue", "cluster_id", "strand_mode",
            "sample_title", "sample_alias")
REQUIRED_ENA = ("run_accession", "study_accession", "sample_accession",
                "experiment_accession", "library_layout", "read_count",
                "base_count", "fastq_ftp", "fastq_bytes")


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing header: {path}")
        rows = list(reader)
    if any(not row for row in rows):
        raise ValueError(f"Empty row in {path}")
    return rows, reader.fieldnames


def safe_id(value, label):
    if not value or not ID.fullmatch(value):
        raise ValueError(f"Missing or unsafe {label}: {value!r}")
    return value


def positive_int(value, label):
    try:
        val = int(value)
    except (ValueError, TypeError):
        raise ValueError(f"Invalid {label}: {value!r}") from None
    if val < 1:
        raise ValueError(f"Invalid {label}: {value!r}")
    return val


def fetch_urls(row, layout):
    urls = (row.get("fastq_ftp") or "").split(";")
    sizes = (row.get("fastq_bytes") or "").split(";")
    hashes = (row.get("fastq_md5") or "").split(";")
    n = 2 if layout == "PE" else 1
    if len(urls) != n or len(sizes) != n:
        raise ValueError(f"{row['run_accession']}: ENA FASTQ mate metadata incomplete")
    result = []
    for i in range(n):
        uri = urls[i]
        if not uri.startswith(("ftp.sra.ebi.ac.uk/", "era-fasp@", "https://")):
            raise ValueError(f"{row['run_accession']}: unexpected ENA FASTQ location {uri!r}")
        if uri.startswith("era-fasp@"):
            raise ValueError("Aspera ENA URL cannot be used as HTTPS")
        if uri.startswith("ftp.sra.ebi.ac.uk/"):
            uri = "https://" + uri
        if not uri.startswith("https://ftp.sra.ebi.ac.uk/"):
            raise ValueError("Expected official ENA FASTQ host")
        positive_int(sizes[i], "FASTQ byte count")
        result.append((uri, hashes[i] if i < len(hashes) else "", sizes[i]))
    return result


def save_tsv(path, fields, records, force=False):
    if path.exists() and not force:
        old, _ = read_tsv(path)
        if old:
            raise ValueError(f"Refusing to overwrite populated manifest: {path} (use --force)")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def scaffold(ena, target, force):
    rows = []
    for record in ena:
        rows.append(dict(run_accession=record["run_accession"], include="NO",
                         unit_id="", maternal_parent="", paternal_parent="", stage="",
                         tissue="", cluster_id="", strand_mode="",
                         sample_title=record.get("sample_title", ""),
                         sample_alias=record.get("sample_alias", "")))
    save_tsv(target, CURATION, rows, force=force)
    print(f"Review {target}: mark only placental RNA runs YES and fill ALL cross metadata.")


def build(ena, curated, dest, force):
    curated_rows, headers = read_tsv(curated)
    if not set(CURATION).issubset(headers):
        raise ValueError("Curation table missing required columns")
    by_run = {r["run_accession"]: r for r in ena}
    if len(by_run) != len(ena):
        raise ValueError("Duplicate run accession in ENA data")
    seen = set()
    sources, libraries, units = [], {}, {}
    for selected in curated_rows:
        run = safe_id(selected["run_accession"], "run accession")
        if run in seen or run not in by_run:
            raise ValueError(f"Duplicate/unknown curated run: {run}")
        seen.add(run)
        if selected["include"] not in {"YES", "NO"}:
            raise ValueError(f"{run}: include must be YES or NO")
        if selected["include"] == "NO":
            continue
        unit = safe_id(selected["unit_id"], "unit ID")
        maternal, paternal = selected["maternal_parent"], selected["paternal_parent"]
        if {maternal, paternal} != PARENTS:
            raise ValueError(f"{run}: explicit Line_A/Line_B reciprocal parents required")
        if selected["tissue"] != "pericardial_sac":
            raise ValueError(f"{run}: include only fetal placental pericardial sac, not embryo/follicle")
        stage = safe_id(selected["stage"], "developmental stage")
        strand = selected["strand_mode"]
        if strand not in {"fr", "rf", "unstranded"}:
            raise ValueError(f"{run}: strand_mode must be checked against the actual library")
        cluster = selected["cluster_id"] or "NA"
        item = by_run[run]
        layout = item["library_layout"]
        if layout not in {"SINGLE", "PAIRED"}:
            raise ValueError(f"{run}: unsupported ENA layout {layout!r}")
        layout = "PE" if layout == "PAIRED" else "SE"
        fq = fetch_urls(item, layout)
        unit_record = dict(unit_id=unit, cluster_id=cluster, maternal_parent=maternal,
                           paternal_parent=paternal, stage=stage, tissue=selected["tissue"])
        lib_record = dict(sample_id=unit, unit_id=unit, layout=layout, strand_mode=strand)
        if unit in units and (units[unit] != unit_record or libraries[unit] != lib_record):
            raise ValueError(f"{run}: conflicting metadata within biological unit {unit}")
        units[unit] = unit_record
        libraries[unit] = lib_record
        sources.append(dict(run_id=run, sample_id=unit,
                            study_accession=item["study_accession"],
                            sample_accession=item["sample_accession"],
                            experiment_accession=item["experiment_accession"],
                            library_name=item.get("library_name", ""),
                            layout=layout, read_count=str(positive_int(item["read_count"], "read_count")),
                            base_count=str(positive_int(item["base_count"], "base_count")),
                            url_1=fq[0][0], md5_1=fq[0][1], bytes_1=fq[0][2],
                            url_2=fq[1][0] if layout == "PE" else "NA",
                            md5_2=fq[1][1] if layout == "PE" else "NA",
                            bytes_2=fq[1][2] if layout == "PE" else "NA"))
    if len(units) < 3:
        raise ValueError("At least 3 curated biological units required for this cohort")
    directions = {r["maternal_parent"] for r in units.values()}
    if directions != PARENTS:
        raise ValueError("Both reciprocal directions required")
    if len({r["stage"] for r in units.values()}) != 1:
        raise ValueError("Separate mixed developmental stages into different cohort configurations")
    outputs = ((SOURCES, sources, "sources"), (LIBRARIES, list(libraries.values()), "libraries"),
               (UNITS, list(units.values()), "units"))
    # Verify all targets before writing any output.
    for _, _, suffix in outputs:
        target = dest / f"{NAME}_{suffix}.tsv"
        if target.exists() and not force and read_tsv(target)[0]:
            raise ValueError(f"Refusing to overwrite populated manifest: {target}")
    for columns, rows, suffix in outputs:
        save_tsv(dest / f"{NAME}_{suffix}.tsv", columns, rows, force=force)
    print(f"Prepared {len(sources)} RNA runs across {len(units)} reviewed placental units.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runinfo", type=pathlib.Path, required=True, help="ENA read_run TSV")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--scaffold", type=pathlib.Path, help="Write curation table for manual review")
    mode.add_argument("--build", type=pathlib.Path, help="Use reviewed curation table")
    parser.add_argument("--output-dir", type=pathlib.Path, help="Working cohort manifest directory")
    parser.add_argument("--force", action="store_true", help="Overwrite populated target manifests")
    args = parser.parse_args()
    ena, cols = read_tsv(args.runinfo)
    if not set(REQUIRED_ENA).issubset(cols):
        raise ValueError("RunInfo is missing required ENA fields")
    if not ena:
        raise ValueError("No runs found in ENA metadata")
    if args.scaffold:
        scaffold(ena, args.scaffold, args.force)
    else:
        if not args.output_dir:
            parser.error("--build requires --output-dir")
        build(ena, args.build, args.output_dir, args.force)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        print(f"Fish manifest preparation failed: {exc}", file=sys.stderr)
        sys.exit(1)
