#!/usr/bin/env python3
"""Offline manifest/build consistency check for the maize CyVerse inputs."""

import csv
from pathlib import Path


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


root = Path(__file__).resolve().parent.parent
resources = read(root / "manifests/maize_B73_Ki11_resources.tsv")
markers = read(root / "manifests/maize_B73_Ki11_marker_sources.tsv")
assert len(resources) == 13, len(resources)
assert len(markers) == 10, len(markers)
prefix = "https://data.cyverse.org/dav-anon/iplant/projects/commons_repo/curated/Qi_Sun_Zea_mays_haplotype_map_2018/hmp321_imputed"
for chrom in range(1, 11):
    row = next(x for x in resources if x["resource_id"] == f"hmp321_imputed_agpv3_chr{chrom}")
    assert row["url"] == f"{prefix}/merged_flt_c{chrom}.imputed.vcf.gz"
    assert row["checksum_url"] == f"{prefix}/MD5.txt"
    assert row["build"] == "AGPv3" and row["method"] == "cyverse_md5"
    assert row["md5"] == "NA" and row["bytes"] == "NA"
    source = next(x for x in markers if x["source_id"] == row["resource_id"])
    assert source["path"] == row["destination"] and source["build"] == "AGPv3"
    assert source["processed_path"].endswith(f"hmp321_imputed_agpv4_chr{chrom}.markers.tsv")
chain = next(x for x in resources if x["resource_id"] == "AGPv3_to_B73_RefGen_v4_chain")
assert chain["bytes"] == "5950740" and chain["md5"] == "51c5eafdccc3d17a8555e899a7dd5e2a"
assert chain["method"] == "https"
print("PASS: CyVerse 10 chromosome URLs, publisher MD5 workflow, AGPv3 markers and pinned chain")
