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
    # Legacy manifest values remain supported, but no catalog is consulted.
    assert row["build"] == "AGPv3" and row["method"] in {"cyverse", "cyverse_md5"}
    assert row["bytes"] == "NA"
    source = next(x for x in markers if x["source_id"] == row["resource_id"])
    assert source["path"] == row["destination"] and source["build"] == "AGPv3"
    assert source["processed_path"].endswith(f"hmp321_imputed_agpv4_chr{chrom}.markers.tsv")
chain = next(x for x in resources if x["resource_id"] == "AGPv3_to_B73_RefGen_v4_chain")
assert chain["bytes"] == "5950740"
assert chain["method"] == "https"

# B73/W22 reuses the exact Ki11 VCF inputs and chain, but produces its own
# parent-specific marker outputs. No CyVerse genotype fetch is configured.
w22 = read(root / "manifests/maize_B73_W22_marker_sources.tsv")
w22_resources = read(root / "manifests/maize_B73_W22_resources.tsv")
assert len(w22) == 10, len(w22)
assert not any("hmp321" in row["resource_id"] for row in w22_resources)
for chrom in range(1, 11):
    name = f"hmp321_imputed_agpv3_chr{chrom}"
    ki11 = next(row for row in markers if row["source_id"] == name)
    source = next(row for row in w22 if row["source_id"] == name)
    assert source["path"] == ki11["path"]
    assert source["build"] == "AGPv3"
    assert source["processed_path"] != ki11["processed_path"]
    assert source["processed_path"].endswith(
        f"maize_B73_W22/reference/hmp321_imputed_agpv4_chr{chrom}.markers.tsv"
    )

ki11_config = (root / "config/maize_B73_Ki11.example.yaml").read_text()
w22_config = (root / "config/maize_B73_W22.example.yaml").read_text()
ki11_chain = next(line.strip() for line in ki11_config.splitlines()
                  if line.strip().startswith("liftover_chain:"))
assert ki11_chain in w22_config
assert "markers:\n  mode: panzea" in w22_config
assert "alternate_sample: W22" in w22_config
assert "duplicate_samples: []" in w22_config
print("PASS: CyVerse 10 chromosome URLs, shared B73/W22 HapMap3 markers, and chain")
