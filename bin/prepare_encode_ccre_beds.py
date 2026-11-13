#!/usr/bin/env python3
"""Split an ENCODE per-biosample cCRE bigBed into BED inputs for human-context.

ENCODE v2 biosample cCRE files (bed9 + "ccRE Group" + "Classification") list
every registry element and label its activity in that biosample. Elements
labelled Low-DNase are inactive there and are dropped; each active group
(PLS, pELS, dELS, DNase-H3K4me3, DNase-only, CTCF-only, ...) becomes its own
gzipped BED with its own evidence_type, so peak overlaps stay separable.

Writes to a new --out-dir:
  <ACCESSION>_<group>.bed.gz   0-based BED, chrom/start/end/cCRE ID
  beds_manifest.tsv            regulatory.human_context.beds manifest
  summary.tsv                  elements per group (written and dropped)

Reads the bigBed directly (no UCSC tools needed). Python 3.6+, standard
library only. Never overwrites an existing --out-dir.
"""
import argparse
import csv
import gzip
import struct
import sys
import zlib
from collections import Counter, OrderedDict
from pathlib import Path

INACTIVE = {"Low-DNase"}


def read_bigbed(path):
    """Yield (chrom, start, end, rest_fields) for every record in a bigBed."""
    data = Path(path).read_bytes()
    if struct.unpack_from("<I", data, 0)[0] != 0x8789F2EB:
        raise ValueError("Not a little-endian bigBed file: {}".format(path))
    (_version, _zooms, chrom_tree, _full_data, full_index, _fields, _defined,
     _autosql, _summary, uncompress_buf, _ext) = struct.unpack_from("<HHQQQHHQQIQ", data, 4)

    magic, _block, key_size, _val, _n, _r = struct.unpack_from("<IIIIQQ", data, chrom_tree)
    if magic != 0x78CA8C91:
        raise ValueError("Corrupt bigBed chromosome tree")
    chroms = {}

    def chrom_walk(offset):
        leaf, _, count = struct.unpack_from("<BBH", data, offset)
        offset += 4
        for _ in range(count):
            key = data[offset:offset + key_size].rstrip(b"\0").decode()
            offset += key_size
            if leaf:
                chrom_id, _size = struct.unpack_from("<II", data, offset)
                chroms[chrom_id] = key
                offset += 8
            else:
                chrom_walk(struct.unpack_from("<Q", data, offset)[0])
                offset += 8
    chrom_walk(chrom_tree + 32)

    if struct.unpack_from("<I", data, full_index)[0] != 0x2468ACE0:
        raise ValueError("Corrupt bigBed R-tree index")
    blocks = []

    def index_walk(offset):
        leaf, _, count = struct.unpack_from("<BBH", data, offset)
        offset += 4
        for _ in range(count):
            if leaf:
                blocks.append(struct.unpack_from("<QQ", data, offset + 16))
                offset += 32
            else:
                index_walk(struct.unpack_from("<Q", data, offset + 16)[0])
                offset += 24
    index_walk(full_index + 48)

    for block_offset, block_size in sorted(blocks):
        raw = data[block_offset:block_offset + block_size]
        buf = zlib.decompress(raw) if uncompress_buf else raw
        pos = 0
        while pos < len(buf):
            chrom_id, start, end = struct.unpack_from("<III", buf, pos)
            pos += 12
            stop = buf.index(b"\0", pos)
            rest = buf[pos:stop].decode().split("\t") if stop > pos else []
            pos = stop + 1
            yield chroms[chrom_id], start, end, rest


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--bigbed", required=True, type=Path, help="ENCODE cCRE bigBed (e.g. ENCFF603ENH.bigBed)")
    p.add_argument("--accession", required=True, help="ENCODE annotation accession, e.g. ENCSR057ORE")
    p.add_argument("--out-dir", required=True, type=Path, help="new output directory")
    p.add_argument("--build", default="GRCh38")
    a = p.parse_args(argv)
    if a.out_dir.exists():
        p.error("Refusing to overwrite existing --out-dir: {}".format(a.out_dir))

    out = a.out_dir.resolve()
    out.mkdir(parents=True)
    handles, counts, dropped = OrderedDict(), Counter(), Counter()
    try:
        for chrom, start, end, rest in read_bigbed(a.bigbed):
            # bed9 fields 4-9 come first; the cCRE group is the 10th BED column.
            if len(rest) < 7:
                raise ValueError("Expected bed9 + cCRE group columns; got {} extra fields".format(len(rest)))
            name, group = rest[0], rest[6]
            if group in INACTIVE:
                dropped[group] += 1
                continue
            safe = group.replace("/", "_")
            if safe not in handles:
                handles[safe] = gzip.open(str(out / "{}_{}.bed.gz".format(a.accession, safe)), "wt")
            handles[safe].write("{}\t{}\t{}\t{}\n".format(chrom, start, end, name))
            counts[safe] += 1
    finally:
        for h in handles.values():
            h.close()
    if not counts:
        raise ValueError("No active cCREs found in {}".format(a.bigbed))

    with open(str(out / "beds_manifest.tsv"), "w", newline="") as stream:
        w = csv.writer(stream, delimiter="\t", lineterminator="\n")
        w.writerow(["input_id", "path", "build", "evidence_type", "selection_used"])
        for group in counts:
            w.writerow(["ENCODE_{}_{}".format(a.accession, group),
                        str(out / "{}_{}.bed.gz".format(a.accession, group)),
                        a.build, "ENCODE_cCRE_{}".format(group), "FALSE"])
    with open(str(out / "summary.tsv"), "w", newline="") as stream:
        w = csv.writer(stream, delimiter="\t", lineterminator="\n")
        w.writerow(["ccre_group", "n_elements", "written"])
        for group, n in counts.most_common():
            w.writerow([group, n, "yes"])
        for group, n in dropped.most_common():
            w.writerow([group, n, "no_inactive_in_biosample"])
    for group, n in counts.most_common():
        print("{:<16} {:>8}".format(group, n))
    for group, n in dropped.most_common():
        print("{:<16} {:>8}  (dropped: inactive)".format(group, n))
    print("\nbeds: {}".format(out / "beds_manifest.tsv"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
