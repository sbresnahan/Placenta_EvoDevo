#!/usr/bin/env python3
"""Host-only AGPv3 -> B73 RefGen_v4 parental-SNP extraction using Ensembl chain.

Streams a CyVerse imputed VCF, maps single-base coordinates, and atomically
publishes a small marker TSV. The containerized R stage validates the target
REF bases before accepting these coordinates.
"""

import argparse
import bisect
import csv
import gzip
import os
import re
import tempfile
from collections import Counter, defaultdict

BASES = set("ACGT")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
CHROMS = {str(i) for i in range(1, 11)}


def chrom_name(value):
    return re.sub(r"^chr", "", value, flags=re.IGNORECASE)


def load_chain(path):
    """Index AGPv3 source intervals; retain all alignments to reject ambiguity."""
    blocks = defaultdict(list)
    with gzip.open(path, "rt") as handle:
        fields = None
        t, q = 0, 0
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("chain "):
                x = line.split()
                if len(x) != 13 or x[4] != "+" or x[9] not in ("+", "-"):
                    raise ValueError("Unsupported assembly chain header")
                fields = x
                t, q = int(x[5]), int(x[10])
                continue
            if fields is None:
                raise ValueError("Malformed chain: missing header")
            step = line.split()
            if len(step) not in (1, 3):
                raise ValueError("Malformed chain block")
            size = int(step[0])
            source = chrom_name(fields[2])
            target = chrom_name(fields[7])
            if source in CHROMS and target in CHROMS:
                blocks[source].append((t, t + size, target, q, int(fields[8]), fields[9], int(fields[1])))
            t += size
            q += size
            if len(step) == 3:
                t += int(step[1])
                q += int(step[2])
            else:
                fields = None
    indexes = {}
    for chrom, entries in blocks.items():
        entries.sort(key=lambda row: row[0])
        starts, max_ends = [], []
        maximum = 0
        for entry in entries:
            starts.append(entry[0])
            maximum = max(maximum, entry[1])
            max_ends.append(maximum)
        indexes[chrom] = (entries, starts, max_ends)
    if not indexes:
        raise ValueError("No chromosome chains in AGPv3-to-AGPv4 assembly chain")
    return indexes


def map_snp(index, chrom, pos):
    chrom = chrom_name(chrom)
    if chrom not in index:
        return None
    entries, starts, max_ends = index[chrom]
    pos0 = pos - 1
    j = bisect.bisect_right(starts, pos0) - 1
    hits = set()
    while j >= 0 and max_ends[j] > pos0:
        start, end, target, q_start, q_size, strand, _score = entries[j]
        if start <= pos0 < end:
            offset = q_start + pos0 - start
            dest0 = offset if strand == "+" else q_size - 1 - offset
            hits.add((target, dest0 + 1, strand))
            if len(hits) > 1:
                return "ambiguous"
        j -= 1
    return next(iter(hits)) if hits else None


def parent_call(value, keys, min_gq):
    parts = value.split(":")
    fields = dict(zip(keys, parts))
    gt = fields.get("GT", ".")
    if gt in ("0/0", "0|0"):
        allele = 0
    elif gt in ("1/1", "1|1"):
        allele = 1
    else:
        return None
    # Imputed hapmap data are not raw reads: AD is not a meaningful filter.
    # Enforce GQ when emitted, while retaining sites without that optional tag.
    gq = fields.get("GQ", ".")
    if gq not in ("", "."):
        try:
            if float(gq) < min_gq:
                return None
        except ValueError:
            return None
    return allele


def process(args):
    index = load_chain(args.chain)
    counts = Counter()
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    handle, tmp = tempfile.mkstemp(prefix=".agpv4_marker_", suffix=".tsv", dir=os.path.dirname(os.path.abspath(args.output)))
    try:
        with os.fdopen(handle, "w", newline="") as out, gzip.open(args.source, "rt") as vcf:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(("chrom", "genomic_pos", "ref", "alt", "focal_index", "alternate_index"))
            expected = [args.focal, args.alternate] + args.duplicates
            indexes = None
            for line in vcf:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    header = line.rstrip("\r\n").split("\t")
                    missing = [name for name in expected if name not in header]
                    if missing:
                        raise ValueError("Missing parental/duplicate VCF samples: " + ", ".join(missing))
                    indexes = [header.index(name) for name in expected]
                    continue
                if indexes is None:
                    raise ValueError("VCF header is missing")
                counts["scanned"] += 1
                x = line.rstrip("\r\n").split("\t")
                if len(x) < len(header):
                    raise ValueError("Malformed VCF row with missing sample columns")
                if x[3] not in BASES or x[4] not in BASES or x[3] == x[4]:
                    counts["not_biallelic_snp"] += 1
                    continue
                if x[6] not in ("PASS", ".") or re.search(r"(?:^|;)NI5(?:;|$)", x[7]):
                    counts["filtered"] += 1
                    continue
                keys = x[8].split(":")
                a, b, *dups = (parent_call(x[j], keys, args.min_gq) for j in indexes)
                if a is None or b is None or a == b:
                    counts["not_opposite_homozygous"] += 1
                    continue
                if any(v is not None and v != (a if i == 0 else b) for i, v in enumerate(dups)):
                    counts["duplicate_discordance"] += 1
                    continue
                try:
                    pos = int(x[1])
                except ValueError:
                    raise ValueError("Non-numeric VCF coordinate") from None
                mapped = map_snp(index, x[0], pos)
                if mapped is None:
                    counts["unmapped"] += 1
                    continue
                if mapped == "ambiguous":
                    counts["ambiguous"] += 1
                    continue
                chrom, new_pos, strand = mapped
                ref, alt = (x[3], x[4]) if strand == "+" else (x[3].translate(COMPLEMENT), x[4].translate(COMPLEMENT))
                writer.writerow((chrom, new_pos, ref, alt, a, b))
                counts["lifted"] += 1
        if counts["lifted"] == 0:
            raise ValueError("No callable parental SNPs mapped to AGPv4; inspect VCF samples/FORMAT and chain")
        os.replace(tmp, args.output)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    if args.audit:
        with open(args.audit + ".tmp", "w") as audit:
            audit.write("metric\tcount\n")
            for name in sorted(counts):
                audit.write("%s\t%d\n" % (name, counts[name]))
        os.replace(args.audit + ".tmp", args.audit)
    print("AGPv3->AGPv4 %s: %s markers, %s unmapped, %s ambiguous" % (
        os.path.basename(args.source), counts["lifted"], counts["unmapped"], counts["ambiguous"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "chain", "output", "focal", "alternate"):
        parser.add_argument("--" + name.replace("_", "-"), required=True)
    parser.add_argument("--duplicate", dest="duplicates", action="append", default=[])
    parser.add_argument("--audit")
    parser.add_argument("--min-gq", type=float, default=20)
    args = parser.parse_args()
    if len(args.duplicates) not in (0, 2) or args.min_gq < 0:
        parser.error("Exactly zero or two duplicate samples, and nonnegative min-GQ, required")
    process(args)


if __name__ == "__main__":
    main()
