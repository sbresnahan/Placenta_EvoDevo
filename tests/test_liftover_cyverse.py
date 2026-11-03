#!/usr/bin/env python3
"""Offline synthetic AGPv3/AGPv4 chain tests (no R or external programs)."""

import csv
import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent.parent
    script = root / "bin" / "liftover_cyverse_markers.py"
    with tempfile.TemporaryDirectory(prefix="cyverse lift ") as temp:
        chain = Path(temp, "assembly.chain.gz")
        source = Path(temp, "chr1.vcf.gz")
        output = Path(temp, "lifted.tsv")
        with gzip.open(chain, "wt") as handle:
            handle.write("chain 100 1 100 + 10 15 1 100 + 20 25 1\n5\n\n")
            handle.write("chain 100 1 100 + 30 35 2 100 - 40 45 2\n5\n\n")
        with gzip.open(source, "wt") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB73\tKi11\t282set_B73\t282set_Ki11\n")
            fields = ["1", "11", ".", "A", "G", ".", "PASS", ".", "GT:GQ", "0/0:99", "1/1:99", "0/0:99", "1/1:99"]
            handle.write("\t".join(fields) + "\n")
            fields[1], fields[3], fields[4] = "31", "C", "T"
            handle.write("\t".join(fields) + "\n")
            fields[1], fields[3], fields[4], fields[-1] = "12", "A", "G", "0/0:99"  # discordant duplicate
            handle.write("\t".join(fields) + "\n")
            fields[1], fields[-1], fields[9] = "13", "1/1:99", "0/0:8"  # low confidence
            handle.write("\t".join(fields) + "\n")
            fields[1], fields[9] = "90", "0/0:99"  # unmapped
            handle.write("\t".join(fields) + "\n")
        subprocess.run([sys.executable, str(script), "--source", str(source),
                        "--chain", str(chain), "--output", str(output),
                        "--focal", "B73", "--alternate", "Ki11",
                        "--duplicate", "282set_B73", "--duplicate", "282set_Ki11",
                        "--audit", str(output) + ".audit.tsv"], check=True)
        with open(output, newline="") as handle:
            lines = list(csv.DictReader(handle, delimiter="\t"))
        assert len(lines) == 2, lines
        assert (lines[0]["chrom"], lines[0]["genomic_pos"], lines[0]["ref"], lines[0]["alt"]) == ("1", "21", "A", "G")
        assert (lines[1]["chrom"], lines[1]["genomic_pos"], lines[1]["ref"], lines[1]["alt"]) == ("2", "60", "G", "A")
        assert all(x["focal_index"] == "0" and x["alternate_index"] == "1" for x in lines)
        assert os.stat(str(output) + ".audit.tsv").st_size > 0
    print("PASS: offline CyVerse AGPv3->AGPv4 chain, reverse complement, parental filters")


if __name__ == "__main__":
    main()
