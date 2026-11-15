#!/usr/bin/env python3
"""Prepare PANTRY cis-QTL results for the GUSTO human-context module.

For each ancestry and modality this keeps every significant permutation lead
(qval < --qval in *_top.tsv) plus every conditionally independent signal
(*_independent_top.tsv), from both the grouped and the ungrouped PANTRY runs
where they exist, de-duplicated by (phenotype_id, variant_id). It writes:

  qtl/<ANCESTRY>_<modality>_lead.tsv  pipeline-ready tables (modality, chrom, pos added)
  qtl_inputs.tsv                      regulatory.human_context.qtl_tables manifest
  phenotype_map.tsv                   regulatory.human_context.phenotype_map crosswalk
  summary.tsv                         rows, phenotypes and mapping outcome per table

chrom and pos are parsed from PANTRY's CHROM:POS:REF:ALT variant IDs; effects
are not flipped and REF/ALT roles are not reinterpreted.

Phenotype map, in order of preference (GUSTO IDs come from its event_map):
  event  verified   TSS-linked modalities only (alt_TSS, isoforms,
                    isoform_expression): the PANTRY transcript is a GUSTO
                    transcript of the same gene; mapped to its TSS event
  gene   verified   otherwise, the PANTRY gene ID matches a GUSTO gene exactly,
                    or by unversioned Ensembl ID when that match is unique
  -      ambiguous  unversioned gene ID matches several GUSTO genes
  -      unmapped   no GUSTO gene
Gene-only modalities (expression, stability, alt_polyA, splicing,
intron_retention, RNA_editing) are never linked to a TSS event: their
phenotypes describe a gene, 3' end, intron or site, not a TSS. Gene-level rows
are gene context only, not validation of a TSS phenotype. PANTRY "combined"
pools the other modalities and is not supported, to avoid double counting.

Python 3.6+, standard library only. Never overwrites an existing --out-dir.
"""
import argparse
import csv
import re
import sys
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path

VARIANT = re.compile(r"^([^:]+):([0-9]+):([^:]+):([^:]+)$")
# Per modality: PANTRY run variants in priority order (first kept on duplicates;
# at least one must exist), and whether phenotypes may link to a TSS event.
BOTH = ("ungrouped_cisqtl", "cisqtl")
MODALITIES = OrderedDict([
    ("alt_TSS", (BOTH, "tss")),
    ("isoforms", (BOTH, "tss")),
    ("isoform_expression", (BOTH, "tss")),
    ("expression", (BOTH, "gene")),
    ("stability", (BOTH, "gene")),
    ("alt_polyA", (BOTH, "gene")),
    ("splicing", (BOTH, "gene")),
    ("intron_retention", (BOTH, "gene")),
    ("RNA_editing", (BOTH, "gene")),
])
KEEP = ("pval_nominal", "slope", "slope_se", "af", "ma_samples", "start_distance",
        "pval_beta", "qval", "rank")


def read_tsv(path):
    with open(str(path), newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows:
        raise ValueError("Empty table: {}".format(path))
    return rows


def transcript_of(modality, phenotype):
    """PANTRY transcript for TSS-linked modalities; None for gene-only ones."""
    if MODALITIES[modality][1] != "tss":
        return None
    if "__" not in phenotype:
        raise ValueError("Unexpected {} phenotype ID: {}".format(modality, phenotype))
    rest = phenotype.split("__", 1)[1]
    if modality != "alt_TSS":
        return rest
    if "_upstream_" not in rest:
        raise ValueError("alt_TSS phenotype lacks _upstream_: " + phenotype)
    return rest.split("_upstream_", 1)[1]


def gene_of(modality, phenotype):
    """Gene ID: the whole phenotype (expression, stability) or the GENE__ prefix.

    Novel PANTRY genes (HPLRG_*) are kept; they simply map to no GUSTO gene.
    """
    gene = phenotype.split("__", 1)[0]
    if not gene:
        raise ValueError("Unexpected {} phenotype ID: {}".format(modality, phenotype))
    return gene


def unversioned(identifier):
    return identifier.split(".", 1)[0]


def collect(pantry, ancestry, modality, qval):
    """Significant leads plus independent signals, keyed by (phenotype, variant)."""
    kept = OrderedDict()
    runs = MODALITIES[modality][0]
    found = [r for r in runs if (pantry / "{}_{}_{}_top.tsv".format(ancestry, modality, r)).exists()]
    if not found:
        raise FileNotFoundError("No PANTRY {}_{}_*_top.tsv in {}".format(ancestry, modality, pantry))
    for run in found:
        for kind in ("independent_top", "top"):
            name = "{}_{}_{}_{}.tsv".format(ancestry, modality, run, kind)
            path = pantry / name
            if not path.exists():
                raise FileNotFoundError(path)
            for row in read_tsv(path):
                if kind == "top":
                    if "qval" not in row:
                        raise ValueError("No qval column in " + name)
                    if row["qval"] in ("", "NA") or float(row["qval"]) >= qval:
                        continue
                key = (row["phenotype_id"], row["variant_id"])
                if key in kept:
                    # Same signal in several tables: keep the first row, fill gaps
                    # (e.g. qval from *_top.tsv, rank from *_independent_top.tsv).
                    for col in KEEP:
                        if kept[key][col] == "NA" and row.get(col) not in (None, "", "NA"):
                            kept[key][col] = row[col]
                    kept[key]["source_tables"] += ";" + name
                    continue
                match = VARIANT.match(row["variant_id"])
                if not match:
                    raise ValueError("Unparseable variant_id {} in {}".format(row["variant_id"], name))
                out = OrderedDict([("modality", modality), ("phenotype_id", row["phenotype_id"]),
                                   ("variant_id", row["variant_id"]), ("chrom", match.group(1)),
                                   ("pos", match.group(2))])
                for col in KEEP:
                    out[col] = row.get(col, "NA") or "NA"
                out["source_tables"] = name
                kept[key] = out
    return list(kept.values())


def load_gusto(event_map, window):
    tx_event, tx_gene, genes = {}, {}, set()
    for row in read_tsv(event_map):
        for col in ("transcript_id", "gene_id", "event_id", "window"):
            if col not in row:
                raise ValueError("event_map lacks column " + col)
        genes.add(row["gene_id"])
        tx_gene[row["transcript_id"]] = row["gene_id"]
        if row["window"] == str(window):
            if row["transcript_id"] in tx_event:
                raise ValueError("Transcript in several window-{} events: {}".format(
                    window, row["transcript_id"]))
            tx_event[row["transcript_id"]] = row["event_id"]
    if not tx_event:
        raise ValueError("No window-{} rows in event_map".format(window))
    by_stable = defaultdict(set)
    for gene in genes:
        by_stable[unversioned(gene)].add(gene)
    return tx_event, tx_gene, genes, by_stable


def map_phenotype(modality, phenotype, gusto):
    tx_event, tx_gene, genes, by_stable = gusto
    gene = gene_of(modality, phenotype)
    transcript = transcript_of(modality, phenotype)
    if (transcript in tx_event and
            unversioned(tx_gene[transcript]) == unversioned(gene)):
        return (tx_gene[transcript], tx_event[transcript], "event", "verified",
                "transcript_TSS_event")
    # TSS-linked phenotype whose transcript is not a GUSTO transcript of the gene.
    note = "transcript_not_in_GUSTO;" if transcript else ""
    if gene in genes:
        return gene, "", "gene", "verified", note + "gene_id_exact"
    stable = by_stable.get(unversioned(gene), set())
    if len(stable) == 1:
        return next(iter(stable)), "", "gene", "verified", note + "gene_stable_id_version_differs"
    if len(stable) > 1:
        return "", "", "gene", "ambiguous", note + "gene_stable_id_multiple"
    return "", "", "gene", "unmapped", note + "gene_absent_from_GUSTO"


def write_tsv(path, rows, fields):
    with open(str(path), "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--pantry-dir", required=True, type=Path, help="PANTRY qtl_results directory")
    p.add_argument("--event-map", required=True, type=Path, help="GUSTO paths.event_map")
    p.add_argument("--out-dir", required=True, type=Path, help="new output directory")
    p.add_argument("--qval", type=float, default=0.05, help="top-lead FDR cutoff (default 0.05)")
    p.add_argument("--window", type=int, default=50, help="TSS event window used by te-tss (default 50)")
    p.add_argument("--ancestries", default="EAS,EUR")
    p.add_argument("--modalities", default=",".join(MODALITIES),
                   help="comma-separated (default: all nine supported)")
    a = p.parse_args(argv)
    modalities = a.modalities.split(",")
    for m in modalities:
        if m not in MODALITIES:
            p.error("Unsupported modality {}; choose from {}".format(m, ",".join(MODALITIES)))
    if a.out_dir.exists():
        p.error("Refusing to overwrite existing --out-dir: {}".format(a.out_dir))
    if not 0 < a.qval <= 1:
        p.error("--qval must be in (0, 1]")

    gusto = load_gusto(a.event_map, a.window)
    final = a.out_dir.resolve()
    out = final.with_name(final.name + ".partial")
    if out.exists():
        p.error("Remove the incomplete previous attempt first: {}".format(out))
    (out / "qtl").mkdir(parents=True)
    manifest, summary, phenotypes = [], [], OrderedDict()
    for ancestry in a.ancestries.split(","):
        for modality in modalities:
            rows = collect(a.pantry_dir, ancestry, modality, a.qval)
            if not rows:
                raise ValueError("No significant rows for {} {}".format(ancestry, modality))
            input_id = "PANTRY_{}_{}_lead".format(ancestry, modality)
            path = out / "qtl" / "{}_{}_lead.tsv".format(ancestry, modality)
            write_tsv(path, rows, list(rows[0].keys()))
            manifest.append(OrderedDict([
                ("input_id", input_id),
                ("path", str(final / "qtl" / path.name)), ("build", "GRCh38"),
                ("kind", "lead"), ("coordinate_system", "one_based"), ("ancestry", ancestry)]))
            outcome = Counter()
            for phenotype in OrderedDict.fromkeys(r["phenotype_id"] for r in rows):
                key = (modality, phenotype)
                if key not in phenotypes:
                    gene, event, level, status, basis = map_phenotype(modality, phenotype, gusto)
                    phenotypes[key] = OrderedDict([
                        ("modality", modality), ("phenotype_id", phenotype),
                        ("gene_id", gene or "NA"), ("event_id", event or "NA"),
                        ("mapping_level", level), ("mapping_status", status),
                        ("mapping_basis", basis),
                        ("link_rule", MODALITIES[modality][1]),
                        ("pantry_transcript_id", transcript_of(modality, phenotype) or "NA")])
                m = phenotypes[key]
                outcome[m["mapping_status"] if m["mapping_status"] != "verified"
                        else "verified_" + m["mapping_level"]] += 1
            summary.append(OrderedDict([
                ("input_id", input_id), ("n_rows", len(rows)),
                ("n_phenotypes", sum(outcome.values())),
                ("verified_event", outcome["verified_event"]),
                ("verified_gene", outcome["verified_gene"]),
                ("ambiguous", outcome["ambiguous"]), ("unmapped", outcome["unmapped"])]))
    write_tsv(out / "qtl_inputs.tsv", manifest, list(manifest[0].keys()))
    pmap = list(phenotypes.values())
    write_tsv(out / "phenotype_map.tsv", pmap, list(pmap[0].keys()))
    write_tsv(out / "summary.tsv", summary, list(summary[0].keys()))
    out.rename(final)
    width = max(len(s["input_id"]) for s in summary)
    print("{:<{w}}  rows  phen  event  gene  ambig  unmapped".format("input_id", w=width))
    for s in summary:
        print("{:<{w}}  {n_rows:>4}  {n_phenotypes:>4}  {verified_event:>5}  {verified_gene:>4}"
              "  {ambiguous:>5}  {unmapped:>8}".format(s["input_id"], w=width, **s))
    print("\nqtl_tables:    {}\nphenotype_map: {}".format(final / "qtl_inputs.tsv",
                                                        final / "phenotype_map.tsv"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
