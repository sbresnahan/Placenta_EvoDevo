#!/usr/bin/env python3
"""Synthetic schema/algebra fixtures; not patient or cohort data."""
import csv
import importlib.util
import math
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("compare_cohorts", ROOT / "bin/compare_cohorts.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def tsv(path, data):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=data[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(data)


class CrossCohortTests(unittest.TestCase):
    def test_bh_and_meta(self):
        self.assertEqual(module.bh([0.01, 0.04, None]), [0.02, 0.04, None])
        z = module.fixed_effect([dict(effect_ln=0.3, se_ln=0.2),
                                 dict(effect_ln=0.3, se_ln=0.2)])
        self.assertAlmostEqual(z["effect_ln"], 0.3)
        self.assertEqual(z["heterogeneity_Q"], 0)
        self.assertAlmostEqual(z["heterogeneity_p"], 1)
        self.assertAlmostEqual(module.gamma_q(0.5, 0.5), math.erfc(math.sqrt(0.5)), places=12)

    def test_separate_stages_and_log_scales(self):
        with tempfile.TemporaryDirectory() as d:
            base = Path(d)
            tsv(base / "genes.tsv", [dict(gene_id="G1", status="eligible_conditional_inference",
                                          estimate="1", se="0.5")])
            tsv(base / "orthologs.tsv", [dict(species="maize", gene_id="G1", orthogroup_id="O1")])
            tsv(base / "manifest.tsv", [dict(cohort_id="ki11", species="maize", tissue="endosperm",
                stage="20_DAP", results_path=str(base / "genes.tsv"), effect_scale="log2",
                independence_group="biological_ki11"),
                dict(cohort_id="w22", species="maize", tissue="endosperm",
                stage="14_DAP", results_path=str(base / "genes.tsv"), effect_scale="log2",
                independence_group="biological_w22")])
            module.run(base / "manifest.tsv", base / "orthologs.tsv", base / "out")
            summary = module.rows(base / "out/within_stratum_meta.tsv")
            self.assertEqual(len(summary), 2)
            self.assertEqual({x["status"] for x in summary}, {"single_cohort_descriptive"})
            included = module.rows(base / "out/included_gene_effects.tsv")
            self.assertAlmostEqual(float(included[0]["effect_ln"]), math.log(2))
            with self.assertRaises(ValueError):
                module.run(base / "manifest.tsv", base / "orthologs.tsv", base / "out")

    def test_human_requires_approved_aggregate(self):
        with tempfile.TemporaryDirectory() as d:
            base = Path(d)
            tsv(base / "gene.tsv", [dict(gene_id="H1", status="eligible_conditional_inference",
                estimate="0.4", se="0.2", n_families="8")])
            tsv(base / "orthologs.tsv", [dict(species="human", gene_id="H1", orthogroup_id="OH")])
            spec = dict(cohort_id="GUSTO", species="human", tissue="placenta",
                stage="term", results_path=str(base / "gene.tsv"), effect_scale="ln",
                independence_group="families_gusto", approved_release="false", minimum_cell_size="10")
            tsv(base / "manifest.tsv", [spec])
            with self.assertRaisesRegex(ValueError, "approved_release"):
                module.run(base / "manifest.tsv", base / "orthologs.tsv", base / "out")
            spec["approved_release"] = "true"
            tsv(base / "manifest.tsv", [spec])
            module.run(base / "manifest.tsv", base / "orthologs.tsv", base / "out")
            self.assertEqual(module.rows(base / "out/included_gene_effects.tsv"), [])
            self.assertEqual(module.rows(base / "out/excluded_gene_effects.tsv")[0]["reason"],
                             "human_aggregate_below_approved_cell_floor")


if __name__ == "__main__":
    unittest.main()
