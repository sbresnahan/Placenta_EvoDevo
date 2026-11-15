#!/usr/bin/env python3
"""Conservative cross-cohort POE analysis from *aggregate* gene estimates.

No individual or trio genotypes enter this program. The script never changes
primary model outputs. Cross-species output is descriptive, not a meta-analysis
across developmental stages, tissues, ploidies or promoter architectures.
"""
import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


def rows(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, data, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ValueError(f"Refusing to overwrite existing analysis: {path}")
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data)


def require(record, names, context):
    missing = [key for key in names if key not in record]
    if missing:
        raise ValueError(f"{context}: missing {', '.join(missing)}")


def number(x):
    try:
        n = float(x)
    except (TypeError, ValueError):
        return None
    return n if math.isfinite(n) else None


def bh(values):
    result = [None] * len(values)
    ranked = sorted((p, i) for i, p in enumerate(values) if p is not None)
    bound = 1.0
    for rev_index in range(len(ranked) - 1, -1, -1):
        p, idx = ranked[rev_index]
        bound = min(bound, p * len(ranked) / (rev_index + 1))
        result[idx] = bound
    return result


def gamma_q(a, x):
    """Regularized upper incomplete gamma, for chi-square survival tests."""
    if x <= 0:
        return 1.0
    if x < a + 1:
        total = term = 1 / a
        ap = a
        for _ in range(10000):
            ap += 1
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-14:
                break
        return max(0.0, min(1.0, 1 - total * math.exp(-x + a * math.log(x) - math.lgamma(a))))
    b = x + 1 - a
    c = 1e300
    d = 1 / max(abs(b), 1e-300) * (-1 if b < 0 else 1)
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        if abs(d) < 1e-300:
            d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 1e-14:
            break
    return max(0.0, min(1.0, math.exp(-x + a * math.log(x) - math.lgamma(a)) * h))


def fixed_effect(data):
    weight = [1 / (d["se_ln"] ** 2) for d in data]
    total = sum(weight)
    estimate = sum(w * d["effect_ln"] for w, d in zip(weight, data)) / total
    se = math.sqrt(1 / total)
    q = sum(w * (d["effect_ln"] - estimate) ** 2 for w, d in zip(weight, data))
    df = len(data) - 1
    return dict(effect_ln=estimate, se_ln=se, ci_low_ln=estimate - 1.96 * se,
                ci_high_ln=estimate + 1.96 * se,
                p_value=math.erfc(abs(estimate / se) / math.sqrt(2)),
                heterogeneity_Q=q, heterogeneity_df=df,
                heterogeneity_p=gamma_q(df / 2, q / 2) if df else None,
                I2=max(0.0, (q - df) / q) if q else 0.0)


def run(manifest, mapping, outdir):
    specs = rows(manifest)
    if not specs:
        raise ValueError("Empty cohort manifest")
    needed = ("cohort_id", "species", "tissue", "stage", "results_path",
              "effect_scale", "independence_group")
    for spec in specs:
        require(spec, needed, "manifest")
        if spec["effect_scale"] not in ("log2", "ln"):
            raise ValueError("effect_scale must be log2 or ln")
        if not all(spec[k].strip() for k in needed):
            raise ValueError("Empty required manifest value")
        if spec["species"] == "human":
            if spec.get("approved_release", "").lower() != "true":
                raise ValueError("Human aggregate results require explicit approved_release=true")
            min_cell = number(spec.get("minimum_cell_size"))
            if min_cell is None or min_cell < 10 or int(min_cell) != min_cell:
                raise ValueError("Human aggregate results require an approved integer minimum_cell_size >= 10")
    if len(set(s["cohort_id"] for s in specs)) != len(specs):
        raise ValueError("Each cohort can be represented once per run")

    ortholog = rows(mapping)
    if not ortholog:
        raise ValueError("Orthology table is required; never infer gene matching from symbols")
    genes = defaultdict(set)
    for r in ortholog:
        require(r, ("species", "gene_id", "orthogroup_id"), "ortholog")
        if not all(r[k] for k in ("species", "gene_id", "orthogroup_id")):
            raise ValueError("Blank orthology identity")
        genes[(r["species"], r["gene_id"])].add(r["orthogroup_id"])
    if any(len(v) != 1 for v in genes.values()):
        raise ValueError("A gene maps to multiple orthogroups; curate mapping first")
    orthogroup = {key: next(iter(value)) for key, value in genes.items()}
    paralogs = Counter((species, group) for (species, _), group in orthogroup.items())

    included, excluded = [], []
    for spec in specs:
        p = Path(spec["results_path"])
        result_rows = rows(p)
        if len({d.get("gene_id") for d in result_rows}) != len(result_rows):
            raise ValueError(f"One and only one gene POE test per gene is required: {p}")
        for d in result_rows:
            require(d, ("gene_id", "status", "estimate", "se"), str(p))
            if d.get("effect_type", "gene_parent_of_origin") != "gene_parent_of_origin":
                raise ValueError(f"Comparator accepts only gene-level parent-of-origin tests: {p}")
            record = {"cohort_id": spec["cohort_id"], "species": spec["species"],
                      "tissue": spec["tissue"], "stage": spec["stage"],
                      "gene_id": d["gene_id"], "independence_group": spec["independence_group"]}
            group = orthogroup.get((spec["species"], d["gene_id"]))
            effect = number(d["estimate"])
            se = number(d["se"])
            cell = number(d.get("n_families")) if spec["species"] == "human" else None
            if spec["species"] == "human" and (cell is None or cell < float(spec["minimum_cell_size"])):
                reason = "human_aggregate_below_approved_cell_floor"
            elif d["status"] != "eligible_conditional_inference":
                reason = "primary_model_ineligible"
            elif group is None:
                reason = "no_curated_ortholog"
            elif paralogs[(spec["species"], group)] != 1:
                reason = "within_species_paralogy_not_one_to_one"
            elif effect is None or se is None or se <= 0:
                reason = "invalid_estimate_or_se"
            else:
                reason = None
            if reason:
                excluded.append({**record, "reason": reason})
            else:
                factor = math.log(2) if spec["effect_scale"] == "log2" else 1.0
                included.append({**record, "orthogroup_id": group,
                                 "effect_ln": effect * factor, "se_ln": se * factor})

    # Comparison strata respect species, tissue and stage. In particular,
    # maize 14-DAP W22 and 20-DAP Ki11 must not be pooled as if exchangeable.
    strata = defaultdict(list)
    for d in included:
        strata[(d["species"], d["tissue"], d["stage"], d["orthogroup_id"])].append(d)
    summaries = []
    for (species, tissue, stage, group), data in sorted(strata.items()):
        unique_groups = len({d["independence_group"] for d in data})
        if len(data) > 1 and unique_groups != len(data):
            state, result = "dependent_cohorts_no_pool", fixed_effect(data[:1])
        elif len(data) == 1:
            state, result = "single_cohort_descriptive", fixed_effect(data)
        else:
            state, result = "within_stratum_fixed_effect", fixed_effect(data)
        summaries.append(dict(species=species, tissue=tissue, stage=stage,
                              orthogroup_id=group, cohorts=";".join(sorted(d["cohort_id"] for d in data)),
                              n_cohorts=len(data), status=state, **result))
    sig = bh([r["p_value"] if r["status"] == "within_stratum_fixed_effect" else None
              for r in summaries])
    for d, q in zip(summaries, sig):
        d["q_value"] = q

    # Descriptive stage- and species-level sign agreement; no pooled p values.
    cross = defaultdict(list)
    for d in summaries:
        cross[d["orthogroup_id"]].append(d)
    patterns = []
    for group, data in sorted(cross.items()):
        species_set = sorted({r["species"] for r in data})
        signs = [1 if r["effect_ln"] > 0 else -1 if r["effect_ln"] < 0 else 0 for r in data]
        patterns.append(dict(orthogroup_id=group, n_strata=len(data),
                             n_species=len(species_set), species=";".join(species_set),
                             strata=";".join(f"{r['species']}:{r['tissue']}:{r['stage']}" for r in data),
                             direction_concordant=(all(s == signs[0] for s in signs) and signs[0] != 0),
                             interpretation="descriptive_context_concordance_no_cross_species_p_value"))

    outdir = Path(outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError(f"Refusing to write into nonempty cross-cohort directory: {outdir}")
    write(outdir / "included_gene_effects.tsv", included,
          ("cohort_id", "species", "tissue", "stage", "gene_id", "orthogroup_id", "effect_ln", "se_ln", "independence_group"))
    write(outdir / "excluded_gene_effects.tsv", excluded,
          ("cohort_id", "species", "tissue", "stage", "gene_id", "reason"))
    write(outdir / "within_stratum_meta.tsv", summaries,
          ("species", "tissue", "stage", "orthogroup_id", "cohorts", "n_cohorts", "status",
           "effect_ln", "se_ln", "ci_low_ln", "ci_high_ln", "p_value", "q_value",
           "heterogeneity_Q", "heterogeneity_df", "heterogeneity_p", "I2"))
    write(outdir / "directional_conservation.tsv", patterns,
          ("orthogroup_id", "n_strata", "n_species", "species", "strata", "direction_concordant", "interpretation"))
    print(f"Analyzed {len(included)} aggregate gene estimates; {len(excluded)} excluded; {len(summaries)} strata")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--orthologs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.manifest, args.orthologs, args.out)


if __name__ == "__main__":
    main()
