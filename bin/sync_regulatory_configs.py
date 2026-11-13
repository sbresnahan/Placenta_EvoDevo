#!/usr/bin/env python3
"""Add missing regulatory defaults to cohort private configs, without rewriting YAML.

Standard-library only. Deliberately edits just the supported indented map keys;
never round-trips or reserializes the user's private path/credential settings.
"""

import argparse
import datetime
import difflib
import os
from pathlib import Path
import re
import shutil
import tempfile

ROOT = Path(__file__).resolve().parent.parent
CONFIG_DIR = ROOT / "config"
COHORT_SPECIES = {
    "maize_B73_Ki11": "maize",
    "maize_B73_W22": "maize",
    "mouse_Babr_CAST": "mouse",
    "mouse_FVB_CAST_E12": "mouse",
    "fish_Pprolifica_AB": "fish",
    "GUSTO": "human",
}
COMMON = {
    "bootstrap": "2000",
    "seed": "1147",
    "min_depth": "10",
    "max_pairs_per_gene": "10",
}
FISH = {
    "runtime": "container",
    "image": "/risapps/apptainer/repo/tetools/20240523/tetools-202405023.sif",
    "executable": "apptainer",
    "seed": "1147",
}
HUMAN = {"require_qtl": "true"}
# Conservative *initial* resource requests for the supported Seadragon queues.
# They may need adjustment after measuring full-genome RepeatModeler demand.
# Existing private settings are never overwritten.
FISH_LSF = {
    "fish_repeats": {
        "queue": "long", "cpus": None,
        "memory_value": "64.0", "reservation_value": "64.0", "wall_time": "'120:00'",
    },
}


class UnsafeConfig(ValueError):
    pass


def section(lines, name, indent, lo=0, hi=None, require_map=False):
    """Return [header, next sibling boundary); fail on duplicate/inline keys."""
    if hi is None:
        hi = len(lines)
    pattern = re.compile(r"^" + " " * indent + re.escape(name) + r":(?:\s*(.*?))?\s*$")
    hits = []
    for i in range(lo, hi):
        m = pattern.match(lines[i].rstrip("\n"))
        if m and not lines[i].lstrip().startswith("#"):
            hits.append((i, m.group(1) or ""))
    if len(hits) > 1:
        raise UnsafeConfig("duplicate YAML key: " + name)
    if not hits:
        return None
    start, value = hits[0]
    if require_map and value and not value.startswith("#"):
        raise UnsafeConfig("inline mapping/scalar not supported for: " + name)
    end = hi
    for i in range(start + 1, hi):
        line = lines[i]
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        spaces = len(line) - len(line.lstrip(" "))
        if "\t" in line[:len(line) - len(line.lstrip())]:
            raise UnsafeConfig("tabs in YAML indentation")
        if spaces <= indent:
            end = i
            break
    return start, end


def add_keys(lines, parent_name, parent_indent, children, parent_range=None):
    """Insert only absent simple scalar keys into a block mapping."""
    low, high = parent_range if parent_range is not None else (0, len(lines))
    match = section(lines, parent_name, parent_indent, low, high, require_map=True)
    if match is None:
        raise UnsafeConfig("missing parent: " + parent_name)
    start, stop = match
    for key, value in children.items():
        match = section(lines, key, parent_indent + 2, start + 1, stop)
        if match is None:
            lines.insert(stop, " " * (parent_indent + 2) + f"{key}: {value}\n")
            stop += 1
    return start, stop


def add_child_block(lines, parent_name, child_name, keys):
    parent = section(lines, parent_name, 0, require_map=True)
    if parent is None:
        raise UnsafeConfig("missing " + parent_name)
    start, end = parent
    match = section(lines, child_name, 2, start + 1, end, require_map=True)
    if match is None:
        additions = [f"  {child_name}:\n"]
        additions += [f"    {k}: {v}\n" for k, v in keys.items()]
        lines[end:end] = additions
    else:
        add_keys(lines, child_name, 2, keys, (start + 1, end))


def fish_cpu_slots(lines):
    """Read the effective existing private regulatory.fish.cpus value."""
    reg = section(lines, "regulatory", 0, require_map=True)
    fish = section(lines, "fish", 2, reg[0] + 1, reg[1], require_map=True)
    if fish is None:
        raise UnsafeConfig("missing regulatory.fish")
    cpu = section(lines, "cpus", 4, fish[0] + 1, fish[1])
    if cpu is None:
        raise UnsafeConfig("missing regulatory.fish.cpus")
    raw = lines[cpu[0]].split(":", 1)[1].split("#", 1)[0].strip()
    if not re.fullmatch(r"[0-9]+", raw) or int(raw) < 4:
        raise UnsafeConfig("regulatory.fish.cpus must be an integer >=4")
    return int(raw)


def add_fish_lsf(lines):
    """Fill missing resource keys in existing lsf.stages; preserve custom values."""
    cpu_slots = fish_cpu_slots(lines)
    for stage, defaults in FISH_LSF.items():
        lsf = section(lines, "lsf", 0, require_map=True)
        if lsf is None:
            raise UnsafeConfig("missing lsf")
        stages = section(lines, "stages", 2, lsf[0] + 1, lsf[1], require_map=True)
        if stages is None:
            raise UnsafeConfig("missing lsf.stages")
        start, end = stages
        values = {k: (str(cpu_slots) if v is None else v)
                  for k, v in defaults.items()}
        current = section(lines, stage, 4, start + 1, end, require_map=True)
        if current is None:
            lines[end:end] = [f"    {stage}:\n"] + [
                f"      {k}: {v}\n" for k, v in values.items()]
        else:
            add_keys(lines, stage, 4, values, (start + 1, end))
        # Do not silently generate a job requesting fewer CPU slots than the
        # configured RepeatModeler/RepeatMasker worker uses.
        lsf = section(lines, "lsf", 0, require_map=True)
        stages = section(lines, "stages", 2, lsf[0] + 1, lsf[1], require_map=True)
        job = section(lines, stage, 4, stages[0] + 1, stages[1], require_map=True)
        cpu = section(lines, "cpus", 6, job[0] + 1, job[1])
        raw = lines[cpu[0]].split(":", 1)[1].split("#", 1)[0].strip()
        if not re.fullmatch(r"[0-9]+", raw) or int(raw) < cpu_slots:
            raise UnsafeConfig(f"lsf.stages.{stage}.cpus must be >= "
                               f"regulatory.fish.cpus ({cpu_slots})")


def config_scalar(lines, key):
    pattern = re.compile(r"^" + re.escape(key) + r":\s*([^#\n]+)")
    values = []
    for line in lines:
        match = pattern.match(line)
        if match:
            values.append(match.group(1).strip().strip("'\""))
    if len(values) != 1:
        raise UnsafeConfig("expected exactly one top-level " + key)
    return values[0]


def update_yaml(data, cohort, cpus):
    lines = data.splitlines(keepends=True)
    if not lines or not lines[-1].endswith("\n"):
        raise UnsafeConfig("YAML file must end with a newline")
    if config_scalar(lines, "schema_version") != "2":
        raise UnsafeConfig("expected schema_version: 2")
    if config_scalar(lines, "cohort_id") != cohort or config_scalar(lines, "species") != COHORT_SPECIES[cohort]:
        raise UnsafeConfig("cohort_id/species mismatch; refusing update")
    if section(lines, "regulatory", 0, require_map=True) is None:
        lines += ["\n", "# Optional regulatory analysis. Does not change primary pipeline stages.\n", "regulatory:\n"]
    add_keys(lines, "regulatory", 0, COMMON)
    if cohort == "fish_Pprolifica_AB":
        add_child_block(lines, "regulatory", "fish", {**FISH, "cpus": str(cpus)})
        add_fish_lsf(lines)
    if cohort == "GUSTO":
        add_child_block(lines, "regulatory", "human_context", HUMAN)
    return "".join(lines)


def update_one(cohort, cpus, dry_run, create):
    example = CONFIG_DIR / f"{cohort}.example.yaml"
    private = CONFIG_DIR / f"{cohort}.private.yaml"
    if not example.is_file():
        raise UnsafeConfig("missing example: " + str(example))
    if private.is_symlink():
        raise UnsafeConfig("refusing symlink private config: " + str(private))
    if private.exists():
        if not private.is_file():
            raise UnsafeConfig("refusing nonregular or symlink private config: " + str(private))
        original = private.read_text(encoding="utf-8")
    elif create:
        original = example.read_text(encoding="utf-8")
    else:
        print(f"SKIP {cohort}: private config does not exist (pass --create-missing)")
        return
    result = update_yaml(original, cohort, cpus)
    if not private.exists() and create and cohort == "fish_Pprolifica_AB":
        # The shipped example is runnable at four slots. On first creation only,
        # honor an explicitly requested larger allocation without modifying
        # existing private CPU choices.
        output = result.splitlines(keepends=True)
        reg_start, reg_end = section(output, "regulatory", 0, require_map=True)
        fish_start, fish_end = section(output, "fish", 2, reg_start + 1, reg_end, require_map=True)
        cpu_start, _ = section(output, "cpus", 4, fish_start + 1, fish_end)
        output[cpu_start] = f"    cpus: {cpus}  # Must not exceed allocated LSF slots.\n"
        # A newly created private config may request >4 CPUs. Keep the LSF
        # stages in sync with that explicit initial choice, without modifying
        # any already-existing private settings.
        for stage in FISH_LSF:
            lsf_start, lsf_end = section(output, "lsf", 0, require_map=True)
            stages_start, stages_end = section(output, "stages", 2,
                                              lsf_start + 1, lsf_end, require_map=True)
            job_start, job_end = section(output, stage, 4,
                                         stages_start + 1, stages_end, require_map=True)
            job_cpu, _ = section(output, "cpus", 6, job_start + 1, job_end)
            output[job_cpu] = f"      cpus: {cpus}\n"
        result = "".join(output)
    if result == original and private.exists():
        print(f"UNCHANGED {private}")
        return
    if dry_run:
        print("".join(difflib.unified_diff(
            original.splitlines(keepends=True), result.splitlines(keepends=True),
            fromfile=str(private) + " (before)", tofile=str(private) + " (after)")), end="")
        return
    # Atomic write, backup, and existing chmod: preserve private edits and secrets.
    if private.exists():
        stamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
        backup = Path(str(private) + ".bak." + stamp)
        if backup.exists():
            raise UnsafeConfig("backup collision: " + str(backup))
        shutil.copy2(private, backup)
        print(f"BACKUP {backup}")
    mode = private.stat().st_mode & 0o777 if private.exists() else 0o600
    fd, tmp = tempfile.mkstemp(prefix=f".{cohort}.", suffix=".tmp", dir=str(CONFIG_DIR))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as out:
            out.write(result)
        os.chmod(tmp, mode)
        os.replace(tmp, private)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    print(f"UPDATED {private}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cohorts", nargs="*", help="cohort IDs; default all six")
    p.add_argument("--fish-cpus", type=int, default=4, metavar="N",
                   help="fish repeat CPU slots, >=4 (default: 4); existing settings are preserved")
    p.add_argument("--dry-run", action="store_true", help="display changes without writing")
    p.add_argument("--create-missing", action="store_true",
                   help="create absent private configs from shipped examples")
    args = p.parse_args()
    if args.fish_cpus < 4:
        p.error("--fish-cpus must be >=4")
    cohorts = args.cohorts or list(COHORT_SPECIES)
    try:
        for cohort in cohorts:
            if cohort not in COHORT_SPECIES:
                raise UnsafeConfig("unknown cohort " + cohort)
            update_one(cohort, args.fish_cpus, args.dry_run, args.create_missing)
    except (UnsafeConfig, UnicodeError, OSError) as exc:
        p.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
