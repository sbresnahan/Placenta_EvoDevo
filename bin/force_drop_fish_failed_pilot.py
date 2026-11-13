#!/usr/bin/env python3
"""Explicit fish pilot-QC exclusion; backs up active manifests and private config.

Python 3.6+, standard library only. Does not mark any RDS receipt complete.
After applying, rerun pilot --resume so valid existing pilot receipts are reused.
"""
import argparse
import csv
import datetime
import json
import os
from pathlib import Path
import re
import shutil
import tempfile

PREFIX = 'fish_Pprolifica_AB'
EXCLUDED = 'poeProA6xB3_pla8_stage45'
EXPECTED = {
    'poeProB3xA6_pla2_stage40', 'poeProB3xA6_pla9_stage40',
    'poeProB3xA6_pla1_stage40', 'poeProA6xB3_pla5_stage45',
    'poeProA6xB3_pla1_stage40', EXCLUDED
}


def tsv(path):
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = reader.fieldnames
        if not fields:
            raise ValueError('Missing manifest header: ' + str(path))
        rows = list(reader)
    if not rows:
        raise ValueError('Empty active manifest: ' + str(path))
    return fields, rows


def encode(fields, rows):
    from io import StringIO
    buffer = StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter='\t',
                            lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def write_atomic(path, text):
    fd, name = tempfile.mkstemp(prefix='.fish_force_', dir=str(path.parent))
    try:
        os.chmod(name, path.stat().st_mode & 0o777)
        with os.fdopen(fd, 'w') as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, str(path))
    finally:
        if os.path.exists(name):
            os.unlink(name)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohort-dir', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--exclude', required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--preview', action='store_true')
    mode.add_argument('--apply', action='store_true')
    args = parser.parse_args(argv)

    if args.exclude != EXCLUDED:
        raise ValueError('This audited override is restricted to ' + EXCLUDED)
    cohort = args.cohort_dir.resolve(strict=True)
    config = args.config.resolve(strict=True)
    if cohort.name != PREFIX or config.name != PREFIX+'.private.yaml':
        raise ValueError('Only the fish_Pprolifica_AB cohort and private config are supported')
    paths = {kind:cohort/'manifests'/(PREFIX+'_'+kind+'.tsv')
             for kind in ('sources', 'libraries', 'units')}
    originals = {kind:tsv(path) for kind, path in paths.items()}
    ids = {
        'sources': (r['sample_id'] for r in originals['sources'][1]),
        'libraries': (r['sample_id'] for r in originals['libraries'][1]),
        'units': (r['unit_id'] for r in originals['units'][1])
    }
    for name, iterable in ids.items():
        values = list(iterable)
        if len(values) != 6 or len(set(values)) != 6 or set(values) != EXPECTED:
            raise ValueError(name + ': expected exactly the six original F1 records; nothing changed')
    src_rows = originals['sources'][1]
    lib_rows = originals['libraries'][1]
    unit_rows = originals['units'][1]
    lib_by_sample = {r['sample_id']:r for r in lib_rows}
    for r in unit_rows:
        if lib_by_sample[r['unit_id']]['unit_id'] != r['unit_id']:
            raise ValueError('Library-to-unit linkage inconsistent')
    if any(r['cluster_id'] not in ('NA', '') for r in unit_rows):
        raise ValueError('Unexpected cluster identities; manually re-evaluate')
    if any(r['tissue'] != 'pericardial_sac' for r in unit_rows):
        raise ValueError('Unexpected tissue')
    for r in unit_rows:
        if (r['maternal_parent'],r['paternal_parent']) == ('Line_A','Line_B'):
            pass
        elif (r['maternal_parent'],r['paternal_parent']) == ('Line_B','Line_A'):
            pass
        else:
            raise ValueError('Invalid reciprocal parent assignment')
        if not r['unit_id'].endswith('_'+r['stage']):
            raise ValueError('Unit ID contradicts stage')
    assert sum(r['maternal_parent']=='Line_A' for r in unit_rows)==3
    assert sum(r['maternal_parent']=='Line_B' for r in unit_rows)==3
    after = {kind:[r for r in rows if r['unit_id' if kind=='units' else 'sample_id']!=EXCLUDED]
             for kind, (_,rows) in originals.items()}
    if len(after['sources'])!=5 or len(after['libraries'])!=5 or len(after['units'])!=5:
        raise ValueError('Expected to exclude exactly one record from each manifest')
    units = after['units']
    if sum(u['maternal_parent']=='Line_A' for u in units)!=2 or \
       sum(u['maternal_parent']=='Line_B' for u in units)!=3 or \
       {u['stage'] for u in units}!={'stage40','stage45'}:
        raise ValueError('Unexpected five-F1 subset')
    lines = config.read_text().splitlines(keepends=True)
    policy_lines = [i for i,x in enumerate(lines) if re.match(r'^fish_stage_policy\s*:',x)]
    if len(policy_lines)!=1 or lines[policy_lines[0]].strip()!='fish_stage_policy: quantification_only':
        raise ValueError('Expected one unmodified fish_stage_policy: quantification_only line')
    lines[policy_lines[0]]='fish_stage_policy: force_pooled_unadjusted\n'
    output = {paths[k]:encode(originals[k][0],after[k]) for k in paths}
    output[config]=''.join(lines)
    report={
        'cohort':PREFIX, 'excluded_sample':EXCLUDED,
        'reason':'pilot pseudoalignment 57278/1000000 < 0.20',
        'n_pilot_candidates_after':5, 'n_A_mother':2, 'n_B_mother':3,
        'stages_retained':{'Line_A':['stage40','stage45'],
                           'Line_B':['stage40','stage40','stage40']},
        'policy':'force_pooled_unadjusted',
        'limitations':'exploratory stage-unadjusted; litter/cluster independence unverified',
        'pilot_receipts':'must be validated by pilot --resume; not fabricated'
    }
    print(json.dumps(report,indent=2,sort_keys=True))
    if args.preview:
        print('PREVIEW ONLY: no files modified')
        return
    backup = cohort/'manifests'/('fish_force_exclusion_'+datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ'))
    backup.mkdir(mode=0o700, parents=False, exist_ok=False)
    for kind,path in paths.items():
        shutil.copy2(str(path),str(backup/path.name))
    shutil.copy2(str(config),str(backup/config.name))
    (backup/'exclusion_provenance.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    try:
        for path, content in output.items():
            write_atomic(path,content)
    except Exception:
        for kind,path in paths.items():
            shutil.copy2(str(backup/path.name),str(path))
        shutil.copy2(str(backup/config.name),str(config))
        raise
    print('APPLIED. Backup and provenance: '+str(backup))

if __name__=='__main__':
    try:
        main()
    except (OSError,ValueError,KeyError,AssertionError) as exc:
        raise SystemExit('Fish force pilot exclusion failed: '+str(exc))
