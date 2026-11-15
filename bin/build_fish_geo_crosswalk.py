#!/usr/bin/env python3
"""Crosswalk GSE318355/SRP674561 to Dryad placental aliases with GEO evidence.

No SRA ordering assumptions. All output remains include=NO, including exact matches.
Python >= 3.6, standard library only.
"""
import argparse
import collections
import csv
import gzip
import json
import re
from pathlib import Path

TISSUE = 'pla'
ALIAS = re.compile(r'^poePro(?P<mother>[AB]\d+)x(?P<father>[AB]\d+)_pla(?P<id>\d+)_stage(?P<stage>\d+)$')
FIELDS = ('run_accession', 'geo_accession', 'biosample', 'experiment_accession',
          'source_tissue', 'dryad_alias', 'match_status', 'cross',
          'maternal_parent', 'paternal_parent', 'stage', 'unit_id', 'cluster_id',
          'include', 'strand_mode', 'evidence', 'metadata_title')


def dryad_samples(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        rd = csv.reader(f)
        cols = next(rd)
        if cols and cols[0].startswith('Data S') and not any(cols[1:]):
            cols = next(rd)
    samples = sorted(set(x[:-len('_ref_count')] for x in cols if x.endswith('_ref_count')))
    if len(samples) != 11 or any(not ALIAS.fullmatch(s) for s in samples):
        raise ValueError('Unexpected Dryad placenta sample names (expected 11)')
    for s in samples:
        if s+'_alt_count' not in cols:
            raise ValueError('Missing alt count for '+s)
    return samples


def geo_samples(path):
    if path is None:
        return {}
    op = gzip.open if str(path).endswith('.gz') else open
    result = {}
    sample = None
    with op(str(path), 'rt', encoding='utf-8', errors='replace') as f:
        for raw in f:
            line = raw.rstrip('\r\n')
            if line.startswith('^SAMPLE = '):
                sample = line.split('=', 1)[1].strip()
                if not re.fullmatch(r'GSM\d+', sample) or sample in result:
                    raise ValueError('Invalid or duplicate GEO SAMPLE: '+sample)
                result[sample] = {'title': '', 'values': []}
            elif line.startswith('^'):
                sample = None
            elif sample and line.startswith('!Sample_'):
                k, sep, v = line.partition(' = ')
                if not sep:
                    continue
                if k == '!Sample_title':
                    result[sample]['title'] = v
                if k.startswith(('!Sample_title', '!Sample_source_name_ch',
                                 '!Sample_characteristics_ch', '!Sample_description',
                                 '!Sample_relation', '!Sample_supplementary_file')):
                    result[sample]['values'].append(v)
    return result


def build(runinfo, aliases, geo):
    with runinfo.open(newline='', encoding='utf-8-sig') as f:
        entries = list(csv.DictReader(f))
    if len(entries) != 29 or len(set(x['Run'] for x in entries)) != 29:
        raise ValueError('Expected 29 unique SRA runs')
    if any(x['SRA Study'] != 'SRP674561' for x in entries):
        raise ValueError('Unexpected study in SRA RunInfo')
    rows = []
    for x in sorted(entries, key=lambda a: a['Run']):
        tissue = x['tissue'].strip().lower()
        if tissue not in ('placenta', 'embryo', 'follicle'):
            raise ValueError('Invalid SRA tissue: '+tissue)
        if tissue != 'placenta':
            continue
        gsm = x['Sample Name']
        meta = geo.get(gsm)
        vals = ' | '.join(meta['values']) if meta else ''
        candidates = [alias for alias in aliases if re.search(
            r'(?<![a-zA-Z0-9_])'+re.escape(alias)+r'(?![a-zA-Z0-9_])', vals)]
        status = ('metadata_unavailable' if meta is None else
                  'exact_geo_alias_match' if len(candidates) == 1 else
                  'ambiguous_geo_alias' if len(candidates) > 1 else
                  'no_literal_alias')
        r = {key: '' for key in FIELDS}
        r.update(run_accession=x['Run'], geo_accession=gsm,
                 biosample=x['BioSample'], experiment_accession=x['Experiment'],
                 source_tissue=tissue, match_status=status, include='NO',
                 evidence=('GSE318355_SOFT_literal_alias' if status == 'exact_geo_alias_match'
                           else 'requires_independent_metadata_review'),
                 metadata_title=meta['title'] if meta else '')
        if status == 'exact_geo_alias_match':
            m = ALIAS.fullmatch(candidates[0])
            r.update(dryad_alias=candidates[0], cross=m.group('mother')+'x'+m.group('father'),
                     maternal_parent='Line_'+m.group('mother')[0],
                     paternal_parent='Line_'+m.group('father')[0],
                     stage='stage'+m.group('stage'), unit_id=candidates[0])
        rows.append(r)
    if len(rows) != 11:
        raise ValueError('Expected exactly 11 placental SRA runs')
    assigned = [r['dryad_alias'] for r in rows if r['dryad_alias']]
    if len(set(assigned)) != len(assigned):
        raise ValueError('Duplicate Dryad alias matched to multiple SRA runs')
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runinfo', required=True, type=Path)
    p.add_argument('--datas1', required=True, type=Path)
    p.add_argument('--geo-soft', type=Path)
    p.add_argument('--out-dir', required=True, type=Path)
    a = p.parse_args()
    dest = a.out_dir
    out = dest/'fish_placenta_geo_dryad_crosswalk.tsv'
    report = dest/'fish_placenta_geo_dryad_crosswalk_report.json'
    design = dest/'fish_placenta_dryad_design.tsv'
    if any(x.exists() for x in (out, report, design)):
        p.error('Output exists. Use a new versioned --out-dir.')
    aliases = dryad_samples(a.datas1)
    rows = build(a.runinfo, aliases, geo_samples(a.geo_soft))
    dest.mkdir(parents=True, exist_ok=True)
    with out.open('w', newline='') as f:
        w = csv.DictWriter(f, FIELDS, delimiter='\t', lineterminator='\n')
        w.writeheader();w.writerows(rows)
    with design.open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(('dryad_alias', 'cross', 'maternal_parent', 'paternal_parent', 'stage', 'analysis_role'))
        for s in aliases:
            m = ALIAS.fullmatch(s)
            mother, father, stage = m.group('mother'), m.group('father'), m.group('stage')
            w.writerow((s, mother+'x'+father, 'Line_'+mother[0], 'Line_'+father[0],
                        'stage'+stage, 'reciprocal_F1' if mother[0]!=father[0] else 'within_line_control'))
    by_status = dict(collections.Counter(x['match_status'] for x in rows))
    by_design = dict(collections.Counter((ALIAS.fullmatch(s).group('mother')[0]+'x'+
                                         ALIAS.fullmatch(s).group('father')[0]+'_stage'+
                                         ALIAS.fullmatch(s).group('stage')) for s in aliases))
    summary = dict(n_sra_placenta=len(rows), n_dryad_placenta=len(aliases),
                   n_resolved=sum(x['match_status']=='exact_geo_alias_match' for x in rows),
                   status=by_status, dryad_design=by_design,
                   stage_caveat='F1 samples include stages 40 and 45; direction and stage are imbalanced',
                   production_manifest_ready=False,
                   geo_soft_supplied=bool(a.geo_soft))
    report.write_text(json.dumps(summary, indent=2, sort_keys=True)+'\n')
    print(json.dumps(summary, indent=2, sort_keys=True))
    print('All include=NO. See', out)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError) as e:
        raise SystemExit('Fish GEO crosswalk failed: '+str(e))
