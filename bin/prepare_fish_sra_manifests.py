#!/usr/bin/env python3
"""Prepare F1-only public-cross manifests from the verified fish GEO/Dryad crosswalk.

This script does NOT change the active cohort manifests. It verifies local SRA
receipts, records all 11 placental runs in an audit, and stages exactly six
reciprocal F1 units for subsequent review and activation. Python >=3.6.
"""
import argparse
import csv
import json
import re
from pathlib import Path

CROSS_FIELDS = ('run_accession','geo_accession','biosample','experiment_accession',
    'source_tissue','dryad_alias','match_status','maternal_parent',
    'paternal_parent','stage','unit_id')
SOURCE_FIELDS = ('run_id','sample_id','study_accession','sample_accession',
    'experiment_accession','library_name','layout','read_count','base_count',
    'url_1','md5_1','bytes_1','url_2','md5_2','bytes_2')
LIB_FIELDS = ('sample_id','unit_id','layout','strand_mode')
UNIT_FIELDS = ('unit_id','cluster_id','maternal_parent','paternal_parent','stage','tissue')
AUDIT_FIELDS = ('run_accession','geo_accession','dryad_alias','maternal_parent',
    'paternal_parent','stage','analysis_role','in_f1_manifests',
    'fastq_mates_receipt_verified','strand_mode_evidence')
RUN_ID = re.compile(r'^SRR[0-9]+$')
SAFE = re.compile(r'^[A-Za-z0-9_.-]+$')
ALIAS = re.compile(r'^poePro(?P<m>[AB][0-9]+)x(?P<p>[AB][0-9]+)_pla[0-9]+_stage(?P<s>[0-9]+)$')


def read(path, sep):
    with path.open(newline='', encoding='utf-8-sig') as f:
        table = csv.DictReader(f, delimiter=sep)
        if not table.fieldnames:
            raise ValueError('Missing header: {}'.format(path))
        return list(table), set(table.fieldnames)


def check_fields(cols, required, desc):
    missing = set(required) - cols
    if missing:
        raise ValueError('{} missing columns: {}'.format(desc, ', '.join(sorted(missing))))


def positive_integer(value, label):
    if not re.fullmatch(r'[1-9][0-9]*', str(value or '')):
        raise ValueError('Invalid {}: {!r}'.format(label, value))
    return int(value)


def validate_receipts(raw, run):
    files = [raw / '{}_{}.fastq.gz'.format(run, mate) for mate in (1,2)]
    receipt = raw / '{}_1.fastq.gz.sra.tsv'.format(run)
    if not receipt.is_file():
        raise ValueError('{}: missing SRA receipt {}'.format(run, receipt))
    rows, fields = read(receipt, '\t')
    check_fields(fields, ('run_id','mate','path','bytes','transport'), 'SRA receipt')
    if len(rows) != 2 or {x['mate'] for x in rows} != {'1','2'}:
        raise ValueError('{}: expected exactly two mate receipts'.format(run))
    by_mate = {int(row['mate']): row for row in rows}
    sizes = []
    for mate in (1,2):
        item = by_mate[mate]
        path = files[mate-1]
        size = positive_integer(item['bytes'], 'SRA receipt byte count')
        if item['run_id'] != run or item['transport'] != 'sra' or item['path'] != str(path):
            raise ValueError('{}: receipt/run/transport/path mismatch'.format(run))
        if not path.is_file() or path.stat().st_size != size:
            raise ValueError('{}: FASTQ missing or file size differs from receipt'.format(run))
        sizes.append(size)
    return sizes


def write(path, names, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=names, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def build(crosswalk, runinfo, raw_dir, output_dir, strand_mode):
    if strand_mode not in ('fr','rf','unstranded'):
        raise ValueError('Explicit strand mode is required')
    cross, cols = read(crosswalk, '\t')
    check_fields(cols, CROSS_FIELDS, 'Crosswalk')
    if len(cross) != 11 or len({r['run_accession'] for r in cross}) != 11:
        raise ValueError('Expected 11 distinct placental crosswalk records')
    runs, cols = read(runinfo, ',')
    check_fields(cols, ('Run','Sample Name','BioSample','Experiment','LibraryLayout',
                        'tissue','Bases','AvgSpotLen','SRA Study'), 'SRA RunInfo')
    if len(runs) != 29 or len({r['Run'] for r in runs}) != 29:
        raise ValueError('Expected 29 unique SRA RunInfo records')
    by_run = {r['Run']:r for r in runs}
    if sum(x['tissue'].strip().lower() == 'placenta' for x in runs) != 11:
        raise ValueError('Expected exactly 11 placental SRA runs')
    if {r['run_accession'] for r in cross} != {
        x['Run'] for x in runs if x['tissue'].strip().lower() == 'placenta'}:
        raise ValueError('RunInfo/crosswalk placenta accession mismatch')
    aliases, sources, libraries, units, audit = set(), [], [], [], []
    for r in sorted(cross, key=lambda x: x['run_accession']):
        run = r['run_accession']; entry = by_run[run]
        if not RUN_ID.fullmatch(run) or r['match_status'] != 'exact_geo_alias_match':
            raise ValueError('{}: unresolved or unsafe SRA-to-Dryad match'.format(run))
        alias = r['dryad_alias']; parsed = ALIAS.fullmatch(alias)
        if parsed is None or not SAFE.fullmatch(alias) or alias in aliases or r['unit_id'] != alias:
            raise ValueError('{}: invalid/duplicate Dryad alias or unit ID'.format(run))
        aliases.add(alias)
        maternal = 'Line_' + parsed.group('m')[0]
        paternal = 'Line_' + parsed.group('p')[0]
        stage = 'stage'+parsed.group('s')
        if (r['maternal_parent'],r['paternal_parent'],r['stage']) != (maternal,paternal,stage):
            raise ValueError('{}: crosswalk parents/stage contradict Dryad alias'.format(run))
        if (r['geo_accession'] != entry['Sample Name'] or
            r['biosample'] != entry['BioSample'] or
            r['experiment_accession'] != entry['Experiment'] or
            r['source_tissue'] != 'placenta' or entry['tissue'] != 'placenta' or
            entry['SRA Study'] != 'SRP674561' or entry['LibraryLayout'] != 'PAIRED'):
            raise ValueError('{}: GEO/SRA accession or layout mismatch'.format(run))
        bases = positive_integer(entry['Bases'], run+' bases')
        avg_len = positive_integer(entry['AvgSpotLen'], run+' AvgSpotLen')
        if bases % avg_len:
            raise ValueError('{}: SRA fragment estimate is not integral'.format(run))
        fragments = bases // avg_len
        role = 'reciprocal_F1' if maternal != paternal else 'within_line_control'
        selected = role == 'reciprocal_F1'
        # Verify all 11 before deciding whether to activate any F1s.
        sizes = validate_receipts(raw_dir, run)
        audit.append(dict(run_accession=run,geo_accession=r['geo_accession'],
            dryad_alias=alias,maternal_parent=maternal,paternal_parent=paternal,
            stage=stage,analysis_role=role,in_f1_manifests='YES' if selected else 'NO',
            fastq_mates_receipt_verified='YES',strand_mode_evidence='provisional_pilot_needed'))
        if not selected:
            continue
        # Source URLs are SRA accession pages, NOT ENA FASTQ links. The
        # retrieval.method=sra worker uses run_id instead; bytes are the
        # existing compressed FASTQ sizes from SRA receipts.
        accession_url = 'https://www.ncbi.nlm.nih.gov/sra/'+run
        sources.append(dict(run_id=run,sample_id=alias,study_accession='SRP674561',
            sample_accession=r['biosample'],experiment_accession=r['experiment_accession'],
            library_name=r['geo_accession'],layout='PE',read_count=str(2*fragments),
            base_count=str(bases),url_1=accession_url,md5_1='NA',bytes_1=str(sizes[0]),
            url_2=accession_url,md5_2='NA',bytes_2=str(sizes[1])))
        libraries.append(dict(sample_id=alias,unit_id=alias,layout='PE',strand_mode=strand_mode))
        # Biological cluster identity (shared dam/litter) is unverified.
        # Missing cluster IDs must NOT be replaced with unique synthetic IDs.
        units.append(dict(unit_id=alias,cluster_id='NA',maternal_parent=maternal,
                          paternal_parent=paternal,stage=stage,tissue='pericardial_sac'))
    if len(audit) != 11 or len(sources) != 6 or len(libraries) != 6 or len(units) != 6:
        raise ValueError('Expected 11 audited placenta samples and six F1 samples')
    design = {(u['maternal_parent'],u['paternal_parent'],u['stage']) for u in units}
    if sum(u['maternal_parent']=='Line_A' for u in units)!=3 or sum(
            u['maternal_parent']=='Line_B' for u in units)!=3:
        raise ValueError('Require three F1 biological units per reciprocal cross')
    if {u['stage'] for u in units} != {'stage40','stage45'}:
        raise ValueError('Unexpected F1 developmental-stage design')
    planned = ['fish_Pprolifica_AB_'+v+'.tsv' for v in ('sources','libraries','units')]
    planned += ['fish_sra_manifest_audit.tsv','fish_sra_manifest_report.json']
    if any((output_dir / name).exists() for name in planned):
        raise ValueError('Output already exists; choose a new versioned --output-dir')
    output_dir.mkdir(parents=True, exist_ok=True)
    write(output_dir / planned[0], SOURCE_FIELDS, sources)
    write(output_dir / planned[1], LIB_FIELDS, libraries)
    write(output_dir / planned[2], UNIT_FIELDS, units)
    write(output_dir / planned[3], AUDIT_FIELDS, audit)
    result = dict(n_sra_placenta=11,n_verified_sra_receipts=11,
        n_f1_units=6,n_within_line_controls=5,n_cross_A_mother=3,n_cross_B_mother=3,
        stage_by_maternal_line={line:{stage:sum(u['maternal_parent']==line and u['stage']==stage for u in units)
             for stage in ('stage40','stage45')} for line in ('Line_A','Line_B')},
        strand_mode=strand_mode,strand_mode_evidence='unverified_until_pilot',
        cluster_id_status='unknown_do_not_assume_independent',
        fasta_marker_orientation='Line_A_REF_Line_B_ALT',
        sample_manifest_ready_for_reference=False,
        blocker='mixed stages require explicit quantification-only gate; POE models must remain blocked until stage-aware model and cluster policy are reviewed',
        source_url_note='NCBI SRA accession pages, not direct FASTQ URLs; local byte sizes from SRA receipts')
    (output_dir/planned[4]).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--crosswalk',type=Path,required=True)
    p.add_argument('--runinfo',type=Path,required=True)
    p.add_argument('--raw-dir',type=Path,required=True)
    p.add_argument('--strand-mode',required=True,choices=('fr','rf','unstranded'),
        help='Provisional setting; must be checked by the strandness pilot')
    p.add_argument('--out-dir',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(build(a.crosswalk,a.runinfo,a.raw_dir,a.out_dir,a.strand_mode),indent=2,sort_keys=True))

if __name__ == '__main__':
    try:
        main()
    except (OSError,ValueError,KeyError) as exc:
        raise SystemExit('Fish SRA manifest builder failed: '+str(exc))
