#!/usr/bin/env python3
"""Audit Dryad SCAFFOLD contigs against a local FASTA index and assembly report.

Only literal aliases documented by NCBI assembly_report are accepted. Generates
contig-map for prepare_fish_published_like_markers.py; never writes active markers.
Python >= 3.6, standard library only.
"""
import argparse
import collections
import csv
import json
from pathlib import Path


def read_fai(path):
    out = {}
    with path.open() as f:
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 2 or not cols[1].isdigit() or cols[0] in out:
                raise ValueError('Invalid FASTA index')
            out[cols[0]] = int(cols[1])
    return out


def assembly_report(path):
    if not path:
        return []
    out = []
    header = None
    with path.open() as f:
        for line in f:
            if line.startswith('# Sequence-Name\t'):
                header = line[2:].rstrip('\n').split('\t')
            elif line.startswith('#') or not line.strip():
                continue
            else:
                if header is None:
                    raise ValueError('Missing NCBI assembly-report header')
                row = dict(zip(header, line.rstrip('\n').split('\t')))
                if 'Sequence-Length' not in row or not row['Sequence-Length'].isdigit():
                    raise ValueError('Invalid assembly report length')
                out.append(row)
    return out


def candidate_contigs(path):
    out = collections.defaultdict(lambda: [0, 0])
    with path.open(newline='') as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for k in ('chrom', 'genomic_pos', 'ref', 'alt'):
            if k not in (rdr.fieldnames or []):
                raise ValueError('Missing candidate column: '+k)
        for r in rdr:
            p = int(r['genomic_pos'])
            if p < 1:
                raise ValueError('Invalid candidate position')
            v = out[r['chrom']]
            v[0] += 1
            v[1] = max(v[1], p)
    if not out:
        raise ValueError('Empty candidate panel')
    return out


def resolve(candidates, fai, report):
    indexes = collections.defaultdict(list)
    fields = ('Sequence-Name', 'GenBank-Accn', 'RefSeq-Accn', 'UCSC-style-name')
    for r in report:
        for k in fields:
            s = r.get(k, '')
            if s and s != 'na':
                indexes[s].append(r)
    rows, aliases = [], []
    used = {}
    for dryad, (sites, maximum) in sorted(candidates.items()):
        alias, reason, length = '', 'unresolved_contig_alias', ''
        if dryad in fai:
            alias, reason, length = dryad, 'direct_fasta_name', fai[dryad]
        else:
            hit = indexes.get(dryad, [])
            match = []
            for r in hit:
                targets = set(r.get(k, '') for k in fields)
                targets.discard('na')
                targets.discard('')
                for target in targets:
                    if target in fai and fai[target] == int(r['Sequence-Length']):
                        match.append((target, r))
            match = list({m[0]: m for m in match}.values())
            if len(hit) == 1 and len(match) == 1:
                alias, reason, length = match[0][0], 'assembly_report_len_verified', fai[match[0][0]]
            elif hit and not match:
                reason = 'assembly_report_target_or_length_mismatch'
            elif len(match) > 1 or len(hit) > 1:
                reason = 'ambiguous_assembly_report_alias'
        if alias and maximum > length:
            reason = 'candidate_position_out_of_bounds'
            alias = ''
        if alias and alias in used and used[alias] != dryad:
            reason = 'two_dryad_contigs_map_one_fasta'
            alias = ''
        if alias:
            used[alias] = dryad
            aliases.append(dict(dryad_scaffold=dryad, fasta_contig=alias))
        rows.append(dict(dryad_scaffold=dryad, fasta_contig=alias,
                         candidate_sites=sites, max_candidate_position=maximum,
                         fasta_contig_length=length, status=reason))
    return rows, aliases


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidates', required=True, type=Path)
    p.add_argument('--genome-fai', required=True, type=Path)
    p.add_argument('--assembly-report', type=Path)
    p.add_argument('--out-dir', required=True, type=Path)
    a = p.parse_args()
    dest = a.out_dir
    audit = dest/'fish_contig_alias_audit.tsv'
    aliases = dest/'verified_contig_aliases.tsv'
    report = dest/'fish_contig_alias_report.json'
    if any(x.exists() for x in (audit, aliases, report)):
        p.error('Output exists; choose a new versioned --out-dir.')
    sites = candidate_contigs(a.candidates)
    fidx = read_fai(a.genome_fai)
    rows, mapping = resolve(sites, fidx, assembly_report(a.assembly_report))
    dest.mkdir(parents=True, exist_ok=True)
    with audit.open('w', newline='') as f:
        w = csv.DictWriter(f, ('dryad_scaffold','fasta_contig','candidate_sites',
                                'max_candidate_position','fasta_contig_length','status'),
                           delimiter='\t', lineterminator='\n')
        w.writeheader();w.writerows(rows)
    with aliases.open('w', newline='') as f:
        w = csv.DictWriter(f, ('dryad_scaffold','fasta_contig'), delimiter='\t', lineterminator='\n')
        w.writeheader();w.writerows(mapping)
    n = sum(r['candidate_sites'] for r in rows)
    good = sum(r['candidate_sites'] for r in rows if r['fasta_contig'])
    summary = dict(n_candidate_sites=n, n_contigs=len(rows),
                   n_contigs_mapped=len(mapping), n_sites_mapped=good,
                   fraction_sites_mapped=good/n,
                   status_counts=dict(collections.Counter(r['status'] for r in rows)),
                   ref_base_validation_performed=False,
                   next_step='run prepare_fish_published_like_markers.py --genome with verified_contig_aliases.tsv')
    report.write_text(json.dumps(summary, indent=2, sort_keys=True)+'\n')
    print(json.dumps(summary, indent=2, sort_keys=True))
    if good < n:
        print('WARNING: some candidate contigs are unresolved; do not promote production markers')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError) as err:
        raise SystemExit('Fish assembly audit failed: '+str(err))
