#!/usr/bin/env python3
"""Approximate Zhang et al.'s placental fixed Line A/B SNP panel from Dryad RNA.

The publication reports 55,224 placental SNPs with fixed line differences.
This reconstruction uses DataS1 (placenta), DataS2 (embryo) and DataS3
(maternal follicle), pooling different tissues from the SAME biological fish.
A callable within-line fish needs >=6 combined reads and ZERO opposing
allele reads. At least one Line A fish must support REF and one Line B fish
must support ALT, with no conflicting callable fish. F1 allelic expression is
never used to determine the inherited allele. All selected variants are
RNA-supported candidate markers, not independently verified germline calls.

Python 3.6+, standard library only. Does not overwrite previous outputs.
"""
import argparse
import collections
import csv
import hashlib
import json
import re
from pathlib import Path

CONTROL = re.compile(r'^poePro(A6xA6|B3xB3)_(pla|emb|fol)(\d+)_stage(\d+)$')
UNITS = ('A6xA6_1', 'A6xA6_2', 'A6xA6_7', 'B3xB3_4', 'B3xB3_11')
BASE_FIELDS = ('chrom', 'genomic_pos', 'ref', 'alt', 'focal_index', 'alternate_index')
FIELDS = BASE_FIELDS + ('snp_id', 'gene_id', 'annotation', 'n_line_a_controls',
                        'n_line_b_controls', 'min_control_depth', 'evidence_tier')
CORE = ('SNP_ID', 'Scaffold', 'Position', 'Ref_allele', 'Alt_allele', 'Annotation')


def table(path, tissue):
    f = open(str(path), newline='', encoding='utf-8-sig')
    rdr = csv.reader(f)
    h = next(rdr)
    if h and h[0].startswith('Data S') and not any(h[1:]):
        h = next(rdr)
    if len(h) != len(set(h)) or any(c not in h for c in CORE):
        f.close()
        raise ValueError('Invalid/duplicate Dryad headers: '+str(path))
    controls = []
    for i, col in enumerate(h):
        if not col.endswith('_ref_count'):
            continue
        sample = col[:-10]
        match = CONTROL.fullmatch(sample)
        if not match:
            continue  # Never use reciprocal hybrids to infer marker genotypes.
        group, t, ident, stage = match.groups()
        unit = group+'_'+ident
        if t != tissue or unit not in UNITS or sample+'_alt_count' not in h:
            f.close()
            raise ValueError('Invalid within-line sample: '+sample)
        controls.append((UNITS.index(unit), i, h.index(sample+'_alt_count'), stage))
    if not controls:
        f.close()
        raise ValueError('No within-line controls in '+str(path))
    def records():
        try:
            for row in rdr:
                if not row:
                    continue
                if len(row) != len(h):
                    raise ValueError('Invalid Dryad row width in '+str(path))
                fields = dict(zip(h, row))
                vals = []
                for idx, ci, cj, _ in controls:
                    a, b = row[ci].strip(), row[cj].strip()
                    vals.append((idx, int(a), int(b)) if a.isdigit() and b.isdigit()
                                else (idx, None, None))
                yield fields, vals
        finally:
            f.close()
    return records(), {UNITS[i]: stage for i, _, _, stage in controls}


def write(path, rows, columns):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t',
                           extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def digest(path):
    h = hashlib.sha256()
    with open(str(path), 'rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def fasta_index(path):
    ans = {}
    with open(str(path)) as f:
        for line in f:
            v = line.rstrip('\n').split('\t')
            if len(v) < 5 or v[0] in ans:
                raise ValueError('Invalid FASTA .fai')
            ans[v[0]] = tuple(int(x) for x in v[1:5])
    return ans


def fasta_base(f, fai, chrom, pos):
    if chrom not in fai:
        return None
    length, offset, width, bytesperline = fai[chrom]
    if not 1 <= pos <= length or width <= 0:
        return None
    p = pos - 1
    f.seek(offset+(p//width)*bytesperline+(p%width))
    return f.read(1).decode('ascii').upper()


def alias_map(path):
    if path is None:
        return {}
    out = {}
    with open(str(path), newline='') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            source, dest = r.get('dryad_scaffold'), r.get('fasta_contig')
            if not source or not dest or source in out:
                raise ValueError('Invalid duplicate/empty contig alias')
            out[source] = dest
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('datas1', 'datas2', 'datas3', 'out-dir'):
        p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--min-depth', type=int, default=6)
    p.add_argument('--genome', type=Path)
    p.add_argument('--contig-map', type=Path)
    p.add_argument('--min-genome-match', type=float, default=0.99)
    a = p.parse_args()
    if a.min_depth < 1 or not 0 < a.min_genome_match <= 1:
        p.error('Invalid coverage/reference-match threshold')
    if a.contig_map and not a.genome:
        p.error('Contig aliases require --genome')
    dest = a.out_dir
    candidates = dest/'lineA_lineB.published_like_candidates.tsv'
    markers = dest/'lineA_lineB.published_like.markers.tsv'
    if candidates.exists() or markers.exists():
        raise ValueError('Output exists; choose a new versioned --out-dir')

    stats = collections.Counter()
    records = {}
    stages = {}
    coords = set()
    sources = [('DataS1', a.datas1, 'pla'), ('DataS2', a.datas2, 'emb'),
               ('DataS3', a.datas3, 'fol')]
    for table_name, path, tissue in sources:
        it, current = table(path, tissue)
        for unit, stage in current.items():
            if unit in stages and stages[unit] != stage:
                raise ValueError('Cross-tissue biological ID/stage mismatch: '+unit)
            stages[unit] = stage
        matched = set()
        for r, values in it:
            sid = r['SNP_ID'].strip()
            key = (r['Scaffold'].strip(), r['Position'].strip(),
                   r['Ref_allele'].strip().upper(), r['Alt_allele'].strip().upper())
            if table_name == 'DataS1':
                stats['placental_input_rows'] += 1
                if sid in records:
                    raise ValueError('Repeated DataS1 SNP_ID')
                if key[1].isdigit() and int(key[1]) > 0:
                    coord = (key[0], int(key[1]))
                    if coord in coords:
                        raise ValueError('Repeated genomic coordinate')
                    coords.add(coord)
                records[sid] = [key, r.get('Gene_ID', ''), r['Annotation'],
                                [[0, 0] for _ in UNITS], False]
            else:
                if sid not in records:
                    continue
                if sid in matched:
                    raise ValueError('Repeated SNP_ID in '+table_name)
                matched.add(sid)
                stats[table_name+'_matched_rows'] += 1
                if key != records[sid][0]:
                    records[sid][4] = True
                    stats['cross_tissue_key_mismatch'] += 1
                    continue
            state = records[sid]
            for index, ref, alt in values:
                if ref is None or alt is None:
                    state[4] = True
                    stats['invalid_parental_allele_counts'] += 1
                else:
                    state[3][index][0] += ref
                    state[3][index][1] += alt
    if set(stages) != set(UNITS):
        raise ValueError('Incomplete matched within-line control identities')

    selected = []
    tiers = collections.Counter()
    for sid, (key, gene, effect, counts, invalid) in records.items():
        if invalid:
            stats['excluded_invalid_count_or_key_mismatch'] += 1
            continue
        chrom, position, ref, alt = key
        if (not chrom or not position.isdigit() or int(position) < 1 or
                len(ref) != 1 or len(alt) != 1 or ref == alt or
                ref not in 'ACGT' or alt not in 'ACGT'):
            stats['excluded_invalid_position_or_snp'] += 1
            continue
        a_support = []
        b_support = []
        contradict = False
        for index, (nr, na) in enumerate(counts):
            dp = nr+na
            if dp < a.min_depth:
                continue
            if index < 3:
                if na:
                    contradict = True
                else:
                    a_support.append(dp)
            else:
                if nr:
                    contradict = True
                else:
                    b_support.append(dp)
        if contradict:
            stats['excluded_contradictory_callable_control'] += 1
            continue
        if not a_support or not b_support:
            stats['excluded_insufficient_line_support'] += 1
            continue
        tier = ('replicated_both_lines' if len(a_support) >= 2 and len(b_support) >= 2
                else 'replicated_one_line' if len(a_support) >= 2 or len(b_support) >= 2
                else 'one_individual_each_line')
        tiers[tier] += 1
        selected.append(dict(chrom=chrom, genomic_pos=int(position), ref=ref, alt=alt,
                             focal_index=0, alternate_index=1, snp_id=sid,
                             gene_id=gene, annotation=effect,
                             n_line_a_controls=len(a_support),
                             n_line_b_controls=len(b_support),
                             min_control_depth=min(a_support+b_support),
                             evidence_tier=tier))
    selected.sort(key=lambda r: (r['chrom'], r['genomic_pos']))
    dest.mkdir(parents=True, exist_ok=True)
    write(candidates, selected, FIELDS)
    report = dict(stats)
    report.update(selected_total=len(selected), study_benchmark=55224,
                  difference_from_benchmark=len(selected)-55224,
                  overlap_with_published_panel_not_known=True,
                  evidence_tiers=dict(tiers), min_depth_per_pooled_fish=a.min_depth,
                  parental_orientation='Line_A_REF_0_Line_B_ALT_1',
                  input_sha256={name: digest(path) for name, path, _ in sources},
                  genome_validation='NOT_PERFORMED', production_marker_file=None)
    if a.genome:
        fai_path = Path(str(a.genome)+'.fai')
        if not a.genome.is_file() or not fai_path.is_file():
            raise ValueError('Missing assembly-matched genome FASTA or index (.fai)')
        fai = fasta_index(fai_path)
        aliases = alias_map(a.contig_map)
        checked = []
        accepted = []
        with a.genome.open('rb') as f:
            for r in selected:
                contig = aliases.get(r['chrom'], r['chrom'])
                nucleotide = fasta_base(f, fai, contig, r['genomic_pos'])
                status = ('not_mappable' if nucleotide is None else
                          'ref_match' if nucleotide == r['ref'] else 'ref_mismatch')
                checked.append(dict(snp_id=r['snp_id'], dryad_scaffold=r['chrom'],
                                    fasta_contig=contig, genomic_pos=r['genomic_pos'],
                                    dryad_ref=r['ref'], fasta_base=nucleotide or '',
                                    status=status))
                if status == 'ref_match':
                    x = r.copy()
                    x['chrom'] = contig
                    accepted.append(x)
        write(dest/'lineA_lineB.genome_audit.tsv', checked,
              ('snp_id', 'dryad_scaffold', 'fasta_contig', 'genomic_pos',
               'dryad_ref', 'fasta_base', 'status'))
        audit = collections.Counter(z['status'] for z in checked)
        ratio = len(accepted)/len(selected) if selected else 0
        report['genome_match_rate'] = round(ratio, 6)
        report['genome_audit_counts'] = dict(audit)
        if ratio >= a.min_genome_match:
            if len({(x['chrom'], x['genomic_pos']) for x in accepted}) != len(accepted):
                raise ValueError('Contig aliases collapse different sites')
            write(markers, accepted, BASE_FIELDS)
            report['production_marker_file'] = str(markers)
            report['genome_validation'] = 'REFERENCE_BASE_MATCH_RNA_PROXY'
        else:
            report['genome_validation'] = 'FAILED_REFERENCE_BASE_MATCH'
    with (dest/'lineA_lineB.published_like_report.json').open('w') as out:
        json.dump(report, out, indent=2, sort_keys=True)
        out.write('\n')
    print(json.dumps(report, indent=2, sort_keys=True))
    if a.genome and not report['production_marker_file']:
        raise SystemExit('No production markers: genome validation failed')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as error:
        raise SystemExit('Fish marker preparation failed: '+str(error))
