#!/usr/bin/env python3
"""Conservative RNA-derived Line A/B markers from Zhang et al. Dryad DataS1.

This deliberately excludes F1 expression from allele-assignment decisions.
Output is CANDIDATES by default, NOT validated germline genotypes. Production
markers require successful genomic FASTA/FAI reference verification.
Python >= 3.6, standard library only. Never modifies the Dryad input.
"""
import argparse
import collections
import csv
import hashlib
import json
import re
from pathlib import Path

CORE = ('SNP_ID', 'Scaffold', 'Position', 'Ref_allele', 'Alt_allele', 'Annotation')
SAMPLE = re.compile(r'^poePro([AB][0-9]+)x([AB][0-9]+)_(pla|emb)([0-9]+)_stage([0-9]+)$')
OUTPUT_COLS = ('chrom', 'genomic_pos', 'ref', 'alt', 'focal_index',
               'alternate_index', 'snp_id', 'gene_id', 'annotation',
               'n_line_a_controls', 'n_line_b_controls', 'min_control_depth',
               'evidence')
CODING = ('synonymous_variant', 'missense_variant', 'stop_gained',
          'stop_lost', 'stop_retained_variant', 'start_lost',
          'start_retained_variant', 'coding_sequence_variant')


def table_reader(path):
    handle = open(str(path), 'r', newline='', encoding='utf-8-sig')
    reader = csv.reader(handle)
    header = next(reader)
    if header and header[0].startswith('Data S') and not any(header[1:]):
        header = next(reader)
    if any(key not in header for key in CORE):
        handle.close()
        raise ValueError('Missing DataS1 columns')
    if len(header) != len(set(header)):
        handle.close()
        raise ValueError('Duplicate DataS1 header names')
    def records():
        try:
            for line_no, values in enumerate(reader, 3):
                if not values:
                    continue
                if len(values) != len(header):
                    raise ValueError('Bad row length near line {}'.format(line_no))
                yield dict(zip(header, values))
        finally:
            handle.close()
    return header, records()


def sample_groups(header, tissue="pla"):
    refs = {col[:-len('_ref_count')] for col in header if col.endswith('_ref_count')}
    alts = {col[:-len('_alt_count')] for col in header if col.endswith('_alt_count')}
    if refs != alts:
        raise ValueError('Unpaired allele-count sample columns')
    groups = {'A': [], 'B': [], 'AB': [], 'BA': []}
    for name in sorted(refs):
        match = SAMPLE.fullmatch(name)
        if not match:
            raise ValueError('Unexpected placental sample label: {}'.format(name))
        m, p, obs_tissue = match.groups()[:3]
        if obs_tissue != tissue:
            raise ValueError('Unexpected tissue in sample label: {}'.format(name))
        group = m[0] if m[0] == p[0] else m[0] + p[0]
        groups[group].append(name)
    if tuple(len(groups[k]) for k in ('A', 'B', 'AB', 'BA')) != (3, 2, 3, 3):
        raise ValueError('Expected 3 A controls, 2 B controls, 3 AB and 3 BA hybrids')
    return groups


def counts(row, sample):
    values = (row[sample+'_ref_count'].strip(), row[sample+'_alt_count'].strip())
    if any(not re.fullmatch(r'[0-9]+', v) for v in values):
        raise ValueError('Noninteger/negative allele count at {} {}'.format(row['SNP_ID'], sample))
    return tuple(int(v) for v in values)


def consensus(row, samples, min_depth, minor_fraction):
    """(allele index, count, minimum called depth), None if conflicting/heterozygous.

    Missing controls (<min_depth) are tolerated; every *callable* control must
    agree on the same near-homozygous direction. Min number checked separately.
    """
    orientations, depths = [], []
    for sample in samples:
        ref, alt = counts(row, sample)
        depth = ref + alt
        if depth < min_depth:
            continue
        fraction = alt / depth
        if fraction <= minor_fraction:
            orientations.append(0)
        elif fraction >= 1 - minor_fraction:
            orientations.append(1)
        else:
            return None
        depths.append(depth)
    if not orientations or len(set(orientations)) != 1:
        return None
    return orientations[0], len(orientations), min(depths)



def embryo_corroboration(path, selected, coding, depth=20, minor_fraction=.10):
    """Independent-tissue allelic check, but often the *same embryo*, not an independent unit.

    Only concordant near-homozygous parental-control embryo RNA is counted;
    no hybrid embryo expression is used to assign parental alleles.
    """
    hdr, reader = table_reader(path)
    groups = sample_groups(hdr, tissue='emb')
    by_key = {(x['chrom'], x['genomic_pos'], x['ref'], x['alt']): x
              for x in selected}
    stats = collections.Counter()
    for row in reader:
        pos = row['Position'].strip()
        if not pos.isdigit():
            continue
        key = (row['Scaffold'].strip(), int(pos), row['Ref_allele'].strip().upper(),
               row['Alt_allele'].strip().upper())
        candidate = by_key.get(key)
        if candidate is None:
            continue
        if 'embryo_status' in candidate:
            raise ValueError('Repeated embryo row for candidate: '+str(key))
        try:
            def support(samples, target):
                n, bad = 0, 0
                for sample in samples:
                    ref, alt = counts(row, sample)
                    dp = ref + alt
                    if dp < depth:
                        continue
                    f = alt/dp
                    if (target == 0 and f <= minor_fraction) or (target == 1 and f >= 1-minor_fraction):
                        n += 1
                    else:
                        bad += 1
                return n, bad
            na, bada = support(groups['A'], candidate['focal_index'])
            nb, badb = support(groups['B'], candidate['alternate_index'])
        except ValueError:
            na, nb, bada, badb = 0, 0, 1, 0
            stats['invalid_embryo_controls'] += 1
        candidate['embryo_n_line_a_controls'] = na
        candidate['embryo_n_line_b_controls'] = nb
        candidate['embryo_n_conflicts'] = bada + badb
        candidate['embryo_status'] = ('corroborated' if na >= 1 and nb >= 1 and not bada+badb
                                      else 'conflicting' if bada+badb else 'insufficient_depth')
    for candidate in selected:
        if 'embryo_status' not in candidate:
            candidate.update(embryo_status='not_in_DataS2', embryo_n_line_a_controls=0,
                             embryo_n_line_b_controls=0, embryo_n_conflicts=0)
        stats[candidate['embryo_status']] += 1
    supported = [x for x in selected if x['embryo_status'] == 'corroborated']
    supported_coding = [x for x in coding if x['embryo_status'] == 'corroborated']
    return supported, supported_coding, dict(stats)

def read_fai(path):
    fai = {}
    with open(str(path), 'r') as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 5:
                raise ValueError('Invalid FASTA index line')
            name, length, offset, width, bytes_per_line = fields[:5]
            if name in fai:
                raise ValueError('Duplicate FASTA contig: '+name)
            fai[name] = tuple(map(int, (length, offset, width, bytes_per_line)))
    return fai


def fasta_base(handle, fai, chrom, pos):
    if chrom not in fai:
        return None
    length, offset, width, bytes_per_line = fai[chrom]
    if pos < 1 or pos > length or not width:
        return None
    zero = pos - 1
    handle.seek(offset + (zero // width)*bytes_per_line + zero % width)
    return handle.read(1).decode('ascii').upper()


def mapping(path):
    if not path:
        return {}
    result = {}
    with open(str(path), 'r', newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if not row.get('dryad_scaffold') or not row.get('fasta_contig'):
                raise ValueError('Map requires dryad_scaffold and fasta_contig columns')
            if row['dryad_scaffold'] in result:
                raise ValueError('Duplicate contig alias')
            result[row['dryad_scaffold']] = row['fasta_contig']
    return result


def sha256(path):
    digest = hashlib.sha256()
    with open(str(path), 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields, lineterminator='\n',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--datas1', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--datas2', type=Path, help='Embryo RNA table, for cross-tissue corroboration')
    parser.add_argument('--min-depth', type=int, default=20)
    parser.add_argument('--max-minor-fraction', type=float, default=0.10)
    parser.add_argument('--min-a-controls', type=int, default=2)
    parser.add_argument('--min-b-controls', type=int, default=2)
    parser.add_argument('--genome', type=Path, help='Genome FASTA indexed with samtools faidx')
    parser.add_argument('--contig-map', type=Path, help='TSV: dryad_scaffold fasta_contig; verified aliases only')
    parser.add_argument('--min-genome-match', type=float, default=.99)
    args = parser.parse_args()
    if args.min_depth < 1 or not 0 < args.max_minor_fraction < 0.5 or not 0 < args.min_genome_match <= 1:
        parser.error('Invalid threshold')
    if not 1 <= args.min_a_controls <= 3 or not 1 <= args.min_b_controls <= 2:
        parser.error('Invalid number of independent parental controls')
    if args.contig_map and not args.genome:
        parser.error('--contig-map requires --genome')
    if args.genome and not args.datas2:
        parser.error('Production marker promotion requires --datas2 embryo corroboration')
    # Never leave a previously verified production marker file in place when
    # a later invocation uses different inputs or fails validation.
    if (args.out_dir/'lineA_lineB.fixed_exonic.markers.tsv').exists():
        raise ValueError('Production markers already exist in out-dir; use a new versioned directory')
    header, iterator = table_reader(args.datas1)
    groups = sample_groups(header)
    stats = collections.Counter()
    selected, coding = [], []
    seen = set()
    for row in iterator:
        stats['input_rows'] += 1
        chrom, pos = row['Scaffold'].strip(), row['Position'].strip()
        ref, alt = row['Ref_allele'].strip().upper(), row['Alt_allele'].strip().upper()
        if not pos.isdigit() or int(pos) < 1:
            stats['invalid_or_imprecise_positions'] += 1
            continue
        if len(ref) != 1 or len(alt) != 1 or ref == alt or ref not in 'ACGT' or alt not in 'ACGT':
            stats['invalid_snp_alleles'] += 1
            continue
        key = (chrom, int(pos))
        if key in seen:
            stats['duplicate_position_rows'] += 1
            continue  # Fail closed on duplicates: any duplicates cause whole-run stop below.
        seen.add(key)
        try:
            a = consensus(row, groups['A'], args.min_depth, args.max_minor_fraction)
            b = consensus(row, groups['B'], args.min_depth, args.max_minor_fraction)
        except ValueError:
            # The deposited table has a handful of negative read counts.
            # Such records are invalid as allele counts and cannot qualify.
            stats['invalid_parental_counts'] += 1
            continue
        if not a or a[1] < args.min_a_controls or not b or b[1] < args.min_b_controls:
            stats['control_support_insufficient_or_conflicting'] += 1
            continue
        if a[0] == b[0]:
            stats['same_allele_in_both_lines'] += 1
            continue
        candidate = dict(chrom=chrom, genomic_pos=int(pos), ref=ref, alt=alt,
                         focal_index=a[0], alternate_index=b[0],
                         snp_id=row['SNP_ID'], gene_id=row.get('Gene_ID', ''),
                         annotation=row['Annotation'], n_line_a_controls=a[1],
                         n_line_b_controls=b[1], min_control_depth=min(a[2], b[2]),
                         evidence='RNAseq_near_homozygous_controls_unvalidated_genomic')
        selected.append(candidate)
        if any(term in row['Annotation'].split('&') for term in CODING):
            coding.append(candidate)
    if stats['duplicate_position_rows']:
        raise ValueError('Repeated positions in Dryad table; resolve before producing marker candidates')
    selected.sort(key=lambda x: (x['chrom'], x['genomic_pos']))
    coding.sort(key=lambda x: (x['chrom'], x['genomic_pos']))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.out_dir/'lineA_lineB.RNA_controls_candidates.tsv', selected, OUTPUT_COLS)
    write_tsv(args.out_dir/'lineA_lineB.RNA_controls_coding_candidates.tsv', coding, OUTPUT_COLS)
    production = coding
    if args.datas2:
        embryo_all, embryo_coding, embryo_stats = embryo_corroboration(
            args.datas2, selected, coding, args.min_depth, args.max_minor_fraction)
        extra = OUTPUT_COLS + ('embryo_status', 'embryo_n_line_a_controls',
                               'embryo_n_line_b_controls', 'embryo_n_conflicts')
        write_tsv(args.out_dir/'lineA_lineB.RNA_controls_embryo_corroborated.tsv', embryo_all, extra)
        write_tsv(args.out_dir/'lineA_lineB.RNA_controls_coding_embryo_corroborated.tsv', embryo_coding, extra)
        production = embryo_coding
    else:
        embryo_stats = None
    stats = dict(stats)
    stats.update({'all_candidate_sites': len(selected), 'coding_candidate_sites': len(coding),
                  'a_control_samples': groups['A'], 'b_control_samples': groups['B'],
                  'excluded_from_parent_inference': groups['AB']+groups['BA'],
                  'genomic_validation': 'NOT_PERFORMED', 'production_marker_file': None,
                  'embryo_corroboration_counts': embryo_stats,
                  'coding_embryo_corroborated_sites': len(production) if args.datas2 else None,
                  'min_depth_per_control': args.min_depth,
                  'max_minor_allele_fraction': args.max_minor_fraction,
                  'minimum_line_a_controls': args.min_a_controls,
                  'minimum_line_b_controls': args.min_b_controls,
                  'minimum_genome_match': args.min_genome_match,
                  'DataS1_sha256': sha256(args.datas1),
                  'DataS2_sha256': sha256(args.datas2) if args.datas2 else None})
    if args.genome:
        if not args.genome.is_file() or not Path(str(args.genome)+'.fai').is_file():
            raise ValueError('Genome FASTA and FASTA.fai required (run samtools faidx)')
        fai = read_fai(Path(str(args.genome)+'.fai'))
        aliases = mapping(args.contig_map)
        audits, validated = [], []
        with open(str(args.genome), 'rb') as fasta:
            for row in production:
                name = aliases.get(row['chrom'], row['chrom'])
                base = fasta_base(fasta, fai, name, row['genomic_pos'])
                status = ('unmapped_or_outside_fasta' if base is None else
                          'ref_match' if base == row['ref'] else 'ref_mismatch')
                audits.append({'dryad_scaffold': row['chrom'], 'fasta_contig': name,
                               'position': row['genomic_pos'], 'dryad_ref': row['ref'],
                               'fasta_base': base or '', 'status': status})
                if status == 'ref_match':
                    valid = row.copy()
                    valid['chrom'] = name
                    validated.append(valid)
        write_tsv(args.out_dir/'lineA_lineB.genome_validation.tsv', audits,
                  ('dryad_scaffold', 'fasta_contig', 'position', 'dryad_ref', 'fasta_base', 'status'))
        counts = collections.Counter(a['status'] for a in audits)
        stats['genome_audit_counts'] = dict(counts)
        total = len(audits)
        match_rate = counts['ref_match']/total if total else 0.0
        stats['genome_match_over_all_coding_candidates'] = round(match_rate, 6)
        if match_rate < args.min_genome_match:
            stats['genomic_validation'] = 'FAILED_TOO_FEW_MAPPABLE_MATCHING_REFERENCE_BASES'
            print('NO PRODUCTION MARKERS: genome match {:.2%} below {:.2%}'.format(
                match_rate, args.min_genome_match))
        else:
            keys = [(x['chrom'], x['genomic_pos']) for x in validated]
            if len(set(keys)) != len(keys):
                raise ValueError('Contig aliases introduce duplicate reference coordinates')
            path = args.out_dir/'lineA_lineB.fixed_exonic.markers.tsv'
            write_tsv(path, validated, OUTPUT_COLS[:6])
            stats['genomic_validation'] = 'PASSED_REFERENCE_BASE_CHECK_RNA_DERIVED_PROXY'
            stats['production_marker_file'] = str(path)
    path = args.out_dir/'lineA_lineB.marker_report.json'
    with path.open('w') as output:
        json.dump(stats, output, indent=2, sort_keys=True)
    print(json.dumps(stats, indent=2, sort_keys=True))
    if args.genome and not stats['production_marker_file']:
        raise SystemExit('Marker promotion blocked by genome mismatch; inspect genome_validation.tsv')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as exc:
        raise SystemExit('Fish marker preparation failed: {}'.format(exc))
