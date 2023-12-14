# Import genes from assemblies using co-ordinates defined in https://github.com/oscarlr/VDJbase-Genomics

# The input file is a 'gene file' produced by make_gene_file.py, which uses co-ordinates provided by the
# bed files to determine the sequence of each element of the gene.

# This code:
# Identifies alleles, producing 'snp' names for those not in the reference set
# Provides note on functionality of the allele
# splits out subject name into project/subject
# Sequences listed in the output are in +ve (5' to 3' sense). The sense parameter specifies the sense in the input file. 
# if the sense parameter is specified as +-, sense is taken from the 'sense' column  and both + and - can be present.
# Cigar strings listed in the output are identical to those in the input file and therefore reflect the alignment of the sequence to the reference.
import argparse
import csv
import os
from hashlib import sha256
from operator import itemgetter

from receptor_utils import simple_bio_seq as simple
from receptor_utils.novel_allele_name import name_novel

import read_bed

class GeneParsingException(Exception):
    pass


# find the closest reference allele, using the allele name and gene name provided in the row
def closest_ref_allele(row, refs, coding_sequence):
    chain = row['gene'][3]
    if chain == 'V':
        seq_choices = [k for k in refs['V'].keys() if row['gene'] in k and refs['V-gapped'][k][0] != '.']
    else:
        seq_choices = [k for k in refs[chain].keys() if row['gene'] in k]

    if not seq_choices:
        return None

    scores = [(simple.nt_diff(coding_sequence, refs[chain][seq_name]), seq_name) for seq_name in seq_choices]
    scores.sort(key=itemgetter(0))

    return scores[0][1], refs[chain][scores[0][1]], scores[0][0]


def best_ref_match(closest_imgt_seq, seq):
    # find the best matching reference allele
    diffs = []
    window = 10
    best_pos, best_diff = 0, 999
    max_pos = len(seq) - len(closest_imgt_seq)

    for i in range(max_pos - window, max_pos + window):
        if i < 0:
            continue
        if i + len(closest_imgt_seq) <= len(seq):
            #s1 = seq[i:i + len(closest_imgt_seq)]
            #s2 = closest_imgt_seq
            diff = simple.nt_diff(seq[i:i + len(closest_imgt_seq)], closest_imgt_seq)
            #print(f"1: {i} {len(seq[i:i + len(closest_imgt_seq)])}  {len(closest_imgt_seq)}")
        else:
            #s1 = seq[i:]
            #s2 = closest_imgt_seq[:len(seq[:len(seq[i:])])]
            diff = simple.nt_diff(seq[i:], closest_imgt_seq[:len(seq[:len(seq[i:])])])
            #print(f"2: {i} {len(seq[i:])}  {len(closest_imgt_seq[:len(seq[:len(seq[i:])])])}")

        diffs.append(diff)
        if diff < best_diff:
            best_pos, best_diff = i, diff

    coding_sequence = seq[best_pos:]
    score = len(coding_sequence) - best_diff
    imgt_allele_match = round(100 * score / len(coding_sequence), 2) if coding_sequence else 0

    return coding_sequence, score, imgt_allele_match


def read_refs(locus, ref_dir):
    refs = {}
    refs_inverted = {}

    ref_types = ['V', 'J', 'C']

    if locus in ['IGH', 'TRB', 'TRD']:
        ref_types.append('D')

    for ref_type in ref_types:
        refs[ref_type] = simple.read_fasta(os.path.join(ref_dir, 'Homo_sapiens_%s%s.fasta' % (locus, ref_type)))
        refs_inverted[ref_type] = dict([(v, k) for k, v in refs[ref_type].items()])
    refs['V-gapped'] = simple.read_fasta(os.path.join(ref_dir, 'Homo_sapiens_%sV_gapped.fasta' % locus))

    return refs, refs_inverted


# find a full-length reference allele to use as a template for positional searching
def find_template_allele(row, refs):
    chain = row['gene'][3]
    if chain == 'V':
        seq_choices = [k for k in refs['V'].keys() if row['gene'] in k and refs['V-gapped'][k][0] != '.']
    else:
        seq_choices = [k for k in refs[chain].keys() if row['gene'] in k]

    if seq_choices:
        return seq_choices[0], refs[chain][seq_choices[0]]
    else:
        return None

d_fields = [
    ('5_nonamer', 'D-5_NONAMER'),
    ('5_spacer', 'D-5_SPACER'),
    ('5_heptamer', 'D-5_HEPTAMER'),
    ('exon_1', 'D-REGION'),
    ('3_heptamer', 'D-3_HEPTAMER'),
    ('3_spacer', 'D-3_SPACER'),
    ('3_nonamer', 'D-3_NONAMER'),
]


def process_d_gene(refname, sense, beds, notes, refs, res, row):
    gene = res['genotyper_gene']
    coords = beds[refname][gene]

    for bed_el, row_el in d_fields:
        res[row_el] = simple.reverse_complement(row[row_el]) if sense == '-' else row[row_el]
        res[row_el + '_CIGAR'] = row[row_el + '_CIGAR']
    #    res[row_el + '_CIGAR'] = Cigar(row[row_el + '_CIGAR'])._reverse_cigar() if sense == '-' else row[row_el + '_CIGAR']

    if row['D-5_HEPTAMER'] and row['D-5_HEPTAMER'] != coords['5_HEPTAMER']['seq']:
        notes.append('5 Heptamer does not match reference')

    if row['D-5_NONAMER'] and row['D-5_NONAMER'] != coords['5_NONAMER']['seq']:
        notes.append('5 Nonamer does not match reference')

    if row['D-3_HEPTAMER'] and row['D-3_HEPTAMER'] != coords['3_HEPTAMER']['seq']:
        notes.append('3 Heptamer does not match reference')

    if row['D-3_NONAMER'] and row['D-3_NONAMER'] != coords['3_NONAMER']['seq']:
        notes.append('3 Nonamer does not match reference')

    try:
        ret = closest_ref_allele(row, refs, res['D-REGION'])

        if not ret:
            print(f'no reference allele found for {row["gene"]}')
            return False
        
        (closest_imgt_allele, closest_imgt_seq, score) = ret
        res['closest_imgt_allele'] = closest_imgt_allele
        res['closest_imgt_seq'] = closest_imgt_seq

        # assign a name based on closest sequence in the reference set of alleles for this gene
        # If we had an exact match, use just that reference allele - this gives us a good name for extensions

        if len(res['D-REGION']):
            res['imgt_allele_match'] = 100.0 - round(100 * score / len(res['D-REGION']), 2)
            if closest_imgt_seq in res['D-REGION'] and closest_imgt_allele.split('*')[0] == row['gene']:
                gene_refs = {res['closest_imgt_allele']: refs['D'][closest_imgt_allele]}
            else:
                gene_refs = {name: seq for name, seq in refs['D'].items() if name.split('*')[0] == row['gene']}
            res['vdjbase_allele'], _, strategy = name_novel(res['D-REGION'], gene_refs, v_gene=False)
            notes.append(strategy)

            if len(res['vdjbase_allele']) > 80:
                res['vdjbase_allele'] = \
                    res['vdjbase_allele'].split('*')[0] + '*' \
                    + res['vdjbase_allele'].split('*')[1].split('_')[0] \
                    + '_' + sha256(res['D-REGION'].encode('utf-8')).hexdigest()[-4:]

        res['notes'] = '\\n'.join([x for x in notes if x])

    except GeneParsingException as e:
        print(e)
        return False

    return True


j_fields = [
    ('nonamer', 'J-NONAMER'),
    ('spacer', 'J-SPACER'),
    ('heptamer', 'J-HEPTAMER'),
    ('exon_1', 'J-REGION'),
]

def process_j_gene(refname, sense, beds, notes, refs, res, row):
    gene = res['genotyper_gene']
    coords = beds[refname][gene]

    for bed_el, row_el in j_fields:
        res[row_el] = simple.reverse_complement(row[row_el]) if sense == '-' else row[row_el]
        res[row_el + '_CIGAR'] = row[row_el + '_CIGAR']
    #    res[row_el + '_CIGAR'] = Cigar(row[row_el + '_CIGAR'])._reverse_cigar()  if sense == '-' else row[row_el + '_CIGAR']

    if row['J-HEPTAMER'] and row['J-HEPTAMER'] != coords['HEPTAMER']['seq']:
        notes.append('Heptamer does not match reference')

    if row['J-NONAMER'] and row['J-NONAMER'] != coords['NONAMER']['seq']:
        notes.append('Nonamer does not match reference')

    try:
        ret = closest_ref_allele(row, refs, res['J-REGION'])

        if not ret:
            print(f'no reference allele found for {row["gene"]}')
            return False

        (closest_imgt_allele, closest_imgt_seq, score) = ret
        res['closest_imgt_allele'] = closest_imgt_allele
        res['closest_imgt_seq'] = closest_imgt_seq

        # assign a name based on closest sequence in the reference set of alleles for this gene
        # If we had an exact match, use just that reference allele - this gives us a good name for extensions

        if len(res['J-REGION']):
            res['imgt_allele_match'] = 100.0 - round(100 * score / len(res['J-REGION']), 2)

            if closest_imgt_seq in res['J-REGION'] and closest_imgt_allele.split('*')[0] == row['gene']:
                gene_refs = {res['closest_imgt_allele']: refs['J'][closest_imgt_allele]}
            else:
                gene_refs = {name: seq for name, seq in refs['J'].items() if name.split('*')[0] == row['gene']}

            res['vdjbase_allele'], _, strategy = name_novel(res['J-REGION'], gene_refs, v_gene=False)
            notes.append(strategy)

            if len(res['vdjbase_allele']) > 80:
                res['vdjbase_allele'] = \
                    res['vdjbase_allele'].split('*')[0] + '*' \
                    + res['vdjbase_allele'].split('*')[1].split('_')[0] \
                    + '_' + sha256(res['J-REGION'].encode('utf-8')).hexdigest()[-4:]

        res['notes'] = '\\n'.join([x for x in notes if x])

    except GeneParsingException as e:
        print(e)
        return False

    return True


c_fields = [
    ('REGION', 'C-REGION'),
]

def process_c_gene(refname, sense, beds, notes, refs, res, row):
    for bed_el, row_el in c_fields:
        res[row_el] = simple.reverse_complement(row[row_el]) if sense == '-' else row[row_el]
        res[row_el + '_CIGAR'] = row[row_el + '_CIGAR']

    if not res['C-REGION']:
        return False

    try:
        ret = closest_ref_allele(row, refs, res['C-REGION'])

        if not ret:
            print(f'no reference allele found for {row["gene"]}')
            return False

        (closest_imgt_allele, closest_imgt_seq, score) = ret
        res['closest_imgt_allele'] = closest_imgt_allele
        res['closest_imgt_seq'] = closest_imgt_seq

        # assign a name based on closest sequence in the reference set of alleles for this gene
        # If we had an exact match, use just that reference allele - this gives us a good name for extensions

        res['imgt_allele_match'] = 100.0 - round(100 * score / len(res['C-REGION']), 2)

        if closest_imgt_seq in res['C-REGION'] and closest_imgt_allele.split('*')[0] == row['gene']:
            gene_refs = {res['closest_imgt_allele']: refs['C'][closest_imgt_allele]}
        else:
            gene_refs = {name: seq for name, seq in refs['C'].items() if name.split('*')[0] == row['gene']}
        res['vdjbase_allele'], _, strategy = name_novel(res['C-REGION'], gene_refs, v_gene=False)
        notes.append(strategy)

        if len(res['vdjbase_allele']) > 80:
            res['vdjbase_allele'] = \
                res['vdjbase_allele'].split('*')[0] + '*' \
                + res['vdjbase_allele'].split('*')[1].split('_')[0] \
                + '_' + sha256(res['C-REGION'].encode('utf-8')).hexdigest()[-4:]

        res['notes'] = '\\n'.join([x for x in notes if x])

    except GeneParsingException as e:
        print(e)
        return False

    return True


v_fields = [
    ('NONAMER', 'V-NONAMER'),
    ('SPACER', 'V-SPACER'),
    ('HEPTAMER', 'V-HEPTAMER'),
    ('INTRON', 'V-INTRON'),
    ('EXON_1', 'L-PART1'),
    ('L-PART2', 'L-PART2'),
    ('V-REGION', 'V-REGION'),
    ('UTR', 'V-UTR')
]

# Report too long or too short elements - indicative of indels
def check_v_element_lengths(row, refname, beds):
    notes = []
    for bed_el, row_el in v_fields:
        if row[row_el]:
            coord_length = beds[refname][row['gene']][bed_el]['end'] - beds[refname][row['gene']][bed_el]['start']
            element_length = len(row[row_el])

            if coord_length != element_length:
                notes.append(f"length mismatch in {row_el}: seq {element_length} vs co-ords {coord_length} - likely indel")

    return notes


def process_v_gene(refname, sense, beds, notes, refs, res, row):
    gene = res['genotyper_gene']
    coords = beds[refname][gene]
    notes.extend(check_v_element_lengths(row, refname, beds))

    for bed_el, row_el in v_fields:
        res[row_el] = simple.reverse_complement(row[row_el]) if sense == '-' else row[row_el]
        res[row_el + '_CIGAR'] = row[row_el + '_CIGAR']
    #    res[row_el + '_CIGAR'] = Cigar(row[row_el + '_CIGAR'])._reverse_cigar() if sense == '-' else row[row_el + '_CIGAR']

    if row['V-HEPTAMER'] and row['V-HEPTAMER'] != coords['HEPTAMER']['seq']:
        notes.append('Heptamer does not match reference')

    if row['V-NONAMER'] and row['V-NONAMER'] != coords['NONAMER']['seq']:
        notes.append('Nonamer does not match reference')

    try:
        ret = closest_ref_allele(row, refs, res['V-REGION'])

        if not ret:
            print(f'no reference allele found for {row["gene"]}')
            return False
        
        (closest_imgt_allele, closest_imgt_seq, score) = ret
        res['closest_imgt_allele'] = closest_imgt_allele
        res['closest_imgt_seq'] = closest_imgt_seq
        res['ref_extension'] = (closest_imgt_seq in res['V-REGION']) and (len(closest_imgt_seq) != len(res['V-REGION']))

        # assign a name based on closest sequence in the reference set of alleles for this gene
        # If we had an exact match, use just that reference allele - this gives us a good name for extensions

        if len(res['V-REGION']) > 40:       # ignore exceptionally truncated sequences
            res['imgt_allele_match'] = 100.0 - round(100 * score / len(res['V-REGION']), 2)

            if closest_imgt_seq in res['V-REGION'] and closest_imgt_allele.split('*')[0] == row['gene']:
                gene_refs = {closest_imgt_allele: refs['V-gapped'][res['closest_imgt_allele']]}
            else:
                gene_refs = {name: seq for name, seq in refs['V-gapped'].items() if name.split('*')[0] == row['gene']}

            res['vdjbase_allele'], res['V-REGION-GAPPED'], strategy,  = name_novel(res['V-REGION'].replace('-', ''), gene_refs)
            notes.append(strategy)

            if len(res['vdjbase_allele']) > 80:
                res['vdjbase_allele'] = res['vdjbase_allele'].split('*')[0] + '*' + res['vdjbase_allele'].split('*')[1].split('_')[0] + '_' + sha256(
                    res['V-REGION'].encode('utf-8')).hexdigest()[-4:]

        res['notes'] = '\\n'.join([x for x in notes if x])

    except GeneParsingException as e:
        print(e)
        return False

    return True


# Process one row from the gene file.
def process_gene(beds, refname, sense, refs, row, row_count, writer):
    res = {}

    if sense == '+-':
        sense = row['sense']

    res['sample_name'] = row['sample_name']
    res['project'] = row['project']
    res['subject'] = row['subject']
    res['haplotype'] = row['haplotype']
    res['sense'] = sense

    # Subsequent processing assumes that the allele_sequences will always be in 'canonical' IMGT order, i.e. 5' to 3'
    # CIGAR reflects the sense of the sequence with respect to the reference sequence: + means in agreement, - means in reverse complement
    res['gene_sequence'] = simple.reverse_complement(row['allele_sequence']) if sense == '-' else row['allele_sequence']
    res['gene_sequence_CIGAR'] = row['allele_sequence_CIGAR']
    #res['gene_sequence_CIGAR'] = Cigar(row['allele_sequence_CIGAR'])._reverse_cigar() if sense == '-' else row['allele_sequence_CIGAR']
    res['genotyper_gene'] = row['gene']
    res['genotyper_allele'] = row['allele']
    gene_type = res['genotyper_gene'][3]
    notes = []

    ret = None

    if gene_type == 'V':
        ret = process_v_gene(refname, sense, beds, notes, refs, res, row)
    elif gene_type == 'D':
        ret = process_d_gene(refname, sense, beds, notes, refs, res, row)
    elif gene_type == 'J':
        ret = process_j_gene(refname, sense, beds, notes, refs, res, row)
    elif gene_type == 'C':
        ret = process_c_gene(refname, sense, beds, notes, refs, res, row)

    if ret:
        row_count += 1
        writer.writerow(res)

    return row_count


v_gene_fields = [
    'L-PART1',
    'V-INTRON',
    'L-PART2',
    'V-REGION',
    'V-REGION-GAPPED',
    'V-HEPTAMER',
    'V-SPACER',
    'V-NONAMER',
    'V-UTR',
    ]

d_gene_fields = [
    'D-5_NONAMER',
    'D-5_SPACER',
    'D-5_HEPTAMER',
    'D-REGION',
    'D-3_HEPTAMER',
    'D-3_SPACER',
    'D-3_NONAMER'
    ]

j_gene_fields = [
    'J-NONAMER',
    'J-SPACER',
    'J-HEPTAMER',
    'J-REGION',
]

c_gene_fields = [
    'C-REGION',
]

w_headers = [
    'project',
    'subject',
    'sample_name',
    'haplotype',
    'genotyper_gene',
    'genotyper_allele',
    'vdjbase_allele',
    'closest_imgt_allele',
    'sense',
    'imgt_allele_match',
    'closest_imgt_seq',
    'ref_extension',
    'notes',
    'exon_reads',
    'exon_min_score',
    'gene_sequence',
    'gene_sequence_CIGAR'
     ]


def main():
    parser = argparse.ArgumentParser(description='Derive annotated genes from a VDJ gene list created by make_gene_file.py')
    parser.add_argument('locus', help='locus (e.g. IGH, IGL)')
    parser.add_argument('refname', help='reference assembly name in SAM files (e.g. chr22)')
    parser.add_argument('sense', help='sense of the genes in the reference sequence: - is 5 to 3, + is 3 to 5, +- is both (eg IGK)')
    parser.add_argument('gene_file', help='the gene list (csv)')
    parser.add_argument('ref_seq', help='the assembly reference sequence (FASTA)')
    parser.add_argument('bed_dir', help='pathname to a directory holding one or more BED files containing the gene co-ordinates in the reference sequence')
    parser.add_argument('ref_dir', help='pathname to a directory containing gapped and ungapped reference directories (FASTA')
    parser.add_argument('outfile', help='output file that will contain the annotations (csv)')
    args = parser.parse_args()

    locus = args.locus
    refname = args.refname
    sense = args.sense
    beds = read_bed.read_beds(refname, sense, args.bed_dir, simple.read_single_fasta(args.ref_seq))

    with open(args.gene_file, 'r') as fi, open(args.outfile, 'w', newline='') as fo:
        row_count = 0

        for field_type in [v_gene_fields, d_gene_fields, j_gene_fields, c_gene_fields]:
            for field in field_type:
                w_headers.append(field)
                w_headers.append(field + '_CIGAR')

        reader = csv.DictReader(fi)
        writer = csv.DictWriter(fo, fieldnames=w_headers)
        writer.writeheader()

        refs, refs_inverted = read_refs(locus, args.ref_dir)

        for row in reader:
            row_count = process_gene(beds, refname, sense, refs, row, row_count, writer)


if __name__ == '__main__':
    main()
