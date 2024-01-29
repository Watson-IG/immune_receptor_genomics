import sys
import os
import csv

GENE_WINDOW_SIZE_3 = 0
GENE_WINDOW_SIZE_5 = 0
dir = "../current"
sense = "+-"
ref_seq = "merged_refseq.fasta"

def get_gene_type(label):
    gene = label.split('*')[0]
    if gene[3] == 'J' or gene[3] == 'V':
        return gene[3]
    elif gene[3] == 'D' and '-' in gene and '_' not in gene:
        return 'D'
    else:
        return 'C'

def read_beds(sense, dir, ref_seq):
    beds = {}
    for entry in os.scandir(dir):
        if entry.is_file() and '.bed' in entry.name:
            with open(os.path.join(dir, entry.name), 'r') as fi:
                el_type = entry.name.replace('.bed', '')
                reader = csv.DictReader(fi, delimiter='\t', fieldnames=['chain', 'start', 'end', 'gene', 'sense'])
                for row in reader:
                    if row['start'] and row['end']:
                        refname = row['chain']
                        row['el_type'] = el_type
                        row['start'] = int(row['start'])
                        row['end'] = int(row['end'])
                        row['seq'] = ref_seq[refname][row['start']:row['end']] if refname in ref_seq else ''
                        if refname not in beds:
                            beds[refname] = {}
                        if row['gene'] not in beds[refname]:
                            beds[refname][row['gene']] = {}
                        if row['el_type'] in beds[refname][row['gene']]:
                            beds[refname][row['gene']]['3_' + row['el_type']] = beds[refname][row['gene']][row['el_type']]
                            del beds[refname][row['gene']][row['el_type']]
                            beds[refname][row['gene']]['5_' + row['el_type']] = row
                        else:
                            beds[refname][row['gene']][row['el_type']] = row

    # Sanity checks
    if sense == '-':
        g_s = 'start'
        g_e = 'end'
    else:
        g_s = 'end'
        g_e = 'start'

    for refname in beds:
        for gene in list(beds[refname]):
            gene_type = get_gene_type(gene)
            row = beds[refname][gene]
            if gene_type == 'V':
                complete = True
                for el in ['EXON_1', 'INTRON', 'EXON_2', 'HEPTAMER', 'SPACER', 'NONAMER']:
                    if el not in row:
                        sys.exit(f"element {el} missing from gene {gene} in file {refname}")
                        complete = False

                if complete:
                    if sense == '+-':
                        if row['NONAMER']['start'] > row['HEPTAMER']['start']:
                            g_s = 'end'
                            g_e = 'start'
                        else:
                            g_s = 'start'
                            g_e = 'end'

                    if row['NONAMER'][g_e] != row['SPACER'][g_s]:
                        sys.exit(f"maths problem in {gene}: row['NONAMER'][{g_e}] != row['SPACER'][{g_s}]")
                    if row['SPACER'][g_e] != row['HEPTAMER'][g_s]:
                        sys.exit(f"maths problem in {gene}: row['SPACER'][{g_e}] != row['HEPTAMER'][{g_s}]")
                    if row['HEPTAMER'][g_e] != row['EXON_2'][g_s]:
                        sys.exit(f"maths problem in {gene}: row['HEPTAMER'][{g_e}] != row['EXON_2'][{g_s}]")
                    if row['EXON_2'][g_e] != row['INTRON'][g_s]:
                        sys.exit(f"maths problem in {gene}: row['EXON_2'][{g_e}] != row['INTRON'][{g_s}]")
                    if row['INTRON'][g_e] != row['EXON_1'][g_s]:
                        sys.exit(f"maths problem in {gene}: row['INTRON'][{g_e}] != row['EXON_1'][{g_s}]")

                if 'GENE' in row:
                    row['GENE']['start'] = row['GENE']['start'] - GENE_WINDOW_SIZE_5
                    row['GENE']['end'] = row['GENE']['end'] + GENE_WINDOW_SIZE_3

                row_sense = sense
                if sense == '+-':
                    if row['NONAMER']['start'] > row['HEPTAMER']['start']:
                        row_sense = '+'
                    else:
                        row_sense = '-'

                if row_sense == '-':
                    if 'V-REGION' not in row:
                        row['V-REGION'] = {}
                        row['V-REGION']['start'] = row['EXON_2']['start']
                        row['V-REGION']['end'] = row['EXON_2']['end'] - 11
                    if 'L-PART2' not in row:
                        row['L-PART2'] = {}
                        row['L-PART2']['start'] = row['EXON_2']['end'] - 11
                        row['L-PART2']['end'] = row['EXON_2']['end']
                else:
                    if 'V-REGION' not in row:
                        row['V-REGION'] = {}
                        row['V-REGION']['start'] = row['EXON_2']['start'] + 11
                        row['V-REGION']['end'] = row['EXON_2']['end']
                    if 'L-PART2' not in row:
                        row['L-PART2'] = {}
                        row['L-PART2']['start'] = row['EXON_2']['start']
                        row['L-PART2']['end'] = row['EXON_2']['start'] + 11

            elif gene_type == 'J':
                for el in ['HEPTAMER', 'SPACER', 'NONAMER']:
                    complete = True
                    if el not in row:
                        print(f'element {el} missing from row {row}')
                        complete = False

                if complete:
                    if sense == '+-':
                        if row['NONAMER']['start'] < row['HEPTAMER']['start']:
                            g_s = 'end'
                            g_e = 'start'
                        else:
                            g_s = 'start'
                            g_e = 'end'

                    if row['HEPTAMER'][g_e] != row['SPACER'][g_s]:
                        print(f"maths problem in {gene}: row['HEPTAMER'][{g_e}] != row['SPACER'][{g_s}]")
                    if row['SPACER'][g_e] != row['NONAMER'][g_s]:
                        print(f"maths problem in {gene}: row['SPACER'][{g_e}] != row['NONAMER'][{g_s}]")

                row['GENE'][g_e] -= GENE_WINDOW_SIZE_5
                row['GENE'][g_s] += GENE_WINDOW_SIZE_3

            elif gene_type == 'D':
                for el in ['3_HEPTAMER', '3_SPACER', '3_NONAMER', '5_HEPTAMER', '5_SPACER', '5_NONAMER']:
                    complete = True
                    if el not in row:
                        print(f'element {el} missing from row {row}')
                        complete = False

                    if complete:
                        if row['3_NONAMER'][g_e] != row['3_SPACER'][g_s]:
                            print(f"maths problem in {gene}: row['3_NONAMER'][{g_e}] != row['3_SPACER'][{g_s}]")
                        if row['3_SPACER'][g_e] != row['3_HEPTAMER'][g_s]:
                            print(f"maths problem in {gene}: row['3_SPACER'][{g_e}] != row['3_HEPTAMER'][{g_s}]")
                        if row['5_HEPTAMER'][g_e] != row['5_SPACER'][g_s]:
                            print(f"maths problem in {gene}: row['5_HEPTAMER'][{g_e}] != row['5_SPACER'][{g_s}]")
                        if row['5_SPACER'][g_e] != row['5_NONAMER'][g_s]:
                            print(f"maths problem in {gene}: row['5_SPACER'][{g_e}] != row['5_NONAMER'][{g_s}]")

                row['GENE'][g_s] -= GENE_WINDOW_SIZE_5
                row['GENE'][g_e] += GENE_WINDOW_SIZE_3
    return beds

def write_beds(beds, dir):
    bed_files = {}
    for refname in beds:
        for gene in beds[refname]:
            for el_type in beds[refname][gene]:
                row = beds[refname][gene][el_type]
                if el_type not in bed_files:
                    bed_files[el_type] = {}
                if refname not in bed_files[el_type]:
                    bed_files[el_type][refname] = []
                if 'el_type' in row:
                    del row['el_type']
                if 'seq' in row:
                    del row['seq']
                bed_files[el_type][refname].append(row)
                
    for el_type in bed_files:
        with open(os.path.join(dir, f'{el_type}.bed'), 'w', newline='') as fo:
            writer = csv.DictWriter(fo, delimiter='\t', fieldnames=['chain', 'start', 'end', 'gene', 'sense'])
            for refname in bed_files[el_type]:
                rows = bed_files[el_type][refname]
                rows.sort(key=lambda x: x['start'])
                writer.writerows(rows)

if __name__ == "__main__":
    read_beds(sense, dir, ref_seq)
