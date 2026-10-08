#!/usr/bin/python
"""Define 5 kb gene promoters from a gene BED file.

Input : tab-delimited BED with columns chrom, start, end, gene_id, strand (header lines start with '#')
Output: <input minus its 4-letter extension>5kb_promoters.stranded.GENEBED.bed with columns
        chrom, promoter_start, promoter_end, strand, gene_id, gene_id

A promoter is the 5 kb upstream of the transcription start site on the gene's strand. Where neighbouring genes lie
closer than 5 kb the promoter is shortened (the space between two start sites is split at its mid-point), and a
promoter that would overlap a gene on the same strand is not annotated.

Usage: python promoters_from_gene_bed.py genes.bed
"""
import sys

if len(sys.argv) < 2:
    print('\nNot enough arguments entered.\nUsage: GENE BED\n')
    sys.exit()

print("\nImporting...\n")
import os
import csv
import operator
import itertools
import time

from operator import itemgetter
from itertools import groupby

def asint(s):
    try: return int(s), ''
    except ValueError: return sys.maxsize, s

print('Loading annotations from:\t' + sys.argv[1])

print('\nWriting to outfile:\t' + sys.argv[1][:-4] + '5kb_promoters.stranded.GENEBED.bed')
out_file  = open(sys.argv[1][:-4] + '5kb_promoters.stranded.GENEBED.bed', 'w')

chr_list = []
with open(sys.argv[1]) as GENEBED_archive:
    for ge in GENEBED_archive:
        if str(ge)[0] != '#':
            chr_list.append(ge.split('\t')[0])
chr_set = set(chr_list)

remove_list       = []
same_strand_start = []
for chrom in sorted(chr_set):
    print('\n\nWorking with chromosome ' + chrom + '\n')
    with open(sys.argv[1]) as GENEBED_archive:
        gene_dict = {}
        for gene in csv.reader(GENEBED_archive, delimiter = '\t'):
            if gene[0][0] != '#':
                if gene[0] == chrom:
                    try:
                        if gene[3] in gene_dict:
                            print('\n\t!!!\tBreak during gene library construction due to duplicate genes in annotation')
                            sys.exit()
                        else:
                            gene_dict[gene[3]] = [gene[0], gene[1], gene[2], gene[4], gene[3]]

                    except:
                        print('\n\n--- Error during gene library construction ---\n\n')
                        print('\n'+ str(gene) +'\n')
                        sys.exit()

        sorted_gene_list = [(k, gene_dict[k]) for k in sorted(gene_dict, key=asint)]

        chrom_dict = {}
        for gen in sorted_gene_list:
            if gen[1][1] not in chrom_dict:
                chrom_dict[gen[1][1]] = [gen]
            else:
                if chrom_dict[gen[1][1]][0][1][3] == gen[1][3]:

                    print('Recording gene ' + gen[1][4] + ' with promoter overlap the same strand')
                    same_strand_start.append(chrom_dict[gen[1][1]][0][1][4])
                    same_strand_start.append(gen[1][4])

                    chrom_dict[gen[1][1]].append(gen)
                else:
                    chrom_dict[gen[1][1]].append(gen)

        sorted_chrm_list = [(k, chrom_dict[k]) for k in sorted(chrom_dict, key=asint)]

        count         = 0
        switch        = 0
        stored_start  = 0
        stored_stop   = 0
        stored_chr    = '@'
        stored_strand = '+'
        strand_skip   = '@'
        for gene in sorted_chrm_list:
            print(gene)
            if switch == 0:
                stored_stop = gene[1][0][1][2]
                switch = 1
            print(stored_chr)
            print(stored_start)
            print(stored_stop)
            print(stored_strand)

            if len(gene[1]) > 1:
                print ('\n\n\n\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' + str(len(gene[1])) + ' DUPLICATES!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n\n')
                count = count + 1
                stored_dup_start = '@'
                stored_dup_stop  = '@'
                stored_dup_prom  = '@'
                for d_c, dup in enumerate(gene[1]):
                    print(dup)
                    print('--- On gene =\t' + dup[0] + ' ' + dup[1][3])
                    if d_c > 0:
                        if dup[1][4] in same_strand_start:
                            print('\t' + str(dup[1][4]) + ' is a start site same strand duplicate - mirroring promoter')
                            dup[1].append(stored_dup_prom)

                            if dup[1][3] == '+':
                                if int(dup[1][2]) > int(stored_stop) and int(dup[1][1]) < int(stored_start):
                                    print('\tRemoving duplicate startsite gene 2 as nested within previous gene\n')
                                    remove_list.append(dup[1][4])
                                elif int(dup[1][2]) > int(stored_stop):
                                    stored_stop = int(dup[1][2])

                            elif dup[1][3] == '-':
                                if int(dup[1][2]) < int(stored_start):
                                    print('Gene overlapping - same start site on the negative strand - annotation not recorded\n')
                                    remove_list.append(dup[1][4])

                                elif int(dup[1][2]) > int(stored_start):
                                    print('Previous gene overlapping - same start site on the negative strand - annotation not recorded\n')
                                    remove_list.append(stored_dup_name)
                                    try:
                                        dup[1][5] = int(dup[1][2]) + 5000
                                    except(IndexError):
                                        dup[1].append(int(dup[1][2]) + 5000)
                                    stored_start = int(dup[1][2])

                                else:
                                    print(str(dup[1][4]) + ' is a start site same strand duplicate - mirroring promoter\n')
                                    dup[1].append(stored_dup_prom)

                        else:
                            if dup[1][3] == '-':
                                if int(dup[1][2]) >= int(stored_dup_stop):
                                    dup[1].append(int(dup[1][2]) + 5000)
                                    print('Gene ' + dup[1][4] + ',\tstop:\t' + str(dup[1][2]) + ',\treplaces with ' + dup[1][2] + ' + 5000\n')
                                    stored_start = int(dup[1][2])
                                    stroed_stop  = int(dup[1][1])
                                else:
                                    stored_start = int(stored_dup_stop)
                                    stroed_stop  = int(dup[1][1])
                                    remove_list.append(dup[1][4])
                                    print('Promoter overlapping CDS - promoter not annotated')
                            elif dup[1][3] == '+':
                                if stored_strand == '-':
                                    print('\n\n\tNOODLE\n\n\n')
                                    dist = int((int(dup[1][1])-int(stored_start))/2)

                                    if dist <= 1:
                                        print('\nOverlapping gene - promoter not annotated')
                                        remove_list.append(dup[1][4])
                                        print('\nGene prior to overlapping pair also promoter not annotated due to overlap')
                                        remove_list.append(sorted_chrm_list[count-2][1][0][1][4])
                                        continue

                                    elif dist < 5000 :
                                        print('Duplicate start site gene - calculating mid distance for promoter')
                                        print('Gene ' + dup[1][4] + ',\tstart:\t' + str(dup[1][1]) + ',\treplaces with ' + dup[1][1] + ' - ' + str(int((int(dup[1][1])-int(stored_start))/2)))
                                        dup[1].append(int(int(dup[1][1]) - int((int(dup[1][1])-int(stored_start))/2)))

                                        print('Adjusting previous promoter stop to: ' + str(stored_start) + ' + ' +  str(int((int(dup[1][1])-int(stored_start))/2)))
                                        print(str(sorted_chrm_list[count-2][1][0][1]))
                                        sorted_chrm_list[count-2][1][0][1][5] = int(int(stored_start) + int((int(dup[1][1])-int(stored_start))/2))
                                        print(str(sorted_chrm_list[count-2][1][0][1]))

                                    else:
                                        print('Gene ' + dup[1][4] + ',\tstart:\t' + str(dup[1][1]) + ',\treplaces with ' + dup[1][1] + ' - 5000')
                                        dup[1].append(int(dup[1][1]) - 5000)
                                        print('shittywacket')

                                else:
                                    if int(dup[1][1]) - 5000 >= stored_stop:
                                        dup[1].append(int(dup[1][1]) - 5000)
                                        print('Gene ' + dup[1][4] + ',\tstart:\t' + str(dup[1][1]) + ',\treplaces with ' + dup[1][1] + ' - 5000\n')
                                        stored_start = int(dup[1][1])
                                        stroed_stop  = int(dup[1][2])
                                        print('shattywacket')
                                    else:
                                        dup[1].append(int(stored_stop))
                                        print('Gene ' + dup[1][4] + ',\tstart:\t' + str(dup[1][1]) + ',\treplaces with ' + str(stored_stop) + '\n')
                                        stored_start = int(dup[1][1])
                                        stroed_stop  = int(dup[1][2])

                                if stored_strand == dup[1][3] and int(dup[1][1]) >= int(stored_start):
                                    print('Removing promoter annotation as overlapping previous gene.')
                                    remove_list.append(dup[1][4])

                                if int(dup[1][2]) > int(stored_dup_start):
                                    print('Promoter of gen prior to ' + str(dup[1][4]) + ' overlapping CDS. Promoter not annotated.\n\n\n')
                                    remove_list.append(stored_dup_name)

                                if stored_dup_strand == '-' and int(stored_dup_start) > int(stored_start):
                                    stored_start = stored_dup_start
                                    if int(stored_dup_stop) > int(stored_stop):
                                        stored_stop = stored_dup_stop

                                print('\n' + str(stored_start)  + ' ' + str(stored_stop))
                                print(str(stored_dup_start)  + ' ' + str(stored_dup_stop)+'\n')

                    else:
                        stored_dup_name   = dup[1][4]
                        stored_dup_strand = dup[1][3]
                        if dup[1][0] == stored_chr:
                            if dup[1][3] == stored_strand:
                                if dup[1][3] == '-':
                                    if int(dup[1][1]) >= int(stored_start) + 5000:
                                        print('\nAdding preliminary promoter annotation to dup gene 1\n')
                                        dup[1].append(int(dup[1][2]) + 5000)
                                        stored_dup_start = int(dup[1][2])
                                        stored_dup_stop  = int(dup[1][1])
                                    else:
                                        print('\nGenes closer than 5kb - editting previous promoter\n')
                                        print('\tGene ' + sorted_chrm_list[count-2][1][0][1][4].strip() + ' stop:\t' + str(sorted_chrm_list[count-2][1][0][1][2]) + '\treplaces with ' +  str(dup[1][1]))
                                        print('\t' + str(sorted_chrm_list[count-2]))
                                        sorted_chrm_list[count-2][1][0][1][5] = int(dup[1][1])
                                        print('\t' + str(sorted_chrm_list[count-2]))
                                        print('\nAdding preliminary promoter annotation to dup gene 1\n')
                                        dup[1].append(int(dup[1][2]) + 5000)
                                        stored_dup_start = int(dup[1][2])
                                        stored_dup_stop  = int(dup[1][1])

                                        print('\n')
                                elif dup[1][3] == '+' and int(dup[1][1]) - 5000 <= 0:
                                    dup[1].append(0)
                                    print('Gene ' + dup[1][4] + ' start:\t' + str(dup[1][1]) + '\treplaces with 0')
                                    print('\n')

                                    stored_dup_start = int(dup[1][1])
                                    stored_dup_stop  = int(dup[1][2])
                                    stored_dup_prom  = int(dup[1][5])

                                else:
                                    if int(dup[1][1]) >= int(stored_stop) + 5000:
                                        print('Gene ' + dup[1][4] + ' start:\t' + str(dup[1][1]) + '\treplaces with ' + dup[1][1] + ' - 5000')
                                        dup[1].append(int(dup[1][1]) - 5000)
                                        print('shittywickoot')
                                        print('\n')

                                        stored_dup_start = int(dup[1][1])
                                        stored_dup_stop  = int(dup[1][2])
                                        stored_dup_prom  = int(dup[1][5])

                                    else:
                                        if int(dup[1][1]) > int(stored_start) and int(dup[1][2]) < int(stored_stop):
                                            print('Duplicate start site nested within previous gene - promoter not annotated')
                                            remove_list.append(dup[1][4])

                                            stored_dup_start = int(dup[1][1])
                                            stored_dup_stop  = int(dup[1][2])
                                        elif int(dup[1][1]) < int(stored_start):
                                            print('Duplicate start site overlaps previous gene on same strand - promoter not annotated')
                                            remove_list.append(dup[1][4])

                                            stored_dup_start = int(dup[1][1])
                                            stored_dup_stop  = int(dup[1][2])
                                        else:
                                            print('Gene ' + dup[1][4] + ' start:\t' + str(dup[1][1]) + '\treplaces with ' + dup[1][1] + ' - ' + str(int(dup[1][1]) - int(stored_stop)))
                                            dup[1].append(int(stored_stop))
                                            print('\n')

                                            stored_dup_start = int(dup[1][1])
                                            stored_dup_stop  = int(dup[1][2])
                                            stored_dup_prom  = int(dup[1][5])

                            elif dup[1][3] != stored_strand:
                                if dup[1][3] == '+':
                                    if int(dup[1][1]) < int(stored_start) and int(dup[1][2]) > int(stored_stop):
                                        print('\tRemoving duplicate start site gene 1 as overlapping/nested within previous gene\n')
                                        remove_list.append(dup[1][4])

                                        stored_dup_start = int(dup[1][1])
                                        stored_dup_stop  = int(dup[1][2])

                                        print('\tAlso removing previous gene on -ve strand due to overlap\n\n')
                                        remove_list.append(sorted_chrm_list[count-2][1][0][1][4])
                                    else:
                                        if (stored_start + 10000) > int(dup[1][1]) and stored_start < int(dup[1][2]):
                                            dist = int(dup[1][1]) - stored_start
                                            print('Stored start:\t' + str(stored_start))
                                            print('Querey start:\t' + str(dup[1][1]))
                                            print('The distance between the genes is:\t' + str(dist))

                                            if dist <= 1:
                                                stored_dup_start = int(dup[1][1])
                                                stored_dup_stop  = int(dup[1][2])

                                                print('\nOverlapping gene - promoter not annotated')
                                                remove_list.append(dup[1][4])
                                                print('\nGene prior to overlapping pair also promoter not annotated due to overlap')
                                                remove_list.append(sorted_chrm_list[count-2][1][0][1][4])
                                                continue
                                            else:
                                                print('Allowed promoter region is\t' + str(int(dist/2)))

                                                print('Gene ' + str(sorted_chrm_list[count-2][1][0][1][4]) + ',\tstop:\t' + str(sorted_chrm_list[count-2][1][0][1][5]) + '\treplaces with ' + str(sorted_chrm_list[count-2][1][0][1][2]) + ' + ' + str(int(dist/2)))
                                                sorted_chrm_list[count-2][1][0][1][5] = int(sorted_chrm_list[count-2][1][0][1][2]) + int(dist/2)
                                                print('Gene ' + dup[1][4].strip() + ',\tstart:\t' + str(dup[1][1]) + '\treplaces with ' + str(dup[1][1]) + ' - ' + str(int(dist/2)))
                                                dup[1].append(int(dup[1][1]) - int(dist/2))
                                                stored_dup_start = int(dup[1][1])
                                                stored_dup_stop  = int(dup[1][2])
                                                stored_dup_prom  = int(dup[1][5])
                                                print('\n')

                                        elif int(dup[1][1]) - 5000 <= 0:
                                            dup[1].append(0)
                                            print('Gene ' + dup[1][4].strip() + ',\tstop:\t' + str(dup[1][1]) + ',\treplaces with 0')
                                            stored_dup_start = int(dup[1][1])
                                            stored_dup_stop  = int(dup[1][2])
                                            stored_dup_prom  = int(dup[1][5])

                                            print('\n')

                                        else:
                                            if int(dup[1][1]) - 5000 >= int(stored_stop):
                                                print('Gene ' + dup[1][4] + ' start:\t' + str(dup[1][1]) + '\treplaces with ' + dup[1][1] + ' - 5000')
                                                gene[1][0][1].append(int(dup[1][1]) - 5000)
                                                stored_dup_start = int(dup[1][1])
                                                stored_dup_stop  = int(dup[1][2])
                                                stored_dup_prom  = int(dup[1][5])
                                                print('shittywooocket')
                                            else:
                                                print('\nGenes closer than 5kb - taking shorter promoter\n')
                                                print('Gene ' + dup[1][4] + ' start:\t' + str(dup[1][1]) + '\treplaces with ' + str(stored_stop))
                                                dup[1].append(int(stored_stop))

                                                stored_dup_start = int(dup[1][1])
                                                stored_dup_stop  = int(dup[1][2])
                                                stored_dup_prom  = int(dup[1][5])

                                else:
                                    dup[1].append(int(dup[1][2]) + 5000)
                                    print('Gene ' + dup[1][4] + ',\tstop:\t' + str(dup[1][1]) + ',\treplaces with ' + dup[1][2] + ' + 5000')
                                    print('\n')
                                    stored_dup_start = int(dup[1][2])
                                    stored_dup_stop  = int(dup[1][1])
                                    stored_dup_prom  = int(dup[1][5])

                        elif dup[1][0] != stored_chr:
                            print('First of two start site duplicates is at the start of a new scaffold...\n')
                            if dup[1][3] == '+':
                                if int(dup[1][1]) - 5000 <= 0:
                                    dup[1].append(0)
                                    print('Gene ' + dup[1][4].strip() + ',\tstop:\t' + str(dup[1][1]) + ',\treplaces with 0')
                                    stored_dup_start = int(dup[1][1])
                                    stored_dup_stop  = int(dup[1][2])
                                    stored_dup_prom  = int(dup[1][5])

                                    print('\n')
                                else:
                                    dup[1].append(int(dup[1][1]) - 5000)
                                    print('Gene ' + dup[1][4].strip() + ',\tstop:\t' + str(dup[1][1]) + ',\treplaces with ' + str(dup[1][1]) + ' - 5000\n')
                                    stored_dup_start = int(dup[1][1])
                                    stored_dup_stop  = int(dup[1][2])
                                    stored_dup_prom  = int(dup[1][5])

                            print(gene[1])

            else:
                print('--- On gene =\t' + gene[1][0][0] + ' ' + gene[1][0][1][3])
                print(gene)

                count = count + 1
                if gene[1][0][1][0] == stored_chr:
                    if gene[1][0][1][3] == stored_strand:
                        if gene[1][0][1][3] == '-':
                            if int(gene[1][0][1][2]) == int(stored_start) and int(gene[1][0][1][1]) == int(stored_stop):

                                print('\n' + str(gene))
                                print(stored_start)
                                print('\n\n\nExact duplicate gene on dif strand error -  manually consult...\n\n\n')
                                sys.exit()

                            print('Gene ' + gene[1][0][1][4].strip() + ' stop:\t' + str(gene[1][0][1][2]) + '\treplaces with ' + gene[1][0][1][2] + ' + 5000')
                            gene[1][0][1].append(int(gene[1][0][1][2]) + 5000)

                            if int(gene[1][0][1][1]) >= int(stored_start) + 5000:
                                print ('plink')
                                stored_start = int(gene[1][0][1][2])
                                stored_stop  = int(gene[1][0][1][1])
                            else:

                                print('\nGenes closer than 5kb - editting previous promoter\n')

                                if len(sorted_chrm_list[count-2][1]) == 1:

                                    if int(stored_start) > int(gene[1][0][1][2]) > int(stored_stop):
                                        print('Gene nested or overlapping previous gene on same strand - promoter not annotated\n')
                                        remove_list.append(gene[1][0][1][4])

                                    if int(gene[1][0][1][2]) >= int(stored_start) > int(gene[1][0][1][1]):
                                        print('Gene nested or overlapping previous gene on same strand - promoter not annotated\n')
                                        remove_list.append(gene[1][0][1][4])
                                        print('Previous gene nested or overlapping previous gene on same strand - promoter annotation removed\n')
                                        remove_list.append(sorted_chrm_list[count-2][1][0][1][4])
                                        stored_start = int(gene[1][0][1][2])
                                        stored_stop  = int(gene[1][0][1][1])

                                    else:
                                        print('tral')
                                        if gene[1][0][1][4] not in remove_list:
                                            print('\tGene ' + sorted_chrm_list[count-2][1][0][1][4] + ' stop:\t' + str(sorted_chrm_list[count-2][1][0][1][2]) + '\treplaces with ' +  gene[1][0][1][1])
                                            sorted_chrm_list[count-2][1][0][1][5] = int(gene[1][0][1][1])
                                            if sorted_chrm_list[count-2][1][0][1][4] not in remove_list:
                                                stored_start = int(gene[1][0][1][2])
                                                stored_stop  = int(gene[1][0][1][1])
                                elif len(sorted_chrm_list[count-2][1]) == 2:
                                    if sorted_chrm_list[count-2][1][1][1][3] == gene[1][0][1][3]:
                                        print('\tGene ' + sorted_chrm_list[count-2][1][1][1][4] + ' stop:\t' + str(sorted_chrm_list[count-2][1][1][1][2]) + '\t + replaces with ' +  gene[1][0][1][1])
                                        try:
                                            sorted_chrm_list[count-2][1][1][1][5] = int(gene[1][0][1][1])
                                        except(IndexError):
                                            print('\n\t!!! Missing annotation appended...\n')
                                            sorted_chrm_list[count-2][1][1][1].append(int(gene[1][0][1][1]))
                                    else:
                                        print('\tGene ' + sorted_chrm_list[count-2][1][0][1][4] + ' stop:\t' + str(sorted_chrm_list[count-2][1][0][1][2]) + '\t + replaces with ' +  gene[1][0][1][1])
                                        try:
                                            sorted_chrm_list[count-2][1][0][1][5] = int(gene[1][0][1][1])
                                        except(IndexError):
                                            print('\n\t!!! Missing annotation appended...\n')
                                            sorted_chrm_list[count-2][1][0][1].append(int(gene[1][0][1][1]))

                                    stored_start = int(gene[1][0][1][2])
                                    stored_stop  = int(gene[1][0][1][1])
                            print('\n')
                        elif gene[1][0][1][3] == '+':

                            if stored_start == 0:
                                print('prar')
                                gene[1][0][1].append(0)
                                print('Gene ' + gene[1][0][1][4].strip() + ' start:\t' + str(gene[1][0][1][1]) + '\treplaces with 0')
                                stored_start  = int(gene[1][0][1][1])
                                stored_stop   = int(gene[1][0][1][2])
                                stored_strand = gene[1][0][1][3]
                                print('\n')

                            else:

                                if int(gene[1][0][1][1]) - 5000 >= int(stored_stop):
                                    print('Gene ' + gene[1][0][1][4].strip() + ' start:\t' + str(gene[1][0][1][1]) + '\treplaces with ' + gene[1][0][1][1] + ' - 5000')
                                    gene[1][0][1].append(int(gene[1][0][1][1]) - 5000)
                                    stored_start = int(gene[1][0][1][1])
                                    stored_stop  = int(gene[1][0][1][2])
                                    print('shooottywicket')
                                    strand_skip = '@'
                                elif int(stored_stop) > int(gene[1][0][1][1]) > int(stored_start):
                                    print('Gene nested or overlapping previous gene on same strand - promoter not annotated\n')
                                    remove_list.append(gene[1][0][1][4])
                                else:
                                    print('\nGenes closer than 5kb - taking shorter promoter\n')

                                    print('Gene ' + gene[1][0][1][4].strip() + ' start:\t' + str(gene[1][0][1][1]) + '\treplaces with ' + gene[1][0][1][1] + ' - ' + str(int(gene[1][0][1][1]) - int(stored_stop)))
                                    gene[1][0][1].append(int(stored_stop))
                                    stored_start = int(gene[1][0][1][1])
                                    stored_stop  = int(gene[1][0][1][2])
                                    strand_skip = '@'
                                print('\n')
                            stored_strand = gene[1][0][1][3]

                    elif gene[1][0][1][3] != stored_strand:
                        print ('Current strand: '+ gene[1][0][1][3] + ' | Stored strand: ' + stored_strand)
                        if gene[1][0][1][3] == '+':
                            if int(stored_start + 10000) > int(gene[1][0][1][1]) and stored_start < int(gene[1][0][1][2]):
                                dist = int(gene[1][0][1][1]) - int(stored_start)
                                print('Stored start:\t' + str(stored_start))
                                print('Querey start:\t' + str(gene[1][0][1][1]))
                                print('The distance between the genes is:\t' + str(dist))

                                if dist <= 1:
                                    print('\nOverlapping gene - promoter not annotated')
                                    remove_list.append(gene[1][0][1][4])
                                    print('Also removing previous promoter annotation due to overlap')
                                    remove_list.append(sorted_chrm_list[count-2][1][0][1][4])
                                    stored_start = int(gene[1][0][1][1])
                                    stored_stop  = int(gene[1][0][1][2])
                                    print('\n')
                                    strand_skip  = '@'
                                else:
                                    print(sorted_chrm_list[count-2])
                                    print(str(gene))
                                    print('Allowed promoter region is\t' + str(int(dist/2)))
                                    if len(sorted_chrm_list[count-2][1]) == 1:
                                        print('Gene ' + str(sorted_chrm_list[count-2][1][0][1][4].strip()) + '  stop:  ' + str(sorted_chrm_list[count-2][1][0][1][5]) + '\treplaces with\t' + str(sorted_chrm_list[count-2][1][0][1][2]) + ' + ' + str(int(dist/2)))
                                        print('Gene ' + gene[1][0][1][4] + ' start:  ' + str(gene[1][0][1][1]) + '\treplaces with\t' + str(gene[1][0][1][1]) + ' - ' + str(int(dist/2)))
                                        gene[1][0][1].append(int(gene[1][0][1][1]) - int(dist/2))
                                        print(int(sorted_chrm_list[count-2][1][0][1][2]) + int(dist/2))
                                        sorted_chrm_list[count-2][1][0][1][5] = int(sorted_chrm_list[count-2][1][0][1][2]) + int(dist/2)
                                        stored_start = int(gene[1][0][1][1])
                                        stored_stop  = int(gene[1][0][1][2])
                                        print('\n')
                                    elif len(sorted_chrm_list[count-2][1]) == 2:
                                        if dup[1][3] != sorted_chrm_list[count-2][1][0][1][3]:
                                            try:
                                                print('Gene ' + str(sorted_chrm_list[count-2][1][0][1][4]) + '  stop:  ' + str(sorted_chrm_list[count-2][1][0][1][5]) + '\treplaces with\t' + str(sorted_chrm_list[count-2][1][0][1][2]) + ' + ' + str(int(dist/2)))
                                            except(IndexError):
                                                print('\n\t!!! Missing annotation appended to dup 1...\n')
                                                sorted_chrm_list[count-2][1][0][1].append(int(sorted_chrm_list[count-2][1][0][1][2]) + int(dist/2))
                                        else:
                                            try:
                                                print('Gene ' + str(sorted_chrm_list[count-2][1][1][1][4].strip()) + '  stop:  ' + str(sorted_chrm_list[count-2][1][1][1][5]) + '\treplaces with\t' + str(sorted_chrm_list[count-2][1][1][1][2]) + ' + ' + str(int(dist/2)))
                                            except(IndexError):
                                                print('\n\t!!! Missing annotation appended to dup 2...\n')
                                                sorted_chrm_list[count-2][1][1][1].append(int(sorted_chrm_list[count-2][1][1][1][2]) + int(dist/2))

                                        print('Gene ' + gene[1][0][1][4] + ' start:  ' + str(gene[1][0][1][1]) + '\treplaces with\t' + str(gene[1][0][1][1]) + ' - ' + str(int(dist/2)))
                                        gene[1][0][1].append(int(gene[1][0][1][1]) - int(dist/2))

                                        try:
                                            print(sorted_chrm_list[count-2][1][1][1][5])
                                        except:
                                            print('None stored')

                                        stored_start = int(gene[1][0][1][1])
                                        stored_stop  = int(gene[1][0][1][2])
                                        print('\n')
                                    strand_skip  = '@'
                            elif int(gene[1][0][1][1]) - 5000 <= 0:
                                if int(gene[1][0][1][1]) - int(stored_stop) > 0:
                                    gene[1][0][1].append(0)
                                    print('Gene ' + gene[1][0][1][4].strip() + ', start:\t' + str(gene[1][0][1][1]) + ',\treplaces with 0\n')
                                    stored_start = int(gene[1][0][1][1])
                                    print('\n')
                                else:
                                    remove_list.append(gene[1][0][1][4])
                                    print('Gene overlapping previous on negative strand - promoter not annotated\n')
                                strand_skip  = '@'
                            elif int(stored_start + 10000) > int(gene[1][0][1][1]) and int(stored_start) >= int(gene[1][0][1][2]):
                                print('Gene nested within another - promoter not annotated')
                                remove_list.append(gene[1][0][1][4])
                                strand_skip = stored_strand
                                print('\n')
                            else:
                                gene[1][0][1].append(int(gene[1][0][1][1]) - 5000)
                                print('Gene ' + gene[1][0][1][4].strip() + ', start:\t' + str(gene[1][0][1][1]) + ',\treplaces with ' + gene[1][0][1][1] + ' - 5000')
                                print('pop')
                                stored_start = int(gene[1][0][1][1])
                                stored_stop  = int(gene[1][0][1][2])
                                strand_skip  = '@'
                                print('\n')

                        else:
                            gene[1][0][1].append(int(gene[1][0][1][2]) + 5000)
                            print('Gene ' + gene[1][0][1][4].strip() + ', stop:\t' + str(gene[1][0][1][2]) + ',\treplaces with ' + gene[1][0][1][2] + ' + 5000')
                            if strand_skip == gene[1][0][1][3] and int(stored_start) + 5000 > int(gene[1][0][1][1]):

                                print('\nGene prior to nested genes closer than 5kb - editting previous promoter\n')

                                if len(sorted_chrm_list[count-2][1]) == 1:

                                    print('PUMP\n')
                                    print(stored_start)
                                    print(str(stored_stop) + '\n')
                                    print(str(gene[1][0][1]) + '\n')
                                    if int(stored_start) > int(gene[1][0][1][2]) > int(stored_stop):
                                        print('Gene nested or overlapping previous gene on same strand - promoter not annotated\n')
                                        remove_list.append(gene[1][0][1][4])

                                    if int(gene[1][0][1][2]) > int(stored_start) > int(gene[1][0][1][1]):
                                        print('Previous gene nested or overlapping previous gene on same strand - promoter annotation removed\n')
                                        remove_list.append(sorted_chrm_list[count-2][1][0][1][4])
                                        print(str(remove_list))

                                    else:
                                        print(strand_skip)
                                        print('peen')
                                        if sorted_chrm_list[count-2][1][0][1][4] in remove_list:
                                            print('\tGene ' + sorted_chrm_list[count-3][1][0][1][4] + ' stop:\t' + str(sorted_chrm_list[count-3][1][0][1][2]) + '\treplaces with ' +  gene[1][0][1][1])
                                            print('\t' + str(sorted_chrm_list[count-2]))
                                            print('peep')
                                            strand_skip = '@'
                                        elif sorted_chrm_list[count-3][1][0][1][4] not in remove_list:
                                            print('jeen')
                                            print(sorted_chrm_list[count-4][1][0][1][4])
                                            print(sorted_chrm_list[count-3][1][0][1][4])
                                            print(sorted_chrm_list[count-2][1][0][1][4])
                                            print(sorted_chrm_list[count-1][1][0][1][4])
                                            print('\tGene ' + sorted_chrm_list[count-3][1][0][1][4] + ' stop:\t' + str(sorted_chrm_list[count-3][1][0][1][2]) + '\treplaces with ' +  gene[1][0][1][1])
                                            sorted_chrm_list[count-3][1][0][1][5] = int(gene[1][0][1][1])
                                            strand_skip = '@'
                                        elif sorted_chrm_list[count-4][1][0][1][4] not in remove_list:
                                            print('leen')
                                            print('\tGene ' + sorted_chrm_list[count-3][1][0][1][4] + ' stop:\t' + str(sorted_chrm_list[count-4][1][0][1][2]) + '\treplaces with ' +  gene[1][0][1][1])
                                            sorted_chrm_list[count-4][1][0][1][5] = int(gene[1][0][1][1])
                                            strand_skip = '@'
                                        else:
                                            print('\n\n\nToooo much nesting - script must be manually tailored to the job')
                                            sys.exit()
                                elif len(sorted_chrm_list[count-2][1]) == 2:
                                    if sorted_chrm_list[count-2][1][1][1][3] == gene[1][0][1][3]:
                                        print('\tGene ' + sorted_chrm_list[count-2][1][1][1][4] + ' stop:\t' + str(sorted_chrm_list[count-2][1][1][1][2]) + '\t + replaces with ' +  gene[1][0][1][1])
                                        try:
                                            sorted_chrm_list[count-2][1][1][1][5] = int(gene[1][0][1][1])
                                        except(IndexError):
                                            print('\n\t!!! Missing annotation appended...\n')
                                            sorted_chrm_list[count-2][1][1][1].append(int(gene[1][0][1][1]))
                                    else:
                                        print('\tGene ' + sorted_chrm_list[count-2][1][0][1][4] + ' stop:\t' + str(sorted_chrm_list[count-2][1][0][1][2]) + '\t + replaces with ' +  gene[1][0][1][1])
                                        try:
                                            sorted_chrm_list[count-2][1][0][1][5] = int(gene[1][0][1][1])
                                        except(IndexError):
                                            print('\n\t!!! Missing annotation appended...\n')
                                            sorted_chrm_list[count-2][1][0][1].append(int(gene[1][0][1][1]))

                                print('pootpoot')
                                print('pootpoot')
                            print('\n')
                            stored_start = int(gene[1][0][1][2])
                            stored_stop  = int(gene[1][0][1][1])
                        stored_strand = gene[1][0][1][3]
                elif gene[1][0][1][0] != stored_chr:
                    print(gene[1][0][1][4])
                    if gene[1][0][1][3] == '-':
                        print('Gene ' + gene[1][0][1][4] + ', stop:\t' + str(gene[1][0][1][2]) + ',\treplaces with ' + gene[1][0][1][2] + ' + 5000')
                        gene[1][0][1].append(int(gene[1][0][1][2]) + 5000)
                        stored_start  = int(gene[1][0][1][2])
                        stored_stop   = int(gene[1][0][1][1])
                        stored_strand = gene[1][0][1][3]
                    elif gene[1][0][1][3] == '+' and int(gene[1][0][1][1]) - 5000 <= 0:
                        gene[1][0][1].append(0)
                        print('Gene ' + gene[1][0][1][4] + ', start:\t' + str(gene[1][0][1][1]) + ',\treplaces with 0')
                        stored_start = int(gene[1][0][1][1])
                        stored_stop  = int(gene[1][0][1][2])
                    else:
                        print('Gene ' + gene[1][0][1][4] + ', start:\t' + str(gene[1][0][1][1]) + ',\treplaces with '+ gene[1][0][1][1] + ' - 5000')
                        gene[1][0][1].append(int(gene[1][0][1][1]) - 5000)
                        stored_start = int(gene[1][0][1][1])
                        stored_stop  = int(gene[1][0][1][2])
                        print('shittywicket')

                    stored_chr = gene[1][0][1][0]

        print('\nWriting to:\t' + sys.argv[1][:-4] + '5kb_promoters.stranded.bed')

        for bit in sorted_chrm_list:
            if len(bit[1]) > 1:
                for res in bit[1]:
                    if res[1][4] not in remove_list:
                        print('Keeping: ' + str(res))
                        if res[1][3] == '+':
                            out_file.write(res[1][0] + '\t' + str(int(res[1][5])) + '\t' + str(res[1][1]) + '\t' + str(res[1][3]) + '\t' + str(res[0]) + '\t' + str(res[0]) + '\n')
                        elif res[1][3] == '-':
                            out_file.write(res[1][0] + '\t' + str(int(res[1][2])) + '\t' + str(res[1][5]) + '\t' + str(res[1][3]) + '\t' + str(res[0]) + '\t' + str(res[0]) + '\n')
                    else:
                        print('--- ' + str(res[1][4]) + ' removed' + '\t' + res[1][3] + '\t' + res[1][1] + '\t' + res[1][2])
            else:
                if bit[1][0][1][4] not in remove_list:
                    print('Keeping: ' + str(bit))
                    if bit[1][0][1][3] == '+':
                        out_file.write(bit[1][0][1][0] + '\t' + str(int(bit[1][0][1][5])) + '\t' + str(bit[1][0][1][1]) + '\t' + str(bit[1][0][1][3]) + '\t' + str(bit[1][0][1][4].strip()) + '\t'+ str(bit[1][0][1][4]) + '\n')

                    elif bit[1][0][1][3] == '-':
                        out_file.write(bit[1][0][1][0] + '\t' + str(int(bit[1][0][1][2])) + '\t' + str(bit[1][0][1][5]) + '\t' + str(bit[1][0][1][3]) + '\t' + str(bit[1][0][1][4].strip()) + '\t' + str(bit[1][0][1][4]) + '\n')
                else:
                    print('--- ' + str(bit[1][0][1][4]) + ' removed' + '\t' + str(bit[1][0][1][3]) + '\t' + str(bit[1][0][1][1]) + '\t' + str(bit[1][0][1][2]))
out_file.close()

with open(sys.argv[1][:-4] + 'same_strand+start.out', 'w') as out2:
    for out in same_strand_start:
        if out not in remove_list:
            out2.write(str(out) + '\n')
out2.close()
print('\n\n\nOutput closed\n')

