#!/usr/bin/env python
# coding: utf-8

import sys 
import os 
import pysam
import io
import argparse
import pandas as pd 

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam',  help='BAM input ')
parser.add_argument('-c', '--cpg',  help='CpG sites input')
parser.add_argument('-o', '--output',  help='output file')
parser.add_argument('--id', default='s01',  help='sampleID')

args = parser.parse_args(None if sys.argv[1:] else ['--help'])

###############################################################

def motif_count(sam, site):
    """
    Notation:
    Site position should be carefully taken according to 'region' defination \
    in [pysam]( https://pysam.readthedocs.io/en/latest/glossary.html#term-region)
    
    If site in your input is 1-based, change it to 0-base according to pysam region defination
    only 5' end should be counted as match!
    """
    cgn = 0
    ncg = 0
    
    depth = sam.count(site[0], site[1] - 2, site[2] + 2)
    if depth > 0:
        iter = sam.fetch(site[0], site[1] - 2, site[2] + 2)
        
        for read in iter:
            #if read.mapping_quality < 20 :
            #    continue
            if any([read.is_forward and read.reference_start == site[1], read.is_reverse and read.reference_end == site[2]+1 ] ):
            #if read.is_forward and read.reference_start == site[1]:
                cgn += 1
            if any([read.is_forward and read.reference_start == site[1]-1, read.is_reverse and  read.reference_end == site[2] + 2] ):
            #if read.is_forward and read.reference_start == site[1]-1:
                ncg += 1

    return [cgn, ncg]  

def clevage_site(sam, site):
    """
    Notation:
    Site position should be carefully taken according to 'region' defination \
    in [pysam]( https://pysam.readthedocs.io/en/latest/glossary.html#term-region)
    
    If site in your input is 1-based, change it to 0-base according to pysam region defination
    only 5' end should be counted as match!
    """
    #if not any( [mstrand  == 'F', mstrand  == 'R' ]) :
    #    raise ValueError("Error: strand of CpG site not provided")

    clevage = [None, None] # [clevage_Forward, clevage_Reverse]

    depth = sam.count(site[0], site[1], site[2])

    if depth > 0:
        iter = sam.fetch(site[0], site[1], site[2])
        #iter = sam.fetch(site[0], site[2], site[2] + 1)
        
        depth_F = 0    
        depth_R = 0
        match_F = 0
        match_R = 0

        for read in iter:
            if read.is_forward:
                # for sites in forward strand
                depth_F += 1
                if read.reference_start == site[1]:
                    match_F += 1
            else:
                # for sites in reverse strand
                depth_R += 1
                if read.reference_end == site[1]:
                    match_R += 1
        if depth_F > 0 : clevage[0] = match_F / depth_F
        if depth_R > 0 : clevage[1] = match_R / depth_R

    #clevage.append(depth)
    #if mstrand  == 'F':
    #    return  clevage[0]
    #esle:
    #    return  clevage[1]
    return  clevage

def get_window(site):
    """
    Notation:
    Site position should be carefully taken according to 'region' defination \
    in [pysam]( https://pysam.readthedocs.io/en/latest/glossary.html#term-region)
    
    If site in your input is 1-based, change it to 0-base according to pysam region defination
    """
    #if not any( [mstrand  == 'F', mstrand  == 'R' ]) :
    #    raise ValueError("Error: strand of CpG site not provided")

    if site[1] < 5 :
        raise ValueError("Error: position of CpG site two close to chromsome end")

    chrs = [site[0]] * 12
    starts = range(site[1] - 5, site[1] + 7)
    ends = range(site[2] - 5, site[2] + 7)
    sites = list(zip(chrs, starts, ends ))
        
    # for CpG sites in forward strand
    #if mstrand  == 'F':
    #    chrs = [site[0]] * 12
    #    starts = range(site[1] - 5, site[1] + 7)
    #    ends = range(site[2] - 5, site[2] + 7)
    #    sites = list(zip(chrs, starts, ends ))
    #else:
    #    chrs = [site[0]] * 11
    #    starts = reversed(range(site[2] - 5, site[2] + 6))
    #    ends = reversed(range(site[2] - 4, site[2] + 7))
    #    sites = list(zip(chrs, starts, ends ))
    return sites  
    
    
def clevage_window(sam, site):
    """
    Notation:
    Site position should be carefully taken according to 'region' defination \
    in [pysam]( https://pysam.readthedocs.io/en/latest/glossary.html#term-region)
    
    If site in your input is 1-based, change it to 0-base according to pysam region defination
    """
    
    #if not any( [mstrand  == 'F', mstrand  == 'R' ]) :
    #    raise ValueError("Error: strand of CpG site not provided")

    sites = get_window(site)
    win = [clevage_site(sam, s) for s in sites]
    
    return win


################################################################
# Bam file handler
sam = pysam.AlignmentFile(args.bam, 'rb')
#    '/SlurmDatabase/Clinical/2023/0810_nipt/2023-07-26/s02_alignment/s022_Brecal/EX-02-3061-1A_S1.marked.BQSR.bam', 'rb')

# read in CpG sites
cpg = pd.read_csv(args.cpg, sep="\t", header=None)
#cpg.columns = ['Chr', 'Start', 'End', 'Tumor_ave', 'Normal_ave', 'Strand']
cpg.columns = ['Chr', 'Start', 'End']


###############################################################
CGN = 0
NCG = 0
for idx in cpg.index:
    site = (cpg['Chr'][idx], cpg['Start'][idx], cpg['End'][idx])
    count_s = motif_count(sam, site)
    CGN += count_s[0]
    NCG += count_s[1]

print(args.id, " CpG sites CGN/NCG ratio is:", str(CGN), str(NCG), str(CGN/NCG), sep='\t')

#CGN = 0
#NCG = 0
#for idx in cpgL.index:
#    site = (cpgL['Chr'][idx], cpgL['Start'][idx], cpgL['End'][idx])
#    count_s = motif_count(sam, site)
#    CGN += count_s[0]
#    NCG += count_s[1]
#
#print("Hypo CpG sites CGN/NCG ratio is:", str(CGN), str(NCG), str(CGN/NCG), sep='\t')


