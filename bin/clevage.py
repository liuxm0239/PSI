#!/usr/bin/env python
# coding: utf-8

import sys 
import pysam
import argparse
import pandas as pd 

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam',  help='BAM input ')
parser.add_argument('-c', '--cpg',  help='CpG sites input')
parser.add_argument('-o', '--output',  help='output file')
parser.add_argument('--id', default='s01',  help='sampleID')
parser.add_argument('--hyper', type=float, default=0.7,  help='hyper methylation threadhold')
parser.add_argument('--hypo', type=float, default=0.3,  help='hypo methylation threadhold')

args = parser.parse_args(None if sys.argv[1:] else ['--help'])

###############################################################
# Bam file handler
sam = pysam.AlignmentFile(args.bam, 'rb')
#    '/SlurmDatabase/Clinical/2023/0810_nipt/2023-07-26/s02_alignment/s022_Brecal/EX-02-3061-1A_S1.marked.BQSR.bam', 'rb')

# read in CpG sites
cpg = pd.read_csv(args.cpg, sep="\t", header=None)
cpg.columns = ['Chr', 'Start', 'End', 'Tumor_ave', 'Normal_ave', 'Strand']

# classfy CpG sites by using methylation levels
cpgH = cpg[(cpg['Tumor_ave'] > args.hyper)]
cpgL = cpg[(cpg['Tumor_ave'] < args.hypo)]
###############################################################

def clevage_site(sam, site, mstrand):
    """
    Notation:
    Site position should be carefully taken according to 'region' defination \
    in [pysam]( https://pysam.readthedocs.io/en/latest/glossary.html#term-region)
    
    If site in your input is 1-based, change it to 0-base according to pysam region defination
    only 5' end should be counted as match!
    """
    if not any( [mstrand  == 'F', mstrand  == 'R' ]) :
        raise ValueError("Error: strand of CpG site not provided")

    match = 0
    clevage = None

    depth = sam.count(site[0], site[1], site[2])

    if depth > 0:
        iter = sam.fetch(site[0], site[1], site[2])
        #iter = sam.fetch(site[0], site[2], site[2] + 1)
        
        # for sites in forward strand
        if mstrand  == 'F': 
            depth_F = 0    
            for read in iter:
                if read.is_forward:
                    depth_F += 1
                    if read.reference_start == site[1]:
                        match += 1
            if depth_F > 0 : clevage = match / depth_F
            
        # for sites in reverse strand
        if mstrand  == 'R':
            depth_R = 0
            for read in iter:
                if read.is_reverse :
                    depth_R += 1    
                    if read.reference_end == site[1]:
                        match += 1
            if depth_R > 0 : clevage = match / depth_R
                
    #clevage.append(depth)
    return  clevage


def get_window(site, mstrand):
    """
    Notation:
    Site position should be carefully taken according to 'region' defination \
    in [pysam]( https://pysam.readthedocs.io/en/latest/glossary.html#term-region)
    
    If site in your input is 1-based, change it to 0-base according to pysam region defination
    """
    if not any( [mstrand  == 'F', mstrand  == 'R' ]) :
        raise ValueError("Error: strand of CpG site not provided")

    if site[1] < 5 :
        raise ValueError("Error: position of CpG site two close to chromsome end")
    # for CpG sites in forward strand
    if mstrand  == 'F': 
        chrs = [site[0]] * 11
        starts = range(site[1] - 5, site[1] + 6)
        ends = range(site[2] - 5, site[2] + 6)
        sites = list(zip(chrs, starts, ends ))
    else:
        chrs = [site[0]] * 11
        starts = reversed(range(site[2] - 5, site[2] + 6))
        ends = reversed(range(site[2] - 4, site[2] + 7))
        sites = list(zip(chrs, starts, ends ))
    return sites    


def clevage_window(sam, site, mstrand):
    """
    Notation:
    Site position should be carefully taken according to 'region' defination \
    in [pysam]( https://pysam.readthedocs.io/en/latest/glossary.html#term-region)
    
    If site in your input is 1-based, change it to 0-base according to pysam region defination
    """
    
    if not any( [mstrand  == 'F', mstrand  == 'R' ]) :
        raise ValueError("Error: strand of CpG site not provided")

    # for CpG sites in forward strand
    sites = get_window(site, mstrand)
    win = []
    for s in sites:
        cle = clevage_site(sam, s, mstrand)
        win.append(cle)
    return win


###############################################################
dfh = pd.DataFrame(columns=['-5', '-4', '-3', '-2', '-1', 'C', 'G', '2', '3', '4', '5'])
for idx in cpgH.index:
    site = (cpgH['Chr'][idx], cpgH['Start'][idx], cpgH['End'][idx])
    #if sam.count(site[0], site[1], site[2]) > 0:
        #mstrand = cpgH['Strand'][idx]
    mstrand = 'F'
    win = clevage_window(sam, site, mstrand)
    dfh.loc[len(dfh)] = win
    mstrand = 'R'
    win = clevage_window(sam, site, mstrand)
    dfh.loc[len(dfh)] = win
print( '_'.join(['dfh', args.id]),  "\t".join( [str(i) for i in dfh.mean() ]), sep="\t" )

dfl = pd.DataFrame(columns=['-5', '-4', '-3', '-2', '-1', 'C', 'G', '2', '3', '4', '5'])
for idx in cpgL.index:
    site = (cpgL['Chr'][idx], cpgL['Start'][idx], cpgL['End'][idx])
    #if sam.count(site[0], site[1], site[2]) > 0:
        #mstrand = cpgL['Strand'][idx]
    mstrand = 'F'
    win = clevage_window(sam, site, mstrand)
    dfl.loc[len(dfl)] = win
    mstrand = 'R'
    win = clevage_window(sam, site, mstrand)
    dfl.loc[len(dfl)] = win
    #print(clevage_window(sam, site, mstrand))
print(  '_'.join(['dfl', args.id]),  "\t".join( [str(i) for i in dfl.mean() ]), sep="\t" )

#dfh.to_csv("./Cleavage_proportion_hyper.csv.gz", compression='gzip', chunksize=1000000)
#dfl.to_csv("./Cleavage_proportion_hypo.csv.gz", compression='gzip', chunksize=1000000)
###############################################################
