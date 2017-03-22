######
#given a file of seed sequences in TargetScan format (tab separated: family name, seed, species, additional fields),
# outputs a defined number of ATCG content and CpG content matched kmers for each seed family
######

import sys, subprocess, os, gzip
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import random
import itertools

def readInFile(tscan_file):

    fams_seeds = {}
    fams_specs = {}
    f = open(tscan_file,'r')
    next(f)
    for line in f:
        if (len(line.strip().split())) == 3:
            fam, seed, spec = line.strip().split()
        if (len(line.strip().split())) == 6:
            fam, seed, spec, mir, mirseq, cons = line.strip().split()
        if (len(line.strip().split())) == 7:
            fam, seed, spec, mir, mirseq, cons, mirbase_id = line.strip().split()
        fams_seeds[fam] = seed
        if fam not in fams_specs:
            fams_specs[fam] = []
        fams_specs[fam].append(spec)

    return(fams_seeds,fams_specs)

def matchKmers(fams_seeds, k, n):
    bases = ['A','U','C','G']
    all_kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    all_kmers_stats = {}
    for kmer in all_kmers:
        all_kmers_stats[kmer] = {'A':kmer.count('A'),'U':kmer.count('U'),'C':kmer.count('C'),'G':kmer.count('G'),'CpG':kmer.count('CG')}

    kmers_used = {}
    fams_selected_kmers = {}
    fams_num_matched_kmers = {}
    for fam in fams_seeds:
        fam_a = fams_seeds[fam].count('A')
        fam_u = fams_seeds[fam].count('U')
        fam_c = fams_seeds[fam].count('C')
        fam_g = fams_seeds[fam].count('G')
        fam_cpg = fams_seeds[fam].count('CG')

        #choose all matched kmers
        matched_kmers = []
        for kmer in all_kmers:
            if (all_kmers_stats[kmer]['A']==fam_a) & (all_kmers_stats[kmer]['U']==fam_u) & (all_kmers_stats[kmer]['C']==fam_c) & (all_kmers_stats[kmer]['G']==fam_g) &(all_kmers_stats[kmer]['CpG']==fam_cpg): # if appropriately matched
                if (kmer not in fams_seeds.values()) & (kmer not in kmers_used):
                    matched_kmers.append(kmer)

        fams_num_matched_kmers[fam] = len(matched_kmers)
        #choose n matched kmers
        fams_selected_kmers[fam] = random.sample(matched_kmers,int(n),)
        for kmer in fams_selected_kmers[fam]: # add the kmers to the used pile
            kmers_used[kmer] = 1

    print(min(fams_num_matched_kmers.values()))

    return(fams_selected_kmers)

def writeFile(fams_selected_kmers,fams_seeds,fams_seqs,fams_specs,out_file):

    f = open(out_file,'w')
    for fam in fams_selected_kmers:
        for kmer in fams_selected_kmers[fam]:
            str_toWrite = list()
            str_toWrite[1:(len(kmer)+1)] = kmer
            str_toWrite = ''.join(str_toWrite)
            f.write('\t'.join([fam,kmer,';'.join(list(set(fams_specs[fam])))])+'\n')
    f.close()

def run():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-tscan_file',help='tscan file',default='')
    parser.add_argument('-n',help='number of kmers to generate per family',default='')
    parser.add_argument('-out_file',help='tscan file',default='')

    args = parser.parse_args()

    fams_seeds,fams_specs = readInFile(args.tscan_file)

    fams_selected_kmers = matchKmers(fams_seeds,7,args.n)

    writeFile(fams_selected_kmers,fams_seeds,fams_seqs,fams_specs,args.out_file)


run()
