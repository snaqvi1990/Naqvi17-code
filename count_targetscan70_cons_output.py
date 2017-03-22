# given targetscan_70.pl output,
# count number all gene mirna interactions conserved between human and given species (by taxonomy id)
# can specify format of output: list or number of all conserved interactions by family, or matrix

import sys, subprocess, os, gzip
import numpy as np

def run():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-in_file',help='output of targetscan_70.pl',default='')
    parser.add_argument('-out_file',help='output fo;e',default='')
    parser.add_argument('-gene_info_file',help='Gene_Info.txt from TargetScan',default='')
    parser.add_argument('-list_num_mat',help='list, num, or mat, for format of output file',default='')
    parser.add_argument('-spec_cons',help='Taxonomy ID of species for conservation w/ human',default='')

    in_file, out_file, gene_info_file, list_or_num, spec_cons = sys.argv[1:]

    gene_mir_spec = {}
    gene_mir_cons = {}
    t_id2gene = {}
    allgenes = {}
    allmirs = {}

    f = open(args.gene_info_file,'r')
    next(f)
    for line in f:
        t_id, g_id, gene = line.strip().split()[:3]
        t_id2gene[t_id] = gene
        allgenes[gene] = 1
    f.close()


    f = open(args.in_file,'r')
    next(f)
    for line in f:
        t_id, mir, spec, msa_start, msa_end, utr_start, utr_end, group_num, site_type = line.strip().split()[:9]
        if t_id in t_id2gene:
            gene = t_id2gene[t_id]
        else: continue
        allmirs[mir] = 1
        if (mir[:3] != 'miR') & (mir[:3] != 'let'):
            print(line)
        if gene not in gene_mir_spec:
            gene_mir_spec[gene] = {}
        if mir not in gene_mir_spec[gene]:
            gene_mir_spec[gene][mir] = {}
        if spec not in gene_mir_spec[gene][mir]:
            gene_mir_spec[gene][mir][spec] = []
        if ((msa_start+':'+msa_end) not in gene_mir_spec[gene][mir][spec]):
            gene_mir_spec[gene][mir][spec].append(msa_start+':'+msa_end)
    f.close()

    f = open(args.out_file,'w')
    f.write('gene\tsite_density_chicken_control\n')
    for gene in allgenes:
        if gene not in gene_mir_spec:
            f.write(gene+'\tNA\n')
    	continue
        if gene not in gene_mir_cons:
            gene_mir_cons[gene] = {}
        for mir in gene_mir_spec[gene]:
            if spec_cons != 'all':
                if ('9606' in gene_mir_spec[gene][mir]) & (spec_cons in gene_mir_spec[gene][mir]): #human and other supplied species
                    if bool(set(gene_mir_spec[gene][mir]['9606']) & set(gene_mir_spec[gene][mir][spec_cons])): #if aligned site
                        gene_mir_cons[gene][mir] = 1
            elif spec_cons == 'all':
                if ('9606' in gene_mir_spec[gene][mir]) : #human
                    gene_mir_cons[gene][mir] = 1
        if (args.list_num_mat == 'num'):
            f.write(gene+'\t'+str(len(gene_mir_cons[gene]))+'\n')
        if (args.list_num_mat == 'list'):
            f.write(gene+'\t'+','.join(gene_mir_cons[gene].keys())+'\n')
    f.close()


    if (args.list_num_mat == 'mat'):
        f = open(args.out_file,'w')
        f.write('gene\t'+'\t'.join(allmirs.keys())+'\n')
        for gene in allgenes:
            str_toWrite = gene
            for mir in allmirs.keys():
    	    str_add = '\t0'
                if gene not in gene_mir_spec:
    		str_add = '\t0'
    	    elif mir not in gene_mir_spec[gene]:
                    str_add = '\t0'
    	    elif mir in gene_mir_spec[gene]:
    	        if ('9606' in gene_mir_spec[gene][mir]) & (spec_cons in gene_mir_spec[gene][mir]): #human and other supplied species
                        if bool(set(gene_mir_spec[gene][mir]['9606']) & set(gene_mir_spec[gene][mir][spec_cons])): #if aligned site
                            str_add = '\t1'
    	    str_toWrite = str_toWrite + str_add
    	f.write(str_toWrite+'\n')
        f.close()

run()
