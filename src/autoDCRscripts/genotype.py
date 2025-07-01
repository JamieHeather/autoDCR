
import os
import sys
import json
import collections as coll
import pandas as pd
from . import autoDCRfunctions as fxn

def genotype(in_file, mode):

    # Read in the requisite data and establish the data structures
    df = pd.read_csv(in_file, sep='\t', compression='infer', keep_default_na=False)

    df_headers = [x for x in df]
    if len([x for x in fxn.out_headers]) < len(fxn.out_headers):
        raise IOError("Unable to run genotype: ensure this in file is the product of autoDCRscripts.")
    if mode != 'full':
        if len([x for x in fxn.full_feat_headers if x in df_headers]) == len(fxn.full_feat_headers):
            print("autoDCRscripts annotate 'full' mode output detected, switching mode.")
            mode = 'full'

    fields = ['v_call', 'j_call']
    if mode == 'full':
        fields.append('c_call')

    # Establish a nested dictionary for the different fields
    genes = {}
    for f in fields:
        genes[f] = {}

    for row in df.index:
        row_dat = df.loc[row]
        for f in fields:
            if row_dat[f]:
                bits = row_dat[f].split(',')

                # Disregard inter-gene ambiguous calls
                if len(list(set([x.split('*')[0] for x in bits]))) > 1:
                    print('*', bits)
                    continue

                # Count cumulative fractional use of each gene
                share = 1/len(bits)
                for call in bits:
                    gene, allele = call.split('*')
                    if gene not in genes[f]:
                        genes[f][gene] = {}
                    if allele not in genes[f][gene]:
                        genes[f][gene][allele] = 0
                    genes[f][gene][allele] += share


    with open('test.json', 'w') as out_file:
        json.dump(genes, out_file)
