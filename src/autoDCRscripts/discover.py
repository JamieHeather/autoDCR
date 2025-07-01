import collections as coll
import os
import sys
from time import time
import json
from . import autoDCRfunctions as fxn

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# TODO can trim what;s not needed
def discover(in_path, out_path, species, loci, orientation, barcoding, dont_gzip, data_dir):
    """
    # TODO docstring
    :param in_file:
    :param out_path:
    :param species:
    :param loci:
    :param orientation:
    :param deletion_limit:
    :param cdr3_limit:
    :param dont_translate:
    :param dont_gzip:
    :return:
    """

    # Establish the necessary input data and parameters
    if not os.path.exists(in_path):
        raise IOError(f"Unable to find input file for TCR annotation ({in_path})!")
    # TODO sanity/presence check the input FQ (including a length check - give warning if too short)

    loci = fxn.check_features(loci, 'loci')
    vj_data = fxn.import_tcr_info(species, loci, 'JV', 'nt', data_dir)
    lc_data = fxn.import_tcr_info(species, loci, 'CL', 'nt', data_dir)

    headers = fxn.out_headers  # TODO del/fix?

    dcr_parameters = {'mode': 'discover',
                      'orientation': orientation,
                      'deletion_limit': 0,
                      'cdr3_limit': 0}

    discovery_parameters = {'min_read_#': 10,
                            'min_gene_%': 5}

    # TODO fix all this, just leftovers from annotate
    # # Determine where to save the results
    # if out_path == '[input-file-name].tsv':
    #     out_path = in_path[:in_path.rfind('.')].split('/')[-1] + '.tsv'
    # if not dont_gzip:
    #     out_path += '.gz'

    counts = coll.Counter()
    start = time()

    # Then loop through the input file, analysing TCRs as we go
    with fxn.opener(in_path, 'r') as in_file, fxn.opener(out_path, 'w') as out_file:

        genes = {'v': {}, 'j': {}}
        vjincidence = {}

        for read_id, seq, qual in fxn.readfq(in_file):
            # Pad empty quality scores for FASTA files
            if not qual:
                qual = ' ' * len(seq)

            counts['reads'] += 1

            # Figure out the relevant parts of the read for decombining, then search
            read, read_qual, bc, bc_qual = fxn.sort_read_bits(seq, qual, '')  # TODO add barcoding

            tcr_check = fxn.discover_genes(read, read_qual, vj_data, lc_data, dcr_parameters, headers)


            if tcr_check:

                region_check = 0
                for r in ['v', 'j']:

                    if r + '_call' not in tcr_check.keys():
                        continue
                    region_check += 1

                    bits = tcr_check[r + '_call'].split(',')

                    # Disregard inter-gene ambiguous calls
                    if len(list(set([x.split('*')[0] for x in bits]))) > 1:
                        continue

                    # Count cumulative fractional use of each gene
                    share = 1 / len(bits)
                    for call in bits:
                        gene, allele = call.split('*')
                        if gene not in genes[r]:
                            genes[r][gene] = {}
                        if allele not in genes[r][gene]:
                            genes[r][gene][allele] = 0
                        genes[r][gene][allele] += share

                # And count the frequency with which different V/J alleles are recombined together
                if region_check == 2:
                    v_bits = tcr_check['v_call'].split(',')
                    j_bits = tcr_check['j_call'].split(',')
                    share = 1 / len(v_bits)+len(j_bits)
                    for v in v_bits:
                        vgene, vallele = v.split('*')
                        if vgene not in vjincidence:
                            vjincidence[vgene] = {}
                        if vallele not in vjincidence[vgene]:
                            vjincidence[vgene][vallele] = {}
                        for j in j_bits:
                            jgene, jallele = j.split('*')
                            if jgene not in vjincidence[vgene][vallele]:
                                vjincidence[vgene][vallele][jgene] = {}
                            if jallele not in vjincidence[vgene][vallele][jgene]:
                                vjincidence[vgene][vallele][jgene][jallele] = 0
                            vjincidence[vgene][vallele][jgene][jallele] += share

    # Do some rudimentary haplotyping, by identifying confidently heterozygous J genes
    # TODO and C?
    haplotype_anchors = {}
    alleles = []
    for j in genes['j']:
        # # Apply simple total read (or share), and % frequency thresholds,  # TODO rm as replaced by fxn
        # filt = {allele: value for allele, value in genes['j'][j].items() if value > discovery_parameters['min_read_#']}
        # if len(filt) <= 1:
        #     continue
        # filt = {allele: value/sum(filt.values())*100 for allele, value in filt.items()}
        # filt = {allele: value for allele, value in filt.items() if value > discovery_parameters['min_gene_%']}
        # if len(filt) == 2:
        #     print(j, filt)
        #     haplotype_anchors[j] = filt[j]
        filtered = filter_alleles_usages(genes['j'][j], discovery_parameters)
        if filtered:
            for f in filtered:
                alleles.append('\t'.join([j, f.split('|')[0], str(filtered[f])]))
            if len(filtered) == 2:
                haplotype_anchors[j] = filtered


    # Having established the anchors, now go through the Vs and see if we can haplotype (potentially novel) alleles
    for v in genes['v']:
        filtered = filter_alleles_usages(genes['v'][v], discovery_parameters)
        if filtered:
            for f in filtered:
                alleles.append('\t'.join([v, f.split('|')[0], str(filtered[f])]))
            # print('-'*100, '\n', v)
            # print(filtered)
            # for vallele in filtered:
            #     print(vallele)
            #     if vallele in vjincidence[vgene]:
            #         for jgene in vjincidence[v][vallele]:
            #             if jgene in haplotype_anchors:
            #                 print('\t', jgene)
            #                 print('\t\t', vjincidence[v][vallele][jgene])
            #     else:
            #         print('\t\t-')

    # vjincidence[vgene][vallele][jgene][jallele]


    # print(haplotype_anchors)
    # print(alleles)

    if out_path == '[input-file-name].allele':
        out_prefix = in_path[:in_path.rfind('.')].split('/')[-1]
        out_path = out_prefix + '.allele'

    with open(out_path, 'w') as out_file:
        out_file.write('\n'.join(alleles) + '\n')

    with open(out_prefix + '-Genes.json', 'w') as out_file:
        json.dump(genes, out_file)

    with open(out_prefix + '-VJ.json', 'w') as out_file:
        json.dump(vjincidence, out_file)

    end = time()
    time_taken = end - start
    print("Took", str(round(time_taken, 2)), "seconds")


def filter_alleles_usages(usage_dict, parameters):
    # Apply simple total read (or share), and % frequency thresholds
    filt = {allele: value for allele, value in usage_dict.items() if value > parameters['min_read_#']}
    if len(filt) == 0:
        return
    filt = {allele: value / sum(filt.values()) * 100 for allele, value in filt.items()}
    filt = {allele: value for allele, value in filt.items() if value > parameters['min_gene_%']}
    if len(filt) <= 2:
        return filt
