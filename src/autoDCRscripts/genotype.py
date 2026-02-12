
import os
import sys
import json
import collections as coll
import pandas as pd
from . import autoDCRfunctions as fxn
from . import referencedata as refdat

def genotype(in_file, mode, out_path, ref, subset_loci, allow_partial, read_threshold, data_dir):
    # TODO docstr

    if out_path == fxn.genotype_default:
        out_path = fxn.basename(in_file) + '_genotype.json'
    else:
        if not out_path.lower().endswith('.json'):
            out_path += '.json'

    # Read in the requisite data and establish the data structures
    df = pd.read_csv(in_file, sep='\t', compression='infer', keep_default_na=False)

    df_headers = [x for x in df]
    if allow_partial:
        expected_headers = fxn.out_headers
    else:
        expected_headers = fxn.out_headers + fxn.tag_coverage_headers
    if mode == 'full':
        expected_headers += fxn.full_feat_headers
    if len([x for x in expected_headers if x in df_headers]) < len(expected_headers):
        raise IOError("Unable to run genotype due to failure to detect the correct headers: "
                      "\nensure this in file is the product of autoDCR,"
                      " run with the appropriate options.")

    fields = ['v_call', 'j_call']
    if mode == 'full':
        fields.append('c_call')

    # Also check the relevent reference files are there, and read in relevant data
    loci = fxn.check_features(subset_loci, 'loci')
    ref_dir = os.path.join(data_dir, ref)
    if not os.path.exists(ref_dir):
        raise IOError("Cannot find relevant data directory: check name, and that autodcr refs completed correctly. ")

    pot_ambig_path = os.path.join(ref_dir, ''.join(loci) + '_JV_potential_ambiguity.json')
    if not os.path.exists(pot_ambig_path):
        raise IOError(f"Cannot find the 'potential_ambiguity' file in the reference directory specified. "
                      f"Please re-run autodcr refs and verify its presence. ")

    with open(pot_ambig_path, mode='r') as in_file:
        pot_ambig = json.load(in_file)

    # Establish a nested dictionary for the different fields
    genes = {}
    genes_unambig = {}
    geno_dat = {}

    for f in fields:
        genes[f] = {}
        genes_unambig[f] = {}
        geno_dat[f] = {}

    inter_gene_ambiguity = coll.Counter()

    for row in df.index:
        row_dat = df.loc[row]
        for f in fields:
            if row_dat[f]:

                # Skip rearrangements that didn't sequence to the end of this region
                if not allow_partial:
                    if str(row_dat['missed_nt_' + f[0]]) != '0':
                        continue

                if ',' not in row_dat[f]:
                    gene, allele = row_dat[f].split('*')
                    if gene not in genes_unambig[f]:
                        genes_unambig[f][gene] = coll.Counter()
                    if allele not in genes_unambig[f][gene]:
                        genes_unambig[f][gene][allele] = 0
                    genes_unambig[f][gene][allele] += 1

                bits = row_dat[f].split(',')

                # Disregard inter-gene ambiguous calls
                if len(list(set([x.split('*')[0] for x in bits]))) > 1:
                    for b in bits:
                        inter_gene_ambiguity[b] += 1/len(bits)
                    continue

                # Count cumulative fractional use of each gene
                share = 1/len(bits)
                for call in bits:
                    gene, allele = call.split('*')
                    if gene not in genes[f]:
                        genes[f][gene] = coll.Counter()
                    if allele not in genes[f][gene]:
                        genes[f][gene][allele] = 0
                    genes[f][gene][allele] += share

    # Go through those dictionaries and write out a prospective genotype file
    for f in fields:

        # First establish the unambiguous calls
        for g in genes_unambig[f]:
            sorted_alleles = genes_unambig[f][g].most_common()
            sorted_alleles = [x for x in sorted_alleles if x[1] >= read_threshold]

            if len(sorted_alleles) >= 1:
                # TODO do we need to add an additional check here - maybe for disregarding 2nd alleles below a certain % of the 1st?
                geno_dat[f][g] = {}
                for i in range(min(2, len(sorted_alleles))):
                    allele = sorted_alleles[i]
                    geno_dat[f][g][allele[0]] = allele[1]

        # Then top this up with missing genes from the ambiguous call box if needed
        for g in genes[f]:

            # Skip those falling below a specified read count threshold
            if sum(genes[f][g].values()) <= read_threshold:
                continue

            # Determine if this gene was observed in the unambiguous call; skip if wholly found or no additional info
            if g in geno_dat[f]:
                gene_seen = True
                if len(geno_dat[f][g]) == 2 or geno_dat[f][g].keys() == genes[f][g].keys():
                    continue
            else:
                gene_seen = False

            # If there is an unambiguous call for a given gene, cross-reference with the 'potential ambiguous' lists
            # to throw out ambiguity most parsimoniously explained by trimming of similar germline genes
            if gene_seen:
                tbd = [x for x in genes[f][g] if x not in geno_dat[f][g]]
                for seen in geno_dat[f][g]:
                    ambig_check_name = g + '*' + seen + '|' + f[0].upper()
                    if ambig_check_name in pot_ambig:
                        ambig_check_genes = pot_ambig[ambig_check_name]
                        for specific_pot_ambig in ambig_check_genes:
                            allele_id = specific_pot_ambig.split('|')[0].split('*')[1]
                            if allele_id in tbd:
                                tbd.pop(tbd.index(allele_id))

                if not tbd:
                    continue

            else:
                tbd = [x for x in genes[f][g]]
                if not tbd:
                    continue

            # If there are ambiguous calls not attributable to an unambiguous one, determine if they can be clustered
            if tbd:
                tbd = [x for x in tbd if genes[f][g][x] > read_threshold]
                if not tbd:
                    continue

                # Take the first sorted options within clusters where do-able
                tbd.sort()
                clustered = coll.defaultdict(list)
                clustered_count = coll.Counter()

                for t_item in tbd:
                    if not clustered:
                        clustered[t_item] = []
                        clustered_count[t_item] += genes[f][g][t_item]
                    else:
                        t_name = g + '*' + t_item + '|' + f[0].upper()
                        if t_name in pot_ambig:
                            ambig_check_genes = pot_ambig[t_name]
                            for specific_pot_ambig in ambig_check_genes:
                                allele_id = specific_pot_ambig.split('|')[0].split('*')[1]
                                if allele_id in clustered:
                                    clustered[allele_id].append(t_item)
                                    clustered_count[t_item] += genes[f][g][allele_id]
                                else:
                                    clustered[t_item] = []
                                    clustered_count[t_item] += genes[f][g][t_item]

                # Then add those to the final output genotype file
                if len(clustered) > 0:
                    if g not in geno_dat[f]:
                        geno_dat[f][g] = {}
                    clustered_order = sorted(clustered, key=lambda k: len(clustered[k]), reverse=True)
                    for co in clustered_order:
                        if len(geno_dat[f][g]) < 2:
                            geno_dat[f][g][co] = clustered_count[co]
                        else:
                            # Only take up to 2 alleles!
                            break

    # Do a final check to see if the inter-gene ambiguity was solved in other calls
    # checked_genes = []
    # for iga in inter_gene_ambiguity:
    #     iga_genes = [x.split('*')[0] for x in iga]
    #     for g in iga_genes:
    #         if g in checked_genes:
    #             continue
    #         f = g[3].lower() + '_call'
    #         if g not in geno_dat[f]:
    #             print(iga_genes, '!'*50)
    #         checked_genes.append(g)
    #         # TODO write to a log? do something? FIX

    with open(out_path, 'w') as out_file:
        json.dump(geno_dat, out_file)

    # If requested, produce a filtered reference
    if ref != fxn.genotype_ref_placeholder:
        ref_dir = os.path.join(data_dir, ref)
        if not os.path.exists(ref_dir):
            raise IOError(f"Requested reference directory ('{ref.upper()}') absent from data directory:"
                          f" check spelling (case sensitive) and/orrun 'autoDCR refs' and try again.")

        # Generate a plaintext list of alleles to subset
        alleles_to_keep = []
        for f in fields:
            for g in geno_dat[f]:
                for a in geno_dat[f][g]:
                    alleles_to_keep.append(g + '*' + a)

        # Determine the appropriate FASTA file(s) to subset
        loci = fxn.check_features(subset_loci, 'loci')

        if mode == 'vjcdr3':
            regions = 'JV'
        elif mode == 'full':
            raise IOError("Currently the genotype-specific referencing only works on the default 'vjcdr3' mode.")
        else:
            raise IOError(f"Unexpected mode detected: {mode}.")

        # TODO write out a p
        # Then take a genotype-specific subset of that reference FASTA file
        in_file_path = os.path.join(ref_dir, ''.join(loci) + '_' + regions + '.fasta')
        out_fasta_path = out_path.replace('.json', '.fasta')
        # print(in_file)
        with open(in_file_path, 'r') as in_file, open(out_fasta_path, 'w') as out_file:
            for read_id, seq, null in fxn.readfq(in_file):
                gene_allele, region = read_id.split('|')
                gene, allele = gene_allele.split('*')
                field = region.lower() + '_call'
                if gene in geno_dat[field]:
                    if allele in geno_dat[field][gene]:
                        out_file.write(fxn.fastafy(read_id, seq))

        refdat.get_reference_data('null', 20, loci, 'JV', True, False,
                                  False, out_fasta_path, fxn.basename(out_fasta_path).upper().replace('_', '-'),
                                  data_dir)
        # TODO get rid of unnecessary upper-ing?
