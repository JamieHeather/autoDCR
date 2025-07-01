# -*- coding: utf-8 -*-

import collections as coll
import datetime
import os
import pandas as pd
import re
import shutil
import sys
from IMGTgeneDL import IMGTgeneDL
from time import strftime, localtime
from . import autoDCRfunctions as fxn


# TODO add name option, to allow multiple different references per species!
def get_reference_data(species, slide_len, loci, regions, skip_download, novel, protein, data_dir):
    """
    # TODO docstring
    :param species:
    :param slide_len:
    :param loci:
    :param regions:
    :param skip_download:
    :param novel:
    :param protein:
    :param data_dir:
    :return:
    # TODO docstr
    """

    species = species.upper()
    species_dir = os.path.join(data_dir, species)
    loci = fxn.check_features(loci, 'loci')
    regions = fxn.check_features(regions, 'regions')
    full_regions = [fxn.regions[x] for x in regions]

    regexes = fxn.read_json(os.path.join(data_dir, 'regexes.json'))

    if novel and skip_download:
        raise IOError("'novel' and skip_download' flags detected, which are incompatible: remove one and try again.")

    overlap_denominator = 2
    if protein:
        # In order for the protein version to be run, the corresponding default nt VJ version must have been run
        species_dir_fl = os.listdir(species_dir)
        needed_prefix = fxn.format_ref_prefixes(loci, 'JV')
        for needed_fl in [needed_prefix + '.fasta', needed_prefix + '.translate', 'J-region-motifs.tsv']:
            if needed_fl not in species_dir_fl:
                raise IOError(f"File required for protein references missing ({needed_fl}): "
                              f"please run default nucleotide 'autoDCRscripts refs' for this species first.")

        skip_download = True
        overlap_denominator = 4

    # If needed, download the raw data
    if not skip_download or species not in os.listdir(data_dir):
        print(f"Downloading TCR data for species '{species}' via IMGTgeneDL...")
        download_reference_data(species, data_dir)
        j_motif_data, j_positions = get_j_motif_data(species_dir)

        # If requested, additionally download novel alleles (human only)
        if novel:
            if species != 'HUMAN':
                raise IOError("The 'novel' flag can only be used when species is set to HUMAN.")
            else:
                print("Downloading additional novel alleles...")
                # TODO make novel allele thresholds controllable from input arguments?
                j_motif_data, j_positions = download_novel_alleles(2, 3, species_dir,
                                                                   j_motif_data, j_positions)
                j_motif_data.to_csv(os.path.join(species_dir, 'J-region-motifs.tsv'), sep = '\t', index=False)
    else:
        print("Skipping data download, proceeding to autoDCRscripts file generation...")
        j_motif_data, j_positions = get_j_motif_data(species_dir)



    # TODO make a flag to allow addition of 'full' regions, without over-writing! maybe to get all at once?

    # Then parse the FASTA files for the selected loci, generating the desired files
    # TODO re-add version?
    log_str = 'Running autoDCRscripts refs on:\t ' + \
              datetime.datetime.today().date().isoformat() + ' ' + datetime.datetime.today().time().isoformat() + r
    log_str += 'Sliding window length used:\t' + str(slide_len) + r + lb
    if slide_len != 20 and not protein:
        print("NB: non-default sliding window length used.")

    # Run here for non-specific (all in one file)
    tiles = coll.defaultdict(list)
    slides = coll.defaultdict(list)
    slides_count = coll.Counter()
    full_genes = coll.defaultdict()
    functionalities = coll.defaultdict()
    partials = coll.defaultdict()

    # TODO make possible for custom outpaths - particularly w.r.t. making donor-specific folders (to do discover/personalised annotate in one)
    out_prefix = os.path.join(species_dir, fxn.format_ref_prefixes(loci, regions))
    if protein:
        out_prefix += '_AA'

    for locus in loci:
        locus_file = os.path.join(species_dir, 'TR' + locus + '.fasta')

        with (open(locus_file, 'r') as in_file):
            for readid, seq, qual in fxn.readfq(in_file):
                region = readid.split('~')[1]
                gene, func, partial = get_gene(readid)
                locus = gene[2]
                # Append the region type to the gene ID, to allow e.g. L/V disambiguation or non-std names
                gene += '|' + fxn.regions[region]

                if region in full_regions and \
                        ((locus in loci) or (locus == 'D' and 'A' in loci)):

                    if protein:
                        seq = translate_germlines(gene, region, seq, j_motif_data)  # TODO

                    full_genes[gene] = seq.upper()

                    window = full_genes[gene]
                    functionalities[gene] = func
                    partials[gene] = partial

                    # Take sliding windows across the full length genes, overlapping by half their length
                    for i in range(0, len(window) - slide_len + 1, int(slide_len / overlap_denominator)):

                        slide = window[i:i + slide_len]
                        if len(slide) == slide_len:
                            tiles[gene].append(slide)
                            slides[slide].append(gene)
                            slides_count[slide] += 1

    # Annotate all tags with all genes they appear in (regardless of window!)
    for tag in slides:
        for gene in full_genes:
            if tag in full_genes[gene]:
                if gene not in slides[tag]:
                    slides[tag].append(gene)

    # Determine the 'jump', and start to prepare the final out strings to be written
    # Jump = how many nucleotides are there between the start of the tag and the relevant recombined edge of the gene
    # I.e. V jumps = position in gene running 5'-3', J jumps = negative values relative to 3' edge
    jumps = {}
    all_tags = list(slides.keys())
    all_tags.sort()
    out_str = []
    manually_verify = False
    for tag in all_tags:

        # If only one gene/allele covered by a tag, take its jump
        # If multiple, output a list of each corresponding jump
        # Note that many will be the same if for multiple alleles of same gene
        relevant_genes = slides[tag]

        if len(relevant_genes) == 0:
            raise IOError("Somehow a tag has been generated that covers no genes: " + tag)

        elif len(relevant_genes) == 1:
            tag_index = full_genes[relevant_genes[0]].index(tag)
            jumps = str(get_jump(full_genes, relevant_genes[0], tag_index))
            tag_names = relevant_genes[0]

        else:
            relevant_genes.sort()
            tag_indices = [full_genes[x].index(tag) for x in relevant_genes]
            jumps = ','.join([str(get_jump(full_genes, relevant_genes[x], tag_indices[x]))
                              for x in range(len(tag_indices))])
            tag_names = ','.join(relevant_genes)

        # Add to output, ensuring each tag has the information of all genes covered and their respective jumps
        out_str.append('\t'.join([tag, jumps, tag_names]))

    # Sanity check that this contains both (and only) properly named V & J genes
    all_genes = list(full_genes.keys())
    all_genes.sort()

    log_str += ('Number genes covered:\t' + str(len(all_genes)) + r +
                'Number tags made:\t' + str(len(out_str)) + r + lb)

    # Then write out tags, and begin to look for translation parameters
    with open(out_prefix + '.tags', 'w') as out_file:
        out_file.write('\n'.join(out_str))

    if regions == ['J', 'V']:
        if not protein:
            tag_str, log_str = find_junction_residues(all_genes, full_genes, functionalities,
                                                      partials, log_str, regexes, j_positions)

            # Finally output the tag and log files, containing the details of this specific run, and the compiled FASTA reads
            with open(out_prefix + '.translate', 'w') as out_file:
                out_file.write(tag_str)

            with open(out_prefix + '.log', 'w') as out_file:
                out_file.write(log_str)

        else:
            shutil.copy(out_prefix.replace('_AA', '') + '.translate', out_prefix + '.translate')

    with open(out_prefix + '.fasta', 'w') as out_file:
        for g in full_genes:
            out_file.write(fastafy(g, full_genes[g]))


def download_novel_alleles(study_threshold, donor_threshold, species_dir, j_motif_df, j_pos):
    """
    # TODO docstr, point out it's from stitchrdl v0.3.0
    :param study_threshold: int of # of studies a given allele must be found in to be retained
    :param donor_threshold: int of # of donors a given allele must be found in to be retained, summed across all studies
    This function grabs novel alleles from the repo where I collate them, and reads them into the additional-genes file
    """
    import requests

    novel_repo_url = 'https://api.github.com/repos/JamieHeather/novel-tcr-alleles/contents/'
    summary_prefix = 'novel-TCR-alleles-'

    # Identify the current novel-tcr-allele file
    response = requests.get(novel_repo_url, headers={})
    novel_file_name = ''
    if response.status_code == 200:
        files = response.json()
        matching_files = [x for x in files if x['name'].startswith(summary_prefix)]
        if len(matching_files) == 1:
            novel_file_name = matching_files[0]['name']
            novel_file_url = matching_files[0]['download_url']

    if not novel_file_name:
        raise IOError("Unable to locate a suitable summary novel TCR allele file name in the GitHub repository. ")

    # If found, download it
    tsv_path = os.path.join(species_dir, novel_file_name)
    response = requests.get(novel_file_url, headers={})

    if response.status_code == 200:
        with open(tsv_path, 'wb') as out_file:
            out_file.write(response.content)
    else:
        raise IOError("Failed to download the novel TCR allele data. ")

    # Then read in and whittle down to those entries that meet select criteria
    novel = pd.read_csv(tsv_path, sep='\t')
    novel_out_fasta = os.path.join(species_dir, novel_file_name.replace('.tsv', '.fasta'))
    additional_motifs = []

    with open(novel_out_fasta, 'w') as out_file:
        for row in novel.index:
            row_bits = novel.loc[row]
            if pd.isna(row_bits['Notes']):
                notes = ''
            else:
                notes = str(row_bits['Notes'])
            func = '?'

            if 'Stop codon' in notes:
                func = 'P'

            # Disregard alleles found in fewer than the threshold # of studies
            if row_bits['Number-Datasets-In'] < study_threshold:
                continue
            # Disregard alleles found in fewer than the threshold # of donors
            elif row_bits['Number-Donors-In'] < donor_threshold:
                continue
            # If it's already added to IMGT, skip
            elif isinstance(row_bits['IMGT-ID'], str):
                continue
            # Same for if it's a shorter version
            elif 'Shorter version of IMGT' in notes:
                continue

            # Determine whether the allele has a valid name
            allele_id = ''
            if row_bits['Standard-ID'].startswith('TR'):
                allele_id = row_bits['Standard-ID']
            else:
                # If not, borrow one of the names from the paper(s) where it was discovered
                studies = [x for x in novel if '-Name' in x]
                for study in studies:
                    if not pd.isna(row_bits[study]):
                        allele_id = row_bits[study] + '-' + study.replace('-Name', '')
                        if allele_id[:2] != 'TR' or '*' not in allele_id:
                            allele_id = row_bits['Gene'] + '*' + allele_id

            if not allele_id:
                raise IOError("Unable to determine a possible allele ID for sequence with data:\n" + str(row_bits))

            else:
                header_bits = [novel_file_name.replace('.tsv', ''), allele_id, 'Homo sapiens', func,
                                   '', '', str(len(row_bits['Ungapped-Sequence'])) + ' nt',
                                   '', '', '', '', '', '', '', notes, '~' + region_key[allele_id[3]]]
                header = '|'.join(header_bits)
                out_file.write(fxn.fastafy(header, row_bits['Ungapped-Sequence']))

                # Fill in the J motif files as we go
                if header_bits[-1] == '~JOINING':
                    motif_row = IMGTgeneDL.determine_j_motifs(header_bits, row_bits['Ungapped-Sequence'])
                    motif_row[-1] = int(motif_row[-1])
                    j_motif_df.loc[len(j_motif_df)] = motif_row
                    j_pos[allele_id] = motif_row[-1]

    # Having written out the selected novel alleles FASTA, filter those reads into the relevant locus files
    with (open(novel_out_fasta, 'r') as in_file,
          open(os.path.join(species_dir, 'TRA.fasta'), 'a') as tra_out_file,
          open(os.path.join(species_dir, 'TRB.fasta'), 'a') as trb_out_file,
          open(os.path.join(species_dir, 'TRG.fasta'), 'a') as trg_out_file,
          open(os.path.join(species_dir, 'TRD.fasta'), 'a') as trd_out_file):
        for header, read, null in fxn.readfq(in_file):
            bits = header.split('|')
            locus = bits[1][2]
            vars()['tr' + locus.lower() + '_out_file'].write(fxn.fastafy(header, read))
            if locus == 'D':
                if bits[1][3] == 'V':
                    tra_out_file.write(fxn.fastafy(header, read))

    # TODO here also need to update the J region motifs file
    # TODO also the j_pos variable
    # Finally return the updated J motif data structures

    return j_motif_df, j_pos


def download_reference_data(dl_species, data_dir):
    """
    :param dl_species:
    :param data_dir:
    :return:
    """
    # TODO docstring
    dl_species = dl_species.upper()

    # Use IMGTgeneDL's 'stitchr' mode to download the relevant data for this species
    # Use a mocked up sys.argv to provide the necessary CLI inputs to argparse
    sysargvbackup = sys.argv
    sys.argv = ['IMGTgeneDL', '-m', 'stitchr', '-n', '-s', dl_species]  # Replace with actual arguments

    try:
        IMGTgeneDL.main()  # Replace with how you would normally call the package code

    except Exception:
        raise IOError("Unable to download TCR data via IMGTgeneDL!")

    finally:
        # Restore the backed up sys.argv
        sys.argv = sysargvbackup

    # Assuming that got downloaded appropriately, move it to the data directory
    if dl_species in os.listdir(os.getcwd()):
        # Delete older entries
        if dl_species in os.listdir(data_dir):
            print("Note: a directory for this species exists: overwriting.")
            shutil.rmtree(os.path.join(data_dir, dl_species))
        shutil.move(os.path.join(os.getcwd(), dl_species), os.path.join(data_dir, dl_species))
    else:
        raise IOError(f"Expected IMGTgeneDL directory for {dl_species} not detected after download!")

    # TODO add novel allele DL functionality?


def get_gene(string):
    """
    :param string: a string containing an IMGT gene name - likely a FASTA header, or something derived from it
    :return: the IMGT gene/allele, and its corresponding functionality and partiality information
    """
    bits = string.split('|')
    return bits[1], bits[3], bits[13]


def fastafy(gene, seq_line):
    """
    :param gene: Gene symbol, extracted from the read id
    :param seq_line: Total protein primary sequence, extracted from input FASTA/generated by in silico splicing
    :return: An output-compatible FASTA entry ready for writing to file
    """
    return ">" + gene + "\n" + seq_line + "\n"


def translate(seq):
    """
    :param seq: Nucleotide sequence
    :return: Translated nucleotide sequence, with overhanging 3' (J) residues trimmed
    """
    protein = ""
    # Trim sequence length to a multiple of 3 (which in this situation should just be removing the terminal J residue)
    difference = len(seq) % 3
    if difference != 0:
        seq = seq[:-difference]
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        protein += fxn.codons[codon]
    return protein


def find_cys(regex, search_str):
    """
    :param regex: str describing a regular expression to search for
    :param search_str: str of a sequence to search within
    :return: the index of the putative conserved Cys reg based on the last instance of that given regex
    """
    hits = [search_str.rfind(x) + len(x) - 1 for x in re.findall(regex, search_str)]
    if hits:
        return hits[-1]
    else:
        return


def get_jump(full_genes, gene_name, tag_number):
    """
    Determine the jump value for proper translation of a rearrangement using a particular gene
    :param gene_name: IMGT gene name, as caught in the 'relevant_genes' variable
    :param tag_number: position of the index of that tag in its source gene
    :return:
    """
    if gene_name[3] in ['V', 'L']:
        return tag_number
    elif gene_name[3] in ['J', 'C']:
        return -abs(tag_number - len(full_genes[gene_name]))
    else:
        raise IOError("Unrecognised gene type detected: " + gene_name)


def log_functionality(funct_dict, allele_name):
    """
    :param funct_dict: dict of strs of functionalities of each gene
    :param allele_name: str of the full gene*allele name in question
    :return: string describing the functionality of the allele in question
    """
    if allele_name not in funct_dict:
        return "Error: gene not found in 'funct_dict' dictionary! "

    elif funct_dict[allele_name] in ['F', '(F)', '[F]']:
        return "Gene is recorded as functional, and expected to be able to make working TCRs. "

    elif funct_dict[allele_name] == 'ORF':
        return "Gene is recorded as an ORF, and thus might not be able to make working TCRs. "

    elif funct_dict[allele_name] in ['P', '(P)', '[P]']:
        return "Gene is recorded as a pseudogene, and thus expected to not make working TCRs. "

    else:
        return "Error: gene has an indeterminate functionality (" + funct_dict[allele_name] + "). "


def log_partiality(part_dict, allele_name, allele_type):
    """
    :param part_dict: dict of str of partial sequence information (where available)
    :param allele_name: str of name of allele (gene*allele)
    :param allele_type: str of 'V' or 'J'
    :return: str (if gene == V which is partial in 3') or nothing
    """
    if allele_type == 'V':
        if "3'" in part_dict[allele_name]:
            return "Gene sequence is also partial in 3', so conserved C was potentially lost. "
        else:
            return ''
    else:
        return ''


def find_motif(regex_list, allele_name, type_gene, allele_seq, func_str, j_pos):
    """
    :param regex_list: list of strings describing regular expressions to search for
    :param allele_name: str of name (gene*allele) of allele in question
    :param type_gene: str specifying the gene type ('V' or 'J')
    :return: aa_seq, index of the (first detected) putative conserved motif, motif used, and the residue
    # TODO update docst
    """
    allele = allele_name.split('|')[0]

    # Vs can just be directly translated, as V-REGIONS start in-frame
    if type_gene == 'V':
        aa_seq = translate(allele_seq)
        for regex in regex_list:
            # Return the 3' most hit
            hits = [aa_seq.rfind(x) + len(x) - 1 for x in re.findall(regex, aa_seq)]
            if hits:
                return aa_seq, hits[-1], regex, aa_seq[hits[-1]]

        return [''] * 4

    # Js are trickier, requiring determination of the right frame, and conserved motifs lying in different registers
    # Start at -11/-10 positions for TRAJ/TRBJ amino acid sequences - defines the 'base' expected motif site
    elif type_gene == 'J':
        if func_str == '(F)':
            # nt_start += 1
            modifier = 1
        else:
            modifier = 0

        # Get AA seq. NB: end of the J needs to be in frame, not the start (after deleting the base that splices to C)
        modulo = (len(allele_seq) - 1 + modifier) % 3
        aa_seq = translate(allele_seq[modulo:])

        backups = [aa_seq, '', '', '']
        backup_dist = 0

        for regex in regex_list:
            # Now we want the 5' index of the hits (the F position)
            hits = [-(len(aa_seq) - aa_seq.index(x)) for x in re.findall(regex, aa_seq)]

            if len(hits) == 1:
                return aa_seq, hits[0], regex, aa_seq[hits[0]]

            elif len(hits) > 1:
                expected_hit = [x for x in hits if x == j_pos[allele]]
                if expected_hit:
                    return aa_seq, expected_hit[0], regex, aa_seq[expected_hit[0]]

                # Keep tabs of the 'least worst' alternative if no single good hit in the right place found
                for hit in hits:
                    if abs(j_pos[allele] - hit) > backup_dist:
                        backup_dist = abs(j_pos[allele] - hit)
                        backups[1:] = [hit, regex, aa_seq[hit]]

        # If we made it here we didn't get a good hit
        if backup_dist != 0:
            return backups
        else:
            return [''] * 4

    else:
        raise IOError("Unexpected gene type detected (" + type_gene + ")! ")


def find_junction_residues(all_gene_list, full_gene_dict, func_dict, partial_dict, log_text, regex_dict, j_pos):
    """
    :param all_gene_list: list of all gene*alleles input
    :param full_gene_dict: dict of full nt seqs of gens
    :param func_dict: dict of strs of recorded (predicted) functionality of each allele
    :param partial_dict: dict of strings of partiality (i.e. whether gene partial in 5'/3'/both, if applicable)
    :param log_text: str of text for log file, to be appended to
    :return: two strings to write out: that for the translate file, and the updated log file contents
    # TODO update docstr
    """

    positions = coll.defaultdict(list)
    motifs = coll.defaultdict(list)

    # If there's only one C residue in the C terminus of the V, take that as the C terminal residue.
    # NB the residue immediately prior to the conserved C is usually an L or an F...
    # ... and the residue before that is almost always a Y
    all_gene_list.sort()

    for ga in all_gene_list:
        alelle = ga.split('|')[0]
        log_text += '\n' + ga + '\t'

        gene_type = ga[3]
        if gene_type not in ['V', 'J']:
            raise IOError("Unexpected gene type (not V/J) detected: " + gene_type)

        func = log_functionality(func_dict, ga)
        part = log_partiality(partial_dict, ga, gene_type)

        # First try high-confidence motifs
        aa, position, motif, residue = find_motif(regex_dict[gene_type]['1'], ga,
                                                  gene_type, full_gene_dict[ga], func_dict[ga], j_pos)
        if motif:
            log_text += 'Found a high-confidence motif: "' + motif + '". '
            positions[ga] = position
            motifs[ga] = residue
            confidence = 'high'

        # If that fails, try less confident
        else:
            aa, position, motif, residue = find_motif(regex_dict[gene_type]['2'], ga,
                                                      gene_type, full_gene_dict[ga], func_dict[ga], j_pos)
            if motif:
                log_text += 'Found a lower-confidence motif: "' + motif + '". '
                positions[ga] = position
                motifs[ga] = residue
                confidence = 'low'
            else:
                log_text += 'Did not find a conserved residue motif - unable to translate. '
                positions[ga] = ''
                motifs[ga] = ''
                confidence = 'null'

        if positions[ga]:
            # If the detected C is far from the end, or gene is partial in 3', try to use the prototypical gene's call
            if gene_type == 'V':
                sanity_len_check = len(aa) - position > c_sanity_len

                if ((sanity_len_check or part) and confidence == 'low') or confidence == 'null':
                    proto = ga.split('*')[0] + '*01'
                    if (ga != proto) and (proto in positions):
                        positions[ga] = positions[proto]
                        motifs[ga] = motifs[proto]
                        log_text += ("Unable to find a high-confidence motif (close to the V region end) and/or gene"
                                     " is partial in it's 3 - using ") + proto + "'s. "

            # Also record whether the J motif is at the expected position relative to the 3' end of the gene
            elif gene_type == 'J':
                if position != j_pos[alelle]:
                    log_text += 'J gene motif not at the expected location (' + str(j_pos[alelle]) + '). '

        log_text += func + part

    all_vs = [x for x in all_gene_list if 'V' in x.split('*')[0]]
    all_vs.sort()

    all_js = [x for x in all_gene_list if 'J' in x.split('*')[0]]
    all_js.sort()

    all_vj = all_vs + all_js
    if len(all_vj) != len(all_gene_list):
        raise IOError('V+J gene count != all gene count! ')

    tag_text = []
    for gene in all_vj:
        tag_text.append('\t'.join([gene, str(positions[gene]), motifs[gene]]))

    return '\n'.join(tag_text), log_text


def translate_germlines(gene_id, region_type, nt_seq, j_motifs):  # TODO docstr
    # RETURNS A TRANSLATED SEQUENCE, ENSURING JS IN THE RIGHT FRAME
    if region_type == 'VARIABLE':
        return fxn.translate(nt_seq)

    elif region_type == 'JOINING':
        # Use the existing motif determination to translate the J in the correct frame
        j_match = j_motifs.loc[j_motifs['J gene'] == gene_id.split('|')[0]]

        if len(j_match) == 1:
            for f in range(3):
                j_seq = fxn.translate(nt_seq[f:])
                if len(j_seq) - j_seq.find(j_match.iloc[0]['Motif']) == abs(j_match.iloc[0]['Position']):
                    return j_seq

            # If this line is reached, the 'right' motif was
            raise IOError(f"Unable to determine correct reading frame for {gene_id}.")

        else:
            raise IOError(f"Inappropriate number of matches for J gene: ({gene_id}).")

    else:
        raise IOError(f"Unexpected gene region type (not V/J): ({region_type}).")


def get_j_motif_data(species_dir_path):
    # TODO docstr
    j_motif_df = pd.read_csv(os.path.join(species_dir_path, 'J-region-motifs.tsv'), sep='\t')
    j_pos = {}
    for row in j_motif_df.index:
        row_dat = j_motif_df.iloc[row]
        j_pos[row_dat['J gene']] = row_dat['Position']
    return j_motif_df, j_pos



region_key = {'L-': 'LEADER', 'V-': 'VARIABLE', 'J-': 'JOINING', 'EX': 'CONSTANT', 'CH': 'CONSTANT',
              'V': 'VARIABLE', 'J': 'JOINING', 'C': 'CONSTANT'}

c_sanity_len = 15

r = '\n'  # Carriage return, for log text
lb = '-----\t-----' + r  # Line break, for log text
