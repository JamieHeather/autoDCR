import argparse
import gzip
import datetime
import json
import os
import sys
import textwrap
import collections as coll
from acora import AcoraBuilder
from time import time


def format_ref_prefixes(in_loci, in_regions):
    """
    :param in_loci: list of str, detailing loci covered (from A/B/G/D)
    :param in_regions: list of str, detailing regions covered (VJ/CL)
    :return: flattened str of loci and regions
    """
    return ''.join(in_loci) + '_' + ''.join(in_regions)


def get_gene(string):
    """
    :param string: a string containing an IMGT gene name - likely a FASTA header, or something derived from it
    :return: the IMGT gene/allele
    """
    return string.split('|')[1]


def build_trie(sequence_list):
    """
    :param sequence_list: A list of DNA sequences (tags) to search for
    :return: An Aho-Corasick trie, which can be used to search sequences downstream with all of the listed tags
    """

    trie_builder = AcoraBuilder()
    for sequence in sequence_list:
        trie_builder.add(sequence)

    return trie_builder.build()


def dna_check(possible_dna):
    """
    :param possible_dna: A sequence that may or may not be a plausible DNA (translatable!) sequence
    :return: True/False
    """
    return set(possible_dna.upper()).issubset({'A', 'C', 'G', 'T', 'N'})


def opener(file_name, open_mode):
    """
    :param file_name: Path to file to open
    :param open_mode: 'open' mode parameter, e.g. 'a', 'w', 'r'
    :return: The correct opener function (open or gzip.open)
    """
    if file_name.endswith('.gz'):
        return gzip.open(file_name, open_mode + 't')
    else:
        return open(file_name, open_mode)


def import_tcr_info(species, in_loci, in_regions, mol_type, data_dir):
    """
    # TODO docstring
    """

    out_data = {'genes': {}, 'tag_genes': {}, 'tag_order': {},
                'loci': in_loci, 'regions': in_regions}

    # Check that the necessary files are present and correct
    species_dir = os.path.join(data_dir, species)
    if not os.path.exists(species_dir):
        raise IOError(f"Requested species directory ('{species.upper()}') absent from data directory:"
                      f" run 'autoDCRscripts refs' and try again.")

    file_prefix = format_ref_prefixes(in_loci, in_regions)
    if mol_type == 'aa':
        file_prefix += '_AA'

    expected_files = [file_prefix + x for x in required_ref_files]
    matching_files = [x for x in os.listdir(species_dir) if x in expected_files]

    if len(matching_files) != len(expected_files):
        raise IOError(f"Insufficient reference data for '{in_loci}' in species-specific data directory ({species_dir}).\n"
                      f"Run 'autoDCRscripts refs' and try again.")

    out_data['data_prefix'] = os.path.join(species_dir, file_prefix)
    # Get full gene sequences
    in_file_path = out_data['data_prefix'] + '.fasta'
    with open(in_file_path, 'r') as in_file:
        for fasta_header, fasta_seq, quality in readfq(in_file):
            out_data['genes'][fasta_header] = fasta_seq.upper()

    # And also the tag information...
    in_file_path = out_data['data_prefix'] + '.tags'
    with open(in_file_path, 'r') as in_file:
        count = 0
        for line in in_file:
            tag_bits = line.rstrip().split('\t')
            out_data['tag_genes'][tag_bits[0]] = tag_bits[2]
            out_data['tag_order'][tag_bits[0]] = count

            count += 1

    # ... and built into Aho-Corasick tries, which will be used to search the data
    out_data['trie'] = build_trie(out_data["tag_genes"].keys())

    # Also generate stop-region proximal tags for constant regions, to later identify full-length ORFs
    if 'C' in in_regions:
        out_data['end_tags'] = {}
        for c_gene in [x for x in out_data['genes'] if x.endswith('|C')]:

            # Presume constant regions are recorded properly, i.e. lacking the first codon's J-contributed nt
            translation = translate(out_data['genes'][c_gene][2:])
            first_stop = translation.find('*')
            if first_stop == -1:
                out_data['end_tags'][c_gene] = out_data['genes'][c_gene][-20:]
            else:
                stop_nt = (first_stop * 3) + 2
                out_data['end_tags'][c_gene] = out_data['genes'][c_gene][stop_nt-17:stop_nt+3]

    return out_data


def import_translate_info(ref_data):
    """
    # TODO fix docstring
    Establishes global dictionaries which contain all the information required for translating TCRs:
    - trans_pos : Position of conserved CDR3 junction defining residues { gene_name: index }
    - trans_res : Identity of conserved CDR3 junction defining residues { gene_name: residue_character }
    :param input_arguments: CLI input arguments
    :return: Nothing, all dictionaries are set up as global objects
    """

    ref_data['trans_pos'] = {}
    ref_data['trans_res'] = {}
    if not os.path.exists(ref_data['data_prefix'] + '.translate'):
        raise IOError(f"Missing .translate file in species-specific data directory "
                      f"({ref_data['data_prefix'] + '.translate'}).\n"
                      f"Run 'autoDCRscripts refs' and try again.")

    with open(ref_data['data_prefix'] + '.translate', 'r') as in_file:
        for line in in_file:
            bits = line.rstrip().split('\t')
            ref_data['trans_pos'][bits[0]] = int(bits[1])
            ref_data['trans_res'][bits[0]] = bits[2]

    return ref_data

def readfq(fastx_file):
    """
    readfq(file):Heng Li's Python implementation of his readfq function
    https://github.com/lh3/readfq/blob/master/readfq.py
    :param fastx_file: opened file containing fastq or fasta reads
    :yield: read id, read sequence, and (where available) read quality scores
    """

    last = None  # this is a buffer keeping the last unprocessed line
    while True:
        if not last:  # the first record or a record following a fastq
            for l in fastx_file:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break

        if not last:
            break

        name, seqs, last = last[1:], [], None  # This version takes the whole line (post '>')
        for l in fastx_file:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])

        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break

        else:  # this is a fastq record
            sequence, leng, seqs = ''.join(seqs), 0, []
            for l in fastx_file:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(sequence):  # have read enough quality
                    last = None
                    yield name, sequence, ''.join(seqs)  # yield a fastq record
                    break

            if last:  # reach EOF before reading enough quality
                yield name, sequence, None  # yield a fasta record instead
                break


def rev_comp(seq):
    """
    :param seq: Any DNA string, composed of just the four bases (upper or lower case)
    :return: The reverse complement of that DNA string (maintaining same case)

    """
    return seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]


def sort_read_bits(whole_read, whole_qual, barcode_details):
    """
    # TODO fix docstring
    :param whole_read: The entire read, from start to finish
    :param whole_qual: The quality scores of that entire reads
    :param input_arguments: The argparse input arguments, to determine whether/how to chop up the read
    :return: 4 str: the TCR-containing portion of the read, it's quality, and the barcode-containing read/its quality
    """

    if barcode_details:
        # TODO all this needs updating
        pass
        # tcr_read = whole_read[input_arguments['bclength']:].upper()
        # tcr_qual = whole_qual[input_arguments['bclength']:]
        # bc_read = whole_read[:input_arguments['bclength']].upper()
        # bc_qual = whole_qual[:input_arguments['bclength']]
        # return tcr_read, tcr_qual, bc_read, bc_qual
    else:
        return whole_read.upper(), whole_qual, '', ''

def tidy_output(output_dict, header_list):
    # TODO docstr

    # Remove the '|X' disambiguation strings from the ends of gene calls
    for field in [x for x in output_dict if x.endswith('_call')]:
        if output_dict[field]:
            output_dict[field] = ','.join([x.split('|')[0] for x in output_dict[field].split(',')])

    # And ensure all necessary fields for plotting are included
    for field in [x for x in header_list if x not in list(output_dict.keys())]:
        output_dict[field] = ''

    return output_dict


def today(date_format=''):
    """
    :return: Today's date, in one of two formats, defaulting to ISO
    """
    date_today = datetime.datetime.today().date().isoformat()
    if date_format == 'ncbi':
        date_today = datetime.datetime.strptime(date_today, '%Y-%m-%d').strftime('%d-%b-%Y').upper()
    return date_today


def pad_spaces(str_to_pad, len_to_pad):
    """
    :param str_to_pad: str to pad with spaces
    :param len_to_pad: int of length of padding required
    :return: str_to_pad prefixed with len_to_pad number of spaces (for format matching)
    """
    return str_to_pad + ' ' * (len_to_pad - len(str_to_pad))


def output_genbank(sequence_name, full_sequence, description, topology, features, save_dir_path,
                   species, numbering, plot_multi=True, division='SYN',
                   accession='.', version='.', keywords='.', source='.',
                   authors='.', title='Direct Submission', journal='.'):
    """
    Save an annotated GenBank (.gb) file for easy visualisation of DNA features (e.g. stitched TCRs)
    :param sequence_name: str of sequence name, for saving
    :param full_sequence: str of full DNA sequence
    :param description: str of description to include in output file
    :param topology: str of topology ('linear' or 'circular')
    :param features: lists of tuples [(str, str, [list of str], [list of int]], detailing target output features in GB:
        1) name, 2) type (if known), 3) DNA sequence(s) (allowing multi-part features), 4) index location(s) (optional)
    :param save_dir_path: str of path to folder in which genbank files should be saved
    :param species: str of species to go in 'organism' field
    :param numbering: boolean, whether to incrementally number regions of the same type
    :param plot_multi: boolean, whether to plot all instances of multi-mapping sequences given without indexes (def=0)
    :param division: GenBank 'division', as per https://www.ncbi.nlm.nih.gov/genbank/samplerecord/#GenBankDivisionB
    :param accession: optional str field (can be left as default '.')
    :param version: optional str field (can be left as default '.')
    :param keywords: optional str field (can be left as default '.')
    :param source: optional str field (can be left as default '.')
    :param authors: optional str field (can be left as default '.')
    :param title: optional str field (can be left as default 'Direct submission')
    :param journal: optional str field (can be left as default '.')
    :return: [nothing] just saves the appropriate file
    """
    pad_len1 = 12
    pad_len2 = 21

    seq_len = str(len(full_sequence))
    out_str = [pad_spaces('LOCUS', pad_len1) + sequence_name + '         ' + seq_len +
               ' bp    DNA    ' + topology + '    ' + division + '    ' + today('ncbi'),
               pad_spaces('DEFINITION', pad_len1) +
               description,

               pad_spaces('ACCESSION', pad_len1) + accession,
               pad_spaces('VERSION', pad_len1) + version,
               pad_spaces('KEYWORDS', pad_len1) + keywords,
               pad_spaces('SOURCE', pad_len1) + source,
               pad_spaces('  ORGANISM', pad_len1) + species,
               pad_spaces('REFERENCE', pad_len1) + '1 (bases 1 to ' + seq_len + ')',
               pad_spaces('  AUTHORS', pad_len1) + authors,
               pad_spaces('  TITLE', pad_len1) + title,
               pad_spaces('  JOURNAL', pad_len1) + journal,

               pad_spaces('FEATURES', pad_len2) + 'Location/Qualifiers',
               pad_spaces('     source', pad_len2) + '1..' + str(len(full_sequence)),
               pad_spaces('', pad_len2) + '/organism="' + species + '"']

    # Plot the specified input features
    if features:
        upper_seq = full_sequence.upper()
        for feature in features:
            index_count = 1  # Used to number features if numbering

            # Determine feature type
            if feature[1]:
                f_type = feature[1]
                if not f_type.endswith('_'):
                    f_type += '_'
            else:
                f_type = 'misc_feature'

            # Determine feature location
            num_seqs = len(feature[2])
            num_is = len(feature[3])
            if (num_seqs > 0 and num_is > 0) and num_seqs != num_is:
                raise IOError("Provided feature entry has an inappropriate number of index locations detected: "
                              + str(feature))

            for subfeatx in range(num_seqs):
                found = False
                f_name = feature[0].rstrip()
                if numbering:
                    f_name += ' #' + "{:02d}".format(index_count)

                # If index information is provided, use it
                if num_is > 0:
                    index_seq = upper_seq[feature[3][subfeatx]:feature[3][subfeatx] + len(feature[2][subfeatx])]
                    if index_seq == feature[2][subfeatx].upper():
                        out_str.append(genbank_write_feature(f_name, feature[3][subfeatx],
                                                                 len(feature[2][subfeatx]), f_type))
                        found = True

                # Otherwise look for sequence matches
                if not found:
                    hits = [x for x in find_all_substr(feature[2][subfeatx].upper(), upper_seq)]
                    if len(hits) == 1:

                        out_str.append(genbank_write_feature(f_name, hits[0],
                                                             len(feature[2][subfeatx]), f_type))
                    elif plot_multi:
                        for hit in hits:
                            out_str.append(genbank_write_feature(f_name + ' (multiple matches)', hit,
                                                                 len(feature[2][subfeatx]), f_type))

                index_count += 1

    out_str.append('ORIGIN')
    out_str.append(ncbi_format_dna(full_sequence, pad_len1))
    out_str.append('//')

    with open(os.path.join(save_dir_path, sequence_name) + '.gb', 'w') as out_file:
        out_file.write('\n'.join(out_str))


def find_all_substr(substr, full_seq):
    """
    :param substr: str, substring to look for in ...
    :param full_seq: str, longer string in which to look for instances of the provided substring
    :yield: indexes of all matches
    NB: case sensitive
    """
    match = full_seq.find(substr)
    while match != -1:
        yield match
        match = full_seq.find(substr, match + 1)


def genbank_write_feature(feat_name, feat_loc, feat_len, feat_type='misc_feature', space_len=21):
    """
    :param feat_name: str, complete name of DNA feature (including any numbering / qualifiers etc)
    :param feat_loc: str, describing 1-indexed position range of the feature in format 'X..Y'
    :param feat_len: int, length of the DNA feature
    :param feat_type: str, type of DNA feature if known
    :param space_len: int, length of spacer used to left indent feature text in GenBank file
    :return: str of formatted GenBank entry for a given feature
    """

    feat_str = (pad_spaces('     ' + feat_type, space_len)
                + str(feat_loc + 1) + '..' + str(feat_loc + feat_len) + '\n' +
                pad_spaces('', space_len) + '/label="' + feat_name + '"'
                )
    return feat_str


def ncbi_format_dna(in_str, pad_len):
    """
    :param in_str: str of DNA sequence
    :param pad_len: int of padding required
    :return: str of DNA sequence in NCBI format (indented/space-separated 10-mers)
    """
    dna_text = ''
    dna_chunks = [in_str[x:x+10] for x in range(0, len(in_str), 10)]
    line_chunks = [dna_chunks[x:x+6] for x in range(0, int(len(in_str)/6), 6)]
    pos = 1
    for line in line_chunks:
        if line:
            dna_text += ' '.join([num_pad(pos, pad_len)] + line) + '\n'
            pos += 60
        else:
            break

    return dna_text


def num_pad(number, len_to_pad):
    """
    :param number: int describing a number to be converted to str and prefixed with spaces
    :param len_to_pad: int of number of spaces to pad
    :return: str of 'number' prefixed with len_to_pad number of spaces (for format matching)
    """
    str_num = str(number)
    return ' ' * (len_to_pad - len(str_num)) + str_num


def find_cdr3(tcr, ref_data, parameters):
    """
    :param tcr: dict of TCR featured detected in a given read
    :param ref_data: dict of reference data
    :param parameters: dict of autoDCRscripts-related parameters
    :return: Same tcr dict, with additional corresponding translated/CDR3 properties (if found)
    """

    # Need to check for presence of both jump values
    # Can lack one whilst still having a detected rearrangement if say multiple Vs/Js in one read

    if 'v_jump' in tcr and 'j_jump' in tcr:

        # If multiple genes detected, just arbitrarily pick one to use for translation parameters
        if ',' in tcr['v_call']:
            translate_v = eg_gene(tcr['v_call'])
        else:
            translate_v = tcr['v_call']

        if ',' in tcr['j_call']:
            translate_j = eg_gene(tcr['j_call'])
        else:
            translate_j = tcr['j_call']

        full_inferred = ref_data['genes'][translate_v][:tcr['v_jump']] \
                                  + tcr['inter_tag_seq'] + ref_data['genes'][translate_j][tcr['j_jump']:]

        if parameters['mol_type'] == 'nt':
            tcr['inferred_full_nt'] = full_inferred
            tcr['inferred_full_aa'] = translate(tcr['inferred_full_nt'])

        elif parameters['mol_type'] == 'aa':
            tcr['inferred_full_aa'] = full_inferred

        tcr['junction_aa'] = tcr['inferred_full_aa'][ref_data['trans_pos'][translate_v]:
                                                     ref_data['trans_pos'][translate_j] + 1]

        # Assume productivity, and remove it as applicable, checking for the various required parameters
        # Note that it's possible for multiple different reasons to be non-productive to be true, in different combos
        tcr['productive'] = 'T'
        tcr['vj_in_frame'] = 'T'

        # Check whether it's feasibly in frame
        if (len(tcr['inferred_full_nt']) - 1) % 3 != 0 and parameters['mol_type'] == 'nt':
            tcr['productive'] = 'F'
            tcr['vj_in_frame'] = 'F'

        # Check for stop codons...
        if '*' in tcr['inferred_full_aa']:
            tcr['stop_codon'] = 'T'
            tcr['productive'] = 'F'
        else:
            tcr['stop_codon'] = 'F'

        # Need to account for cases where there is no detectable/valid CDR3
        if tcr['junction_aa']:

            # ... and the conserved V gene residue at the right position...
            if tcr['junction_aa'][0] != ref_data['trans_res'][translate_v]:
                tcr['conserved_c'] = 'F'
            else:
                tcr['conserved_c'] = 'T'

            # ... and same for the J gene...
            if tcr['junction_aa'][-1] != ref_data['trans_res'][translate_j]:
                tcr['conserved_f'] = 'F'
                tcr['productive'] = 'F'
            else:
                tcr['conserved_f'] = 'T'

            # And check whether the CDR3 falls within the expected length range
            if parameters['cdr3_limit'] > 0:
                if len(tcr['junction_aa']) <= parameters['cdr3_limit']:
                    tcr['cdr3_in_limit'] = 'T'
                else:
                    tcr['cdr3_in_limit'] = 'F'
                    tcr['junction_aa'] = ''
                    tcr['productive'] = 'F'

            else:
                tcr['cdr3_in_limit'] = ''

        else:
            tcr['productive'] = 'F'

        # Cryptic splices in leader processing can result in a variety of translation issues, let's catch those
        if tcr['inferred_full_aa'] and not tcr['junction_aa'] and tcr['inter_tag_seq']:
            tcr = attempt_salvage_irregular_cdr3s(tcr, ref_data, translate_v, translate_j)

        if tcr['productive'] == 'F':
            # TODO fix - this will catch odd non-C CDR3s

            if tcr['junction_aa']:

                # Additional check to try to salvage recovery of CDR3s in frame-shifted receptors
                # Works off scanning from conserved F (if present) and scanning N-wards for the next C
                if tcr['stop_codon'] == 'T' and tcr['conserved_f'] == 'F':
                    tcr = attempt_salvage_irregular_cdr3s(tcr, ref_data, translate_v, translate_j)

                if 'non_productive_junction_aa' not in tcr:
                    tcr['non_productive_junction_aa'] = tcr['junction_aa']
                tcr['junction_aa'] = ''

            # TODO else: somehow note that most failures to generate CDR3s are from V region (in)dels?
        else:
            tcr['junction'] = tcr['inferred_full_nt'][
                              ref_data['trans_pos'][translate_v] * 3:
                              (ref_data['trans_pos'][translate_j] * 3) + 2]

    return tcr


def attempt_salvage_irregular_cdr3s(tcr_dict, ref_data, v_translate, j_translate):
    """
    # TODO fix docstr
    Try to find CDR3s in irregular rearrangements (particularly necessary when using irregular V-REGION references)
    :param tcr_dict: Dictionary describing rearrangement under consideration
    :param v_translate: str ID of V gene being used for translation purposes (maybe differ if >1 detected)
    :param j_translate: str ID of V gene being used for translation purposes (maybe differ if >1 detected)
    :return: tcr_dict, hopefully updated with any CDR3s the regular assumptions failed to catch now annotated
    """
    potential_frames = []
    potential_translations = []

    for frame in range(3):
        potential_translation = translate(tcr_dict['inferred_full_nt'][frame:])

        if potential_translation[ref_data['trans_pos'][j_translate]] == ref_data['trans_res'][j_translate]:
            potential_frames.append(frame)
            potential_translations.append(potential_translation)

    # If there's two frames, one which has stops and one which doesn't, keep the one without
    if len(potential_frames) == 2:
        discard = []
        for i in [0, 1]:
            if '*' in potential_translations[i]:
                discard.append(i)
        if len(discard) == 1:
            potential_frames.pop(discard[0])

    # If there's a single frame in which there's a conserved F residue at the appropriate location...
    if len(potential_frames) == 1:
        # ... look upstream until it finds a C (within a typical distance)
        potential_translation = translate(tcr_dict['inferred_full_nt'][potential_frames[0]:])
        j_f_pos = ref_data['trans_pos'][j_translate]
        for i in range(8, 21):
            potential_junction_aa = potential_translation[j_f_pos - i:j_f_pos + 1]
            if potential_junction_aa[0] == ref_data['trans_res'][v_translate]:
                tcr_dict['non_productive_junction_aa'] = potential_junction_aa
                tcr_dict['inferred_full_aa'] += ',' + translate(
                    tcr_dict['inferred_full_nt'][potential_frames[0]:])
                tcr_dict['conserved_f'] = 'T'
                break

    return tcr_dict


def get_deletions(results_dict, ref_data, parameters):
    """
    # TODO fix docstring
    :param results_dict: The dictionary containing all of the TCR details discovered so far
    :param input_arguments: The argparse input arguments
    :return: the same results_dict, with additional information relating to the deletions discovered in the germline V/J
    """
    # TODO filter out cases where end of gene is beyond the limits of the read

    germlines = coll.defaultdict()
    for gene in ['V', 'J']:

        # Take the equivalent un-deleted sequence of that germline gene to compare against the recombined sequence
        full_germline = ref_data['genes'][eg_gene(results_dict[gene.lower() + '_call'])]
        if gene == 'V':
            try:
                tag_start = full_germline.index(results_dict['inter_tag_seq'][:results_dict['v_tag_len']])
                germlines[gene] = full_germline[tag_start:]
                recombined = results_dict['inter_tag_seq'][:len(germlines[gene])]
            except Exception:
                results_dict[gene.lower() + '_deletion_found'] = False
                break

        elif gene == 'J':
            try:
                tag_site = full_germline.index(results_dict['inter_tag_seq'][-results_dict['j_tag_len']:])
                germlines[gene] = full_germline[:tag_site + results_dict['j_tag_len']]
                recombined = results_dict['inter_tag_seq'][-len(germlines[gene]):]
            except Exception:
                results_dict[gene.lower() + '_deletion_found'] = False
                break

        # Need to start at the right place and take sliding windows in the correction direction
        # V genes = Start at 3' and move 5' / J genes = start at 5' and move 3'
        position = {'V': -10, 'J': 0}
        increment = {'V': -1, 'J': 1}
        matched = False
        deletions = 0

        # Starting at the end of the gene, move the sliding window 1 nt away each iteration
        while not matched and deletions < len(recombined):
            if recombined[position[gene]:][:10] == germlines[gene][position[gene]:][:10]:
                matched = True
            else:
                position[gene] += increment[gene]
                deletions += 1

        if matched or deletions < parameters['deletion_limit']:
            results_dict[gene.lower() + '_deletions'] = deletions
            results_dict[gene.lower() + '_deletion_found'] = True
        else:
            results_dict[gene.lower() + '_deletion_found'] = False

    # Use the values determined above to pull out the 'insert' sequence, i.e. the non-template nt between ends of V/J
    if results_dict['v_deletion_found'] and results_dict['j_deletion_found']:
        its = results_dict['inter_tag_seq']
        try:
            # Use the length of the V (minus deletions) up to the 10-mer of post-deleted J sequence found above
            results_dict['insertion'] = its[len(germlines['V']) - results_dict['v_deletions']:its.index(
                germlines['J'][results_dict['j_deletions']:results_dict['j_deletions'] + 10])]
        except Exception:
            results_dict[gene.lower() + '_deletion_found'] = False

    return results_dict


def populate_reads(sequence_read, sequence_qual, requested_orientation):
    # TODO docstr
    out_reads = []
    request = requested_orientation.lower()
    if request in ['both', 'either', 'b']:
        return [{'read': sequence_read, 'quality': sequence_qual, 'rev_comp': 'F'},
                {'read': rev_comp(sequence_read), 'quality': sequence_qual[::-1], 'rev_comp': 'T'}]
    elif request in ['forward', 'fwd', 'f']:
        return [{'read': sequence_read, 'quality': sequence_qual, 'rev_comp': 'F'}]
    elif request in ['reverse', 'rev', 'r']:
        return [{'read': rev_comp(sequence_read), 'quality': sequence_qual[::-1], 'rev_comp': 'T'}]
    else:
        raise IOError(f"Incorrect orientation argument provided: '{requested_orientation}'\n"
                      f"Please specify one of the three options: forward/reverse/both - or f/r/b.")



def dcr(read, quality, ref_data, parameters, tsv_headers):
    # TODO docstr

    # Depending on the mode, perform different combinations of TCR processing tasks
    # The shared functionality is the application of an Aho Corasick trie to search for TCR sub-sequences
    # This is to be applied to the requested read orientations
    strands = populate_reads(read, quality, parameters['orientation'])

    for strand in strands:
        tag_hits = find_tag_hits(strand['read'], ref_data)
        # print(1, tag_hits)  # TODO rm
        if not tag_hits or tag_hits is None:
            continue

        # If suitable to the mode, look for rearrangements (i.e. V/J), including the junction-distal tags
        if parameters['mode'] in ['VJCDR3', 'FULL'] and ref_data['regions'] == 'JV':
            # TODO add the relevant modes
            # TODO also need to make sure the tag hit dicts are turned once successful
            tag_hits = call_rearrangements(strand, tag_hits, ref_data)
            if not tag_hits:
                continue

            # Pad dict with empty values for all columns required in final output document (AIRR community format)
            for field in tsv_headers:
                if field not in tag_hits.keys():
                    tag_hits[field] = ''

            tag_hits = find_cdr3(tag_hits, ref_data, parameters)
            # TODO make modes/options to NOT do deletions but still do recombs?
            tag_hits = get_deletions(tag_hits, ref_data, parameters)

            # If this strand has reached this far, it must be suitably rearranged, so it can be output
            break

        elif parameters['mode'] == 'DISCOVER':# and ref_data['regions'] == 'JV':
            # # TODO just focusing on Vs for now, ensure have full Vs
            # if ref_data['genes'][tag_hits['v_call'].split(',')[0]].startswith(tag_hits['v_matching_tags'][0][0]):#
            # tag_hits = discover_genes(tag_hits, read, ref_data)
            tag_hits['sequence'] = read
            return tag_hits

    if not tag_hits or tag_hits is None:
        return

    # Make sure there are 'x_call' fields present, even if not detected
    if ref_data['regions'] == 'JV':
        fields = ['v', 'j']
    elif ref_data['regions'] == 'CL':
        fields = ['l', 'v', 'j', 'c']
    else:
        print(f"Warning: unexpected region combination ({ref_data['regions']}) detected, "
              f"cannot auto-complete results dict. ")

    for f in fields:
        if f + '_call' not in tag_hits.keys():
            tag_hits[f + '_call'] = ''

    return tag_hits


def get_3p_seq_offset(full_seq, ref_dat, gene, region, last_tag):
    # TODO docstr
    # Find out how much farther 3' after the terminal tag match a region continues / whether it goes further

    tag_len = len(last_tag[0])
    full_ref_seq = ref_dat['genes'][gene + '|' + region.upper()]
    end_ref_seq = full_ref_seq[full_ref_seq.index(last_tag[0]) + tag_len:]
    end_seq = full_seq[last_tag[1] + tag_len:]

    shared_prefix = os.path.commonprefix([end_seq, end_ref_seq])

    return shared_prefix, end_seq[len(shared_prefix):]


def eg_gene(gene_str):
    # TODO docstr
    # Pick the first example gene in an ambiguous list for various purposes
    return gene_str.split(',')[0]



def determine_gene_boundaries(tcr_dict, vj_refdat, lc_refdat):
    # todo docstr

    ref_dict = {'v': vj_refdat, 'j': vj_refdat,
                'l': lc_refdat, 'c': lc_refdat}
    annotations = {}

    for r in ['v', 'j', 'l', 'c']:
        if tcr_dict[r + '_call']:

            # Need to record (for ambiguous calls) which gene sequence is being used
            annotations[r] = eg_gene(tcr_dict[r + '_call'])

            # Need to determine the edges of the regions
            start = tcr_dict[r + '_matching_tags'][0][1]
            end = (tcr_dict[r + '_matching_tags'][-1][1] +
                   len(tcr_dict[r + '_matching_tags'][-1][0]))

            # get_3p_seq_offset(full_seq, ref_dat, gene, region, last_tag)  # TODO here

            # Fill the sequence between the last tag and the end of the region (tags being tiled from 5'-3')
            r_3p, rest_3p = get_3p_seq_offset(tcr_dict['sequence'], ref_dict[r], annotations[r],
                                              r, tcr_dict[r + '_matching_tags'][-1])
            end += len(r_3p)

            # Add those features to the dict for output
            tcr_dict[r + '_start'] = start
            tcr_dict[r + '_end'] = end
            tcr_dict[r + '_seq'] = tcr_dict['sequence'][start:end]

            # Determine the label
            label = tcr_dict[r + '_call']
            if r == 'l':
                # Just label the leader region in general if it's the same as the V (i.e. expected)
                if (tcr_dict['v_call'] == tcr_dict['l_call'] or
                        tcr_dict['v_call'] in tcr_dict['l_call'] or
                        tcr_dict['l_call'] in tcr_dict['v_call']):
                    label = 'leader'
                # Otherwise label with the detected gene(s)
                else:
                    label = tidy_gene_list(tcr_dict['l_call'], resolution='gene') + ' (L)'

            if r == 'c':
                # Often won't have sufficient read length to distinguish different TRxC alleles
                if ',' in tcr_dict['c_call']:
                    label = tidy_gene_list(tcr_dict['c_call'], resolution='gene')

            # 5' partiality checks
            if r in ['v', 'l']:
                if r == 'v' and tcr_dict['v_jump'] != 0 and 1 == 2:
                    label += ' (partial in 5\')'
                else:
                    if not ref_dict[r]['genes'][annotations[r] + '|' + r.upper()].startswith(
                            tcr_dict[r + '_matching_tags'][0][0]):
                        label += ' (partial in 5\')'

            # 3' partiality checks
            if r in ['j', 'c']:
                if not ref_dict[r]['genes'][annotations[r] + '|' + r.upper()].endswith(
                        tcr_dict[r + '_seq'][-20:]):
                    label += ' (partial in 3\')'
            # TODO uncomment and make use of v/j jumps

            tcr_dict[r + '_label'] = label

    return tcr_dict


def discover_genes(read, read_qual, vj_data, lc_data, dcr_parameters, headers):
    # TODO docstring

    tag_search = dcr(read, read_qual, vj_data, dcr_parameters, headers)
    if not tag_search:
        return

    if 'v_call' not in tag_search.keys() or 'j_call' not in tag_search.keys():
        return

    tag_search = find_full_feats(tag_search, lc_data, dcr_parameters)
    # Filter down to those reads which contain whole V/J regions
    if 'v_matching_tags' not in tag_search.keys():
        return

    print([x[2].split('*')[0] for x in tag_search['v_matching_tags']])
    print(tag_search['v_call'])
    #
    if (len(list(set([x[2].split('*')[0] for x in tag_search['v_matching_tags']]))) == 1 or
        'l_call' in tag_search) and 'c_call' in tag_search:

        # TODO fix - this is too stringent, can still look for novel Vs in reads that have only partial Js
        for r in ['v', 'j']:
            if r + '_matching_tags' in tag_search.keys():
                # Only run on cases with unambiguous gene calls (even if allele is ambiguous, which it may well be)
                if len(list(set([x[2].split('*')[0] for x in tag_search[r+'_matching_tags']]))) == 1:
                    tmp = noncontiguous_tag_check(read, tag_search[r + '_matching_tags'],
                                                  10, tag_search[r + '_call'], vj_data)
                    if tmp:
                        tag_search[r + '_mismatch'] = tmp
                        tag_search[r + '_call'] = tag_search[r + '_call'].replace('|', '_'+tag_search[r+'_mismatch']['snp']+'|')

    return tag_search


def variant_position(start_index, ref_seq, var_seq):
    # todo docstr
    if len(ref_seq) != len(var_seq):
        return '!'
    var_site = [x for x in range(len(ref_seq)) if ref_seq[x] != var_seq[x]]
    if len(var_site) == 1:
        var_site = var_site[0]
        return ref_seq[var_site] + str(start_index + 1 + var_site) + var_seq[var_site]
    else:
        return '?'


def noncontiguous_tag_check(sequence, tag_hits, win_len, gene, ref_dat):
    """
    Used in trying to identify potential novel alleles, based on a break in consecutive tag hits
    :param sequence: read being decombined
    :param tag_hits: the list of V/J tag hits produced by find_tag_hits
    :param win_len: length of expected window between adjacent tags
    :param gene: gene used
    :return: dict detailing the detected noncontiguous mismatch, i.e. the tags before/after and sequence in between
    """

    # Then simply find 'missing' tags, based on what tag positions we'd expect to see, given a certain tag window length
    tag_positions = [x[1] for x in tag_hits]
    non_contig = [x for x in range(min(tag_positions), max(tag_positions), win_len) if x not in tag_positions]

    # TODO add check here that all tags map to the same gene?
    if not non_contig:
        return
    else:
        # Find the substring that's not covered by the tags, throwing it out if the tags before/after aren't in sync
        # Find the tags and positions in the read
        tag1 = [x for x in tag_hits if x[1] == non_contig[0] - win_len]
        if tag1:
            tag1 = tag1[0]
            tag2 = [x for x in tag_hits if x[1] == non_contig[-1] + win_len]
            if tag2:
                tag2 = tag2[0]
            else:
                return
        else:
            return

        # And the corresponding parts in the germline gene
        tag1_gl_index = ref_dat['genes'][eg_gene(gene)].index(tag1[0])
        tag2_gl_index = ref_dat['genes'][eg_gene(gene)].index(tag2[0])

        # Otherwise report back all the relevant bracketing and intervening sequence/tag information
        mismatch_dict = {
            'tag1_index': tag1[1],
            'tag2_index': tag2[1],
            'tag1_seq': tag1[0],
            'tag2_seq': tag2[0],
            'tag1_gl_index': tag1_gl_index,
            'tag2_gl_index': tag2_gl_index,
            'intervening': sequence[tag1[1] + len(tag1[0]):tag2[1]],
            'gl_equivalent': ref_dat['genes'][eg_gene(gene)][tag1_gl_index+len(tag1[0]):tag2_gl_index]
        }

        mismatch_dict['snp'] = variant_position(tag1_gl_index+len(tag1[0]),
                                                mismatch_dict['intervening'], mismatch_dict['gl_equivalent'])
        return mismatch_dict


def find_full_feats(tcr_dict, ref_dat, parameters):
    # Todo doctsr
    feat_check = find_tag_hits(tcr_dict['sequence'], ref_dat)
    if feat_check:

        # Account for many V alleles sharing the same leader sequences
        # I.e. presume that the same L was used if the V is unambiguously determined
        if 'l_call' in feat_check:
            if 'v_call' not in tcr_dict:
                # TODO fix
                print(tcr_dict)
                sys.exit()
            if ',' in feat_check['l_call'] and ',' not in tcr_dict['v_call']:
                if tcr_dict['v_call'].split('|')[0] in [x.split('|')[0] for x in feat_check['l_call'].split(',')]:
                    feat_check['l_call'] = tcr_dict['v_call'].replace('|V', '|L')

            # Look for the beginning of the leader sequence
            ref_l = eg_gene(feat_check['l_call'])
            if ref_dat['genes'][ref_l].startswith(feat_check['l_matching_tags'][0][0]):
                feat_check['start_pos'] = feat_check['l_matching_tags'][0][1]
                feat_check['found_start'] = 'T'

        # Similarly can look for the start and end of the constant region
        if 'c_call' in feat_check:
            ref_c = eg_gene(feat_check['c_call'])
            if ref_dat['end_tags'][ref_c] in tcr_dict['sequence']:
                feat_check['stop_pos'] = (tcr_dict['sequence'].index(ref_dat['end_tags'][ref_c])
                                          + len(ref_dat['end_tags'][ref_c]))
                feat_check['found_stop'] = 'T'

        if parameters['mode'] == 'DISCOVER':
            return tcr_dict

        # If both are found, see if they're in frame
        if 'start_pos' in feat_check and 'c_call' in feat_check:
            if ref_dat['genes'][ref_c].startswith(feat_check['c_matching_tags'][0][0]):
                l_to_c = feat_check['c_matching_tags'][0][1] - feat_check['start_pos'] +2
                if l_to_c % 3 == 0:
                    feat_check['lc_in_frame'] = 'T'
                else:
                    feat_check['lc_in_frame'] = 'F'

        if 'found_start' in feat_check and 'found_stop' in feat_check:
            feat_check['orf_len'] = feat_check['stop_pos'] - feat_check['start_pos']
            if feat_check['orf_len'] % 3 == 0:
                feat_check['orf_in_frame'] = 'T'
            else:
                feat_check['orf_in_frame'] = 'F'

        # Add these to the tcr dictionary
        # for field in full_feat_headers:  # TODO remove?
        #     if field in feat_check:
        #         tcr_dict[field] = feat_check[field]
        for field in feat_check:
            tcr_dict[field] = feat_check[field]

    return tcr_dict


def call_rearrangements(strand, trie_results, ref_data):
    if 'v_call' in trie_results and 'j_call' in trie_results:

        trie_results['rev_comp'] = strand['rev_comp']
        trie_results['sequence'] = strand['read']
        for gene_type in ['V', 'J']:
            # Then go through the matched tag and find the relevant hit that uses the tag combination found
            if gene_type in ['V', 'L']:
                tag_order = [x for x in trie_results[gene_type.lower() + '_matching_tags']]
            elif gene_type in ['J', 'C']:
                tag_order = [x for x in trie_results[gene_type.lower() + '_matching_tags']][::-1]

            else:
                raise ValueError(f"Unexpected gene type during 'dcr()' execution: {gene_type}")  # TODO reword better

            # Some genes are uniquely identifiable only through their unique appearance in the intersection of many tags
            # Need to find the most-distal tag that covers that gene/all of those genes
            # Use reference positions from the only OR first listed gene
            call_bits = trie_results[gene_type.lower() + '_call'].split(',')

            for match in tag_order:
                tag_genes_sep = ref_data['tag_genes'][match[0]].split(',')

                if all(g in tag_genes_sep for g in call_bits):

                    trie_results[gene_type.lower() + '_jump'] = ref_data['genes'][call_bits[0]].index(match[0])

                    # Need to add on the J tag length to make the inter-tag sequence include the J tag
                    if gene_type == 'J':
                        trie_results[gene_type.lower() + '_jump'] += len(match[0])

                    break  # Only need to find the outermost jump values

                else:
                    # Of the V|J tags found, there are no alleles represented among all and thus no unambig call
                    pass
                    # TODO if desired not-rearranged-but-TCR-containing reads could be output here

            trie_results[gene_type.lower() + '_tag_seq'] = match[0]
            trie_results[gene_type.lower() + '_tag_len'] = len(match[0])
            trie_results[gene_type.lower() + '_tag_position'] = match[1]

        trie_results['inter_tag_seq'] = strand['read'][trie_results['v_tag_position']:
                                                       trie_results['j_tag_position'] +
                                                       trie_results['j_tag_len']]
        trie_results['inter_tag_qual'] = strand['quality'][trie_results['v_tag_position']:
                                                           trie_results['j_tag_position'] +
                                                           trie_results['j_tag_len']]

        return trie_results

    else:
        return


def find_tag_hits(sequence, ref_data):
    """
    # TODO fix docstring
    Does the actual tag searching
    :param sequence: DNA seq to search for tag hits (i.e. sequences matching tags corresponding to TCR genes)
    :return: defaultdict containing fields corresponding to the most distal/most unique V/J tag hits
    """

    # Search the read using the relevant Aho-Corasick trie, and pull out the genes those tags come from
    # call, tag_index, match, jump, i, len_check = '', '', ['', ''], '', '', ''
    out_dict = coll.defaultdict()
    check = ref_data['trie'].findall(sequence)
    check = [(x[0], x[1], ref_data['tag_genes'][x[0]], ref_data['tag_genes'][x[0]][-1]) for x in check]

    found = 0
    for gene_type in ref_data['regions']:
        specific_check = [x for x in check if x[3] == gene_type]

        if specific_check:
            found += 1

            # Go through and find the minimal set of genes that are featured in the maximum number of tags
            gene_lists = [x[2] for x in specific_check]

            # Go through and find the minimal set of genes that are featured in the maximum number of tags
            gene_counts = coll.Counter()
            for gene_list in [x.split(',') for x in list(gene_lists)]:
                for g in gene_list:
                    gene_counts[g] += 1

            for i in range(len(specific_check), 0, -1):
                hits = [x for x in gene_counts if gene_counts[x] == i]
                if len(hits) > 0:
                    hits.sort()
                    out_dict[gene_type.lower() + '_call'] = ','.join(hits)
                    break

            out_dict[gene_type.lower() + '_number_tags_matches'] = i
            out_dict[gene_type.lower() + '_number_tags_total'] = len(specific_check)
            out_dict[gene_type.lower() + '_matching_tags'] = specific_check

    if len(out_dict) > 1:
        return out_dict
    else:
        return


def get_output_name(tcr_name, tcr_dict):
    # TODO docstr
    # Preferentially use provided names
    if not tcr_name or tcr_name == 'autoDCRscripts-TCR':
        for field in ['v_call', 'j_call', 'junction_aa']:
            if field not in tcr_dict:
                tcr_dict[field] = ''
            if not tcr_dict[field]:
                tcr_dict[field] = 'no_' + field

        return '_'.join(['autoDCRscripts-TCR', eg_gene(tcr_dict['v_call']).replace('*', '-'),
                         eg_gene(tcr_dict['j_call']).replace('*', '-'),
                         tcr_dict['junction_aa']])

    else:
        return tcr_name


def save_json(out_path, to_save, print2stdout=True):
    """
    :param out_path: str, path to save JSON to
    :param to_save: dict to save as JSON
    :param print2stdout: bool, detailing whether to print a confirmation after saving
    :return: nothing
    """
    with open(out_path, 'w') as out_file:
        json.dump(to_save, out_file)

    if print2stdout:
        print('\tSaved to', out_path)


def read_json(in_path):
    """
    :param in_path: str, path to JSON file to read in
    :return: parsed JSON document
    """
    try:
        with open(in_path, 'r') as in_file:
            return json.load(in_file)

    except Exception:
        raise IOError(f"Unable to read in JSON file '{in_path}'.")


def translate(nt_seq):
    """
    :param nt_seq: Nucleotide sequence to be translated
    :return: corresponding amino acid sequence
    """

    aa_seq = ''
    nt_seq = nt_seq.upper()
    for i in range(0, len(nt_seq), 3):
        codon = nt_seq[i:i+3]
        if len(codon) == 3:
            try:
                aa_seq += codons[codon]
            except Exception:
                return ''

    return aa_seq


def check_features(feat_str, feat_type):
    """
    :param feat_str: str, being a description of particular feature to be processed
    :param feat_type: str telling what kind of feature it is, one of 'loci' or 'regions'
    :return: list, with each item in feat_str standardised checked/sorted against pre-coded references
    """

    if feat_type == 'loci':
        feat_list = ['A', 'B', 'G', 'D']
    elif feat_type == 'regions':
        feat_list = ['L', 'V', 'J', 'C']
    else:
        raise IOError("Inappropriate feature type given:", feat_str)

    chars = [x.upper() for x in feat_str]
    if all(x in feat_list for x in chars):
        chars.sort()
        return chars
    else:
        non_listed = '/'.join([x for x in chars if x not in feat_list])
        raise IOError(f"Unexpected {feat_type} characters detected ({non_listed}). "
                      f"Only combinations of {'/'.join(feat_list)} allowed.")


def fastafy(gene, seq_line):
    """
    :param gene: Gene symbol, extracted from the read id
    :param seq_line: Total protein primary sequence, extracted from input FASTA/generated by in silico splicing
    :return: An output-compatible FASTA entry ready for writing to file
    """
    return ">" + gene + "\n" + textwrap.fill(seq_line, 60) + "\n"


def orientation_options():
    return ['forward', 'reverse', 'both']


def protein_mode_check(desired_mode):
    """
    :param desired_mode: str of the requested mode to run
    :return: nothing, but raise an error if this mode isn't allowed for protein decombining
    """
    protein_modes = ['STDOUT', 'RETURN', 'JSON']
    if desired_mode not in protein_modes:
        raise IOError(f"TCR protein sequences can only be analysed in the following modes: {', '.join(protein_modes)}.")

def tidy_gene_list(str_gene_list, delimiter=',', resolution='allele'):
    """
    :param str_gene_list: str, detailing a list of genes or alleles split by a certain identifier
    :param delimiter: str, detailing a specific delimiter used in str_gene_list
    :param resolution: str, being 'allele' (the default) or 'gene'
    :return: str of input list, split out on the delimiter, and sorted (to standardise matching lists)
    """
    gene_list = str_gene_list.split(delimiter)
    if resolution == 'gene':
        gene_list = [x.split('*')[0] for x in gene_list]
    gene_list = list(set(gene_list))
    gene_list.sort()
    return delimiter.join(gene_list)


codons = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
          'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
          'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
          'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
          'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
          'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
          'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
          'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
          'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
          'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
          'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
          'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
          'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
          'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
          'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
          'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}

regions = {'L': 'LEADER', 'LEADER': 'L',
           'V': 'VARIABLE', 'VARIABLE': 'V',
           'J': 'JOINING', 'JOINING': 'J',
           'C': 'CONSTANT', 'CONSTANT': 'C',
            # TODO add '1' (L-PART1) , 'E' (V-EXON) ?
           }


out_headers = ['sequence_id', 'v_call', 'd_call', 'j_call', 'junction_aa', 'duplicate_count', 'sequence',
               'junction', 'rev_comp', 'productive', 'sequence_aa',
               'inferred_full_nt', 'inferred_full_aa', 'non_productive_junction_aa',
               'vj_in_frame', 'stop_codon', 'conserved_c', 'conserved_f', 'cdr3_in_limit',
               'inter_tag_seq', 'inter_tag_qual', 'umi_seq', 'umi_qual',
               'sequence_alignment', 'germline_alignment', 'v_cigar', 'd_cigar', 'j_cigar']

full_feat_headers = ['orf_in_frame', 'orf_len', 'lc_in_frame', 'found_stop', 'found_start', 'c_call', 'l_call']
    # TODO add other fields? start/codon/ORF full/in_frame?

required_ref_files = ['.tags', '.fasta']
