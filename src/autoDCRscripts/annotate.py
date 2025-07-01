import collections as coll
import os
import sys
from time import time
from . import autoDCRfunctions as fxn

def vjcdr3_annotate(mode, in_path, out_path, species, loci, orientation, protein, barcoding,
                    deletion_limit, cdr3_limit, dont_translate, dont_gzip, data_dir):
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

    mode = mode.upper()

    loci = fxn.check_features(loci, 'loci')
    if not protein:
        mol_type = 'nt'
    else:
        mol_type = 'aa'

    reference_data = fxn.import_tcr_info(species, loci, 'JV', mol_type, data_dir)
    reference_data = fxn.import_translate_info(reference_data)
    headers = fxn.out_headers
    if mode == 'FULL':
        extra_refdat = fxn.import_tcr_info(species, loci, 'CL', 'nt', data_dir)
        for field in fxn.full_feat_headers:
            headers.insert(4, field)

    dcr_parameters = {'mode': mode,
                      'orientation': orientation,
                      'deletion_limit': deletion_limit,
                      'cdr3_limit': cdr3_limit,
                      'mol_type': mol_type}
    # TODO add in don't translate?

    # Determine where to save the results
    if out_path == '[input-file-name].tsv':
        out_path = in_path[:in_path.rfind('.')].split('/')[-1] + '.tsv'
    elif not out_path.endswith('.tsv'):
        out_path += '.tsv'
    if not dont_gzip:
        out_path += '.gz'

    counts = coll.Counter()
    start = time()

    # Then loop through the input file, analysing TCRs as we go
    with fxn.opener(in_path, 'r') as in_file, fxn.opener(out_path, 'w') as out_file:

        # Initialise the output file with the header
        out_file.write('\t'.join(headers) + '\n')
        out_str = []
        for read_id, seq, qual in fxn.readfq(in_file):

            # Pad empty quality scores for FASTA files
            if not qual:
                qual = ' ' * len(seq)

            counts['reads'] += 1

            # Figure out the relevant parts of the read for decombining, then search
            read, read_qual, bc, bc_qual = fxn.sort_read_bits(seq, qual, '')  # TODO add barcoding

            # TODO break it down into different functions
                # TODO 1) find tags 2) call rearrangements 3) translate

            tcr_check = fxn.dcr(read, read_qual, reference_data, dcr_parameters, headers)

            if tcr_check:

                # TODO barcoding
                # if input_args['barcoding']:
                #     tcr_check['umi_seq'] = bc
                #     tcr_check['umi_qual'] = bc_qual

                counts['rearrangements'] += 1
                tcr_check['sequence_id'] = read_id




                # # TODO full - l+c search
                if mode == 'FULL':
                    tcr_check = fxn.find_full_feats(tcr_check, extra_refdat, dcr_parameters)

                # Remove in-process gene region labeling
                tcr_check = fxn.tidy_output(tcr_check, headers)

                line_out = '\t'.join([str(tcr_check[x]) for x in headers])

                # TODO move
                # if discover:
                #     if tcr_check['v_mismatches'] or tcr_check['j_mismatches']:
                #         counts['mismatched_germlines'] += 1

                out_str.append(line_out)

                # Bulk write the results out once there's a sufficient chunk (to prevent this getting too big in memory)
                if len(out_str) % 10000 == 0:
                    out_file.write('\n'.join(out_str) + '\n')
                    out_str = []

        # Then write out any leftover calls
        out_file.write('\n'.join(out_str))
        # TODO fix duplicate counts (if desired?)

    end = time()
    time_taken = end - start
    print("Took", str(round(time_taken, 2)), "seconds")
    print("Found", str(counts['rearrangements']), "rearranged TCRs in", str(counts['reads']), "reads "
          "(~" + str(round(counts['rearrangements']/counts['reads'] * 100)) +"%)")

    # if discover:
    #     print("Of these,", str(counts['mismatched_germlines']), "showed discontinuous tag matches "
    #                                                             "and were kept aside for inference of potential new alleles.")

    # TODO sort summary output (maybe into YAML or JSON?)
