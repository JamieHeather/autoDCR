import collections as coll
import os
import sys
from . import autoDCRfunctions as fxn


def cli_annotate(tcr, dcr_mode, output_mode, tcr_name, species, loci, orientation, barcoding,
                 protein, genbank_mode, deletion_limit, cdr3_limit, data_dir):
    # TODO docstring

    dcr_mode, output_mode, genbank_mode = [x.upper() for x in [dcr_mode, output_mode, genbank_mode]]

    # Establish the necessary input data and parameters
    if not protein:
        mol_type = 'nt'
        if not fxn.dna_check(tcr):
            raise IOError("Non-DNA character detected in input TCR string - only A/C/G/T/N characters allowed.")
    else:
        mol_type = 'aa'
        dcr_mode = 'VJCDR3'
        fxn.protein_mode_check(output_mode)


    loci = fxn.check_features(loci, 'loci')
    reference_data = fxn.import_tcr_info(species, loci, 'JV', mol_type, data_dir)
    reference_data = fxn.import_translate_info(reference_data)
    headers = fxn.out_headers
    fields = ['v_call', 'j_call', 'junction_aa']

    if dcr_mode == 'FULL':  # TODO or if it detects a CL dataset in the directory?

        extra_refdat = fxn.import_tcr_info(species, loci, 'CL', mol_type, data_dir)
        fields = ['v_call', 'j_call', 'junction_aa', 'l_call', 'c_call']
        for field in fxn.full_feat_headers:
            headers.insert(4, field)

    # TODO diagnosis mode flag? Output all tag matches (updating, calling it TAGS now)

    dcr_parameters = {'mode': dcr_mode,
                      'orientation': orientation,
                      'deletion_limit': deletion_limit,
                      'cdr3_limit': cdr3_limit,
                      'mol_type': mol_type}

    qual = ' ' * len(tcr)

    # Figure out the relevant parts of the read for decombining, then search
    read, read_qual, bc, bc_qual = fxn.sort_read_bits(tcr, qual, '')  # TODO add barcoding

    tcr_check = fxn.dcr(read, read_qual, reference_data, dcr_parameters, headers)

    if tcr_check:
        name = fxn.get_output_name(tcr_name, tcr_check)

        if dcr_mode == 'FULL':
            tcr_check = fxn.find_full_feats(tcr_check, extra_refdat, dcr_parameters)

        tcr_check = fxn.tidy_output(tcr_check, headers)

        if output_mode not in ['RETURN']:
            print("\nTCR regions detected!")

        if output_mode.upper() == 'STDOUT':
            if not protein:
                if tcr_check['rev_comp'] == 'F':
                    print('\torientation\tforward')
                else:
                    print('\torientation\treverse')

            print(f"\tproductive\t{tcr_check['productive']}")

            for field in fields:  # TODO add in some reporting as to whether certain fields are present?
                if tcr_check[field]:
                    print('\t' + field + '\t' + tcr_check[field])
                else:
                    print('\t' + field + '\tNot detected')

            if output_mode not in ['RETURN']:
                print('')  # Visually pad the output

        elif output_mode.upper() == 'JSON':
            fxn.save_json(name + '.json', tcr_check)


        elif output_mode.upper() == 'GB':
            # Using the GenBank output functions from stitchr v0.3.0

            if genbank_mode.upper() not in ['READ', 'INFERRED', 'TAGS']:
                raise IOError("Unexpected GenBank mode option:", genbank_mode)

            description = 'TCR annotated with autoDCRscripts'  # TODO add version/ref details etc, put in function
            out_path = name   # TODO remove or add teaks?

            if genbank_mode.upper() in ['READ', 'TAGS']:
                full_sequence = tcr_check['sequence']
            elif genbank_mode.upper() in ['INFERRED']:
                full_sequence = tcr_check['inferred_full_nt']

            # TODO this is for 'read' mode, need to sort for inferred (diagnostic = just plot all tags)
            gb_feats = []

            if not protein:
                tcr_check = fxn.determine_gene_boundaries(tcr_check, reference_data, extra_refdat)
                regions = ['v', 'j', 'l', 'c']
            else:
                # TODO fix -- needs the determine_gene_boundaries to plot GB currently, but I haven't/done plan to add CL stuff for proteins!
                regions = ['v', 'j']

            for region in regions:
                if tcr_check[region + '_call']:
                    gb_feats.append((tcr_check[region + '_label'], '',
                                     [tcr_check[region + '_seq']], [tcr_check[region + '_start']]))

                if genbank_mode.upper() == 'TAGS':
                    if region + '_matching_tags' in tcr_check.keys():
                        for tagx in range(len(tcr_check[region + '_matching_tags'])):
                            tag = tcr_check[region + '_matching_tags'][tagx]
                            gb_feats.append((region.upper() + ' tag#' + str(tagx + 1), '', [tag[0]], [tag[1]]))

            # TODO add insert?
            if tcr_check['junction_aa']:
                if not tcr_check['junction_aa'].startswith('no_'):
                    gb_feats.append(('CDR3 junction: ' + tcr_check['junction_aa'], '', [tcr_check['junction']], ''))

            genbank_params = {'sequence_name': out_path, 'full_sequence': full_sequence,
                              'description': description, 'topology': 'linear', 'features': gb_feats,
                              'save_dir_path': './', 'species': species, 'numbering': False,
                              'plot_multi': True, 'division': 'SYN', 'journal': 'autoDCRscripts'}

            fxn.output_genbank(**genbank_params)

        return tcr_check

    else:
        if output_mode not in ['RETURN']:
            print("No TCR detected.")
        return