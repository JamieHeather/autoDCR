Importing autoDCR to other scripts
==================================

If an application were to require V/J/CDR3 (and potentially even constant and/or leader sequence) annotation, there are various ways to import ``autoDCR`` functionality into other Python scripts. Perhaps the simplest way is to import the function which underlies ``autoDCR cli``, as this provides a versatile interface to perform different kinds of TCR annotation. Assuming reference data has been suitably processed, users can do the following:

.. code::

    # Import the necessary bits and bobs
    from autoDCRscripts import commandline as cli
    from autoDCRscripts.main import data_dir

    # Let's pick a couple of TCR sequences
    in_tcr_nt = ('GGGTTTTTCTGCTGTGGGTACGTGAGCAGGAAACATGGAGAAGAATCCTTTGGCAGCCCCATTACTAATCCTCTGGTTTCATCTTGACTGCGTGAGCAGCATACT'
                 'GAACGTGGAACAAAGTCCTCAGTCACTGCATGTTCAGGAGGGAGACAGCACCAATTTCACCTGCAGCTTCCCTTCCAGCAATTTTTATGCCTTACACTGGTCCAG'
                 'ATGGGAAACTGCAAAAAGCCCCGAGGCCTTGTTTGTAATGACTTTAAATGGGGATGAAAAGAAGAAAGGACGAATAAGTGCCACTCTTAATACCAAGGAGGGTTA'
                 'CAGCTATTTGTACATCAAAGGATCCCAGCCTGAAGACTCAGCCACATACCTCTGTGCCTTTACTTATGGAGGAAGCCAAGGAAATCTCATCTTTGGAAAAGGCAC'
                 'TAAACTCTCTGTTAAACCAAATATCCAGAACCCTGACCCTGACGTG')
    in_tcr_aa = ('DSGVTQTPKHLITATGQRVTLRCSPRSGDLSVYWYQQSLDQGLQFLIQYYNGEERAKGNIL'
                 'ERFSAQQFPDLHSELNLSSLELGDSALYFCASSARSGELFFGEGSRLTVLEDLKNVFPPL')

    # Then run the annotation command on them
    out_tcr_nt = cli.cli_annotate(tcr=in_tcr_nt, dcr_mode='full', output_mode='return', tcr_name='some-cool-nt-tcr',
                                  species='HUMAN', loci='AB', orientation='forward', barcoding=False, protein=False,
                                  genbank_mode='read', deletion_limit=30, cdr3_limit=30, data_dir=data_dir)

    out_tcr_aa = cli.cli_annotate(in_tcr_aa, 'vjcdr3', 'return', 'tcr', 'HUMAN', 'AB', 'forward',
                                  False, True, 'read', 30, 30, data_dir)


    # The different components of the detected rearrangements
    out_tcr_nt['v_call']
    'TRAV24*01'
    out_tcr_nt['j_call']
    'TRAJ42*01'
    out_tcr_nt['l_call']
    'TRAV24*01'
    out_tcr_nt['junction_aa']
    'CAFTYGGSQGNLIF'

    out_tcr_aa['v_call']
    'TRBV9*01,TRBV9*03'
    out_tcr_aa['j_call']
    'TRBJ2-2*01'
    out_tcr_aa['junction_aa']
    'CASSARSGELFF'


Things to note
--------------

* The above approach makes use of an additional ``dcr_mode``, 'return'. This is somewhat analogous to the 'json' mode, but instead of saving the output dictionary to a file, it returns it.


