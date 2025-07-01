Annotating TCR protein sequences
================================

While the vast majority of TCR analyses are concerned with gDNA or cDNA sequences, there are some rare applications that require annotation of V/J/CDR3 regions of TCR `protein` sequences. However tools capable of processing such data in a comparable manner as nucleotide data are hard to find (e.g. while `IgBLAST can determine V genes it doesn't seem to be able to annotate CDR3 junctions or J genes <https://www.ncbi.nlm.nih.gov/igblast/>`_). However tasks like re-annotation of structural data that doesn't have rearrangement information contained in its metadata are becoming more common, as people re-purpose such data for informing TCR antigen predictors.

The two major ``autoDCR`` V/J/CDR3 annotation subcommands - ``autoDCR annotate`` and ``autoDCR cli`` are able to process such sequences, through generation of tag tries based off translated TCR protein germline sequences. Note that in order to do so the proper references must have been generated in the correct order, as described in the :ref:`Generating reference data` section.

In order to do so, users need supply the boolean ``--protein / -aa`` flag to their ``autoDCR annotate`` or ``autoDCR cli`` commands. E.g. the following examples, using some sequences taken from PDB TCR-pMHC structures:

.. code::

    # Inspect the file of TCR polypeptide chains
    cat prot-tcrs.fasta

    >3MV7_4|Chain D|alpha chain of the TK3 TCR|Homo sapiens (9606)
    QVTQSPEALRLQEGESSSLNCSYTVSGLRGLFWYRQDPGKGPEFLFTLYSAGEEKEKERLKATLTKKESFLHITAPKPEDSATYLCAVQDLGTSGSRLTFGEGTQLTVNPNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPS
    >3MV7_5|Chain E|beta chain of the TK3 TCR|Homo sapiens (9606)
    DSGVTQTPKHLITATGQRVTLRCSPRSGDLSVYWYQQSLDQGLQFLIQYYNGEERAKGNILERFSAQQFPDLHSELNLSSLELGDSALYFCASSARSGELFFGEGSRLTVLEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD

    >2AK4_4|Chains D, I, N, S[auth T]|SB27 T cell receptor alpha chain|Homo sapiens (9606)
    HMAQKVTQAQTEISVVEKEDVTLDCVYETRDTTYYLFWYKQPPSGELVFLIRRNSFDEQNEISGRYSWNFQKSTSSFNFTITASQVVDSAVYFCALSGFYNTDKLIFGTGTRLQVFPNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESS
    >2AK4_5|Chains E, J, O[auth P], T[auth U]|SB27 T cell receptor beta chain|Homo sapiens (9606)
    HMNAGVTQTPKFQVLKTGQSMTLQCAQDMNHNSMYWYRQDPGMGLRLIYYSASEGTTDKGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASPGLAGEYEQYFGPGTRLTVTEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD

    # Then analyse this with annotate
    autoDCR annotate -in prot-tcrs.fasta -aa

    Looking for TCRs in results/prot-tcrs.fa
    Took 0.0 seconds
    Found 4 rearranged TCRs in 4 reads (~100%)

    # Or with cli
    trb="HMNAGVTQTPKFQVLKTGQSMTLQCAQDMNHNSMYWYRQDPGMGLRLIYYSASEGTTDKGEVPNGYNVSRLNKREFSLRLESAAPSQTSVYFCASPGLAGEYEQYFGPGTRLTVTEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD"
    autoDCR cli $trb --protein

    TCR regions detected!
            productive      yes
            v_call  TRBV6-1*01
            j_call  TRBJ2-7*01
            junction_aa     CASPGLAGEYEQYF


Things to note
--------------

* Currently the protein version of ``autoDCR`` only works in standard 'vjcdr3' mode, not 'full' (i.e. it cannot be used to detect leader or constant regions).






