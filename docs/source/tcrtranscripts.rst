Analysing full TCR transcripts
==============================

In addition to performing regular annotation of the rearranged V and J genes, and the intervening hypervariable CDR3, ``autoDCR`` is able under certain circumstances to look for some additional features that should also be present in a fully mature TCR mRNA: the leader sequence, and the constant gene. If analysing suitably amplified/sequenced (or assembled) transcriptomic data, the leader sequence (the TCR chain signal peptide, containing the leader intron) will be spliced together and found immediately upstream of the V-REGION, while the constant region should be spliced on to the 3' end of the J-REGION.

This is achieved by changing the ``--dcr_mode / -m`` flag from its default value (``vjcdr3``, describing its search for the V, J, and CDR3 sequences) to ``full``. Note that in order to be capable of doing so, additional steps need to be performed when generating the necessary files as described in the :ref:`Generating reference data` section. Doing this generates a second set of tags, specific for the constant and leader regions (C/L), which will be applied in a second round of Aho-Corasick string matching to a read if a rearranged V-J gene is detected.

This can be performed for both ``autoDCR annotate`` and ``autoDCR cli`` (being default for the latter command) as described below:

.. code::

    autoDCR annotate -in some-file -m full

    tcr="TACCGGCTATCTGTACGGGGACAGATACAGAAGACCCCTCCGTCATGCAGCATCTGCCATGAGCATCGGCCTCCTGTGCTGTGCAGCCTTGTCTCTCCTGTGGGCAGGTCCAGTGAATGCTGGTGTCACTCAGACCCCAAAATTCCAGGTCCTGAAGACAGGACAGAGCATGACACTGCAGTGTGCCCAGGATATGAACCATGAATACATGTCCTGGTATCGACAAGACCCAGGCATGGGGCTGAGGCTGATTCATTACTCAGTTGGTGCTGGTATCACTGACCAAGGAGAAGTCCCCAATGGCTACAATGTCTCCAGATCAACCACAGAGGATTTCCCGCTCAGGCTGCTGTCGGCTGCTCCCTCCCAGACATCTGTGTACTTCTGTGCCAGCAGACTGGACAGGGAGTACGAGCAGTACTTCGGGCCGGGCACCAGGCTCACGGTCACAGAGGACCTGAAAAACGTGTTCCCCTCT"

    # Use the default 'full' autoDCR mode to analyse a long TCR read
    autoDCR cli $tcr

    TCR regions detected!
            orientation     forward
            productive      yes
            v_call  TRBV6-5*01
            j_call  TRBJ2-7*01
            junction_aa     CASRLDREYEQYF
            l_call  TRBV6-5*01



Things to note
--------------

* When using 'full' mode with ``annotate``, additional columns will be generated in the output TSV file, which aim to give a wider-inference as to the predicted productivity of a given rearrangement than just looking at the V-J region section.

* Currently 'full' TCR analysis is only valid when analysing nucleotide sequences (and thus cannot be used with the protein analysis capabilities as described in the :ref:`Annotating TCR protein sequences` section).

