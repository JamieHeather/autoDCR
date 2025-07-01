FASTA/Q V/J/CDR3 annotation with autoDCR annotate
=================================================

At its simplest, ``autoDCR`` can be used like any other TCR annotation tool via the ``annotate`` command. When pointed to a particular FASTA or FASTQ file containing rearranged TCR DNA sequences it will loop through and detect V/J gene usage, describe CDR3 sequences that join them, and look for different sequence features from which we infer the productivity of the rearrangement (i.e. is it expected to be successfully transcribed and translated).


.. code:: bash

    autoDCR annotate -in some-file.fastq

``autoDCR`` then produces a tab separated file detailing the rearrangements it detected in the `AIRR Community standard scheme <https://docs.airr-community.org>`_ as proposed in `Vander Heiden et al., 2018 <https://doi.org/10.3389/fimmu.2018.02206>`_, supplemented with some additional custom fields.

Basic parameters
----------------

There are a few optional command line arguments that users can supply to override the default parameters. The basic parameters include:

* ``-o / --out_path``: specify a directory or file name to output the results to. Note that the default is the name of the supplied fastq file (up to the first '.' character) + '.tsv', written to whatever directory ``autoDCR`` is being run in.
* ``-or / --orientation``: specify the orientation of the DNA to search for TCRs in, choosing from forward/reverse/both (or f/r/b), default = both. Can speed up the assignation process if a particular orientation is expected due to the library preparation.
* ``-sp / --species``: name of the species (or at least the name of the species directory, see the :ref:`Generating reference data` section), used to identify the appropriate .fasta/.tags/.translate files.
* ``-dl / --deletion_limit``: upper limit of TCR V/J gene deletions to allow in a functional recombination. Default = 30.
* ``-dt / --dont_tranlate``: toggle off automatic translation and CDR3 detection.
* ``-dz / --dont_gzip``: toggle off automatic gzip compression of produced results.

Advanced options
----------------

``autoDCR`` is currently capable of several additional annotation tasks. For reference, see their specific sections:

* :ref:`Analysing full TCR transcripts`
* :ref:`Annotating TCR protein sequences`
* :ref:`Annotate single TCRs in the command line with autoDCR cli`

Things to note
--------------

* ``autoDCR`` aims to be as precise in its reporting of TCRs as possible, therefore it will list the smallest list of alleles that could have contributed to a given recombination as possible, separated by columns.

* TCRs are deemed to be 'non-productive' if they contain stop codons, don't have in-frame variable and constant domains, and/or lack certain conserved motifs. While these rearrangements might be deemed less interesting, they are useful for some applications, and so ``autoDCR`` will still try to determine their CDR3 junction. Ordinarily ``autoDCR`` determines correct frames and conserved motifs using both the V and J genes, but in cases where there are detected stop codons it will only do in reference to the conserved J region FGXG motif. If it finds one, it will then just look for the nearest upstream V region Cys residue. In these circumstances it will report two comma-separated full inferred translations in the ``inferred_full_aa`` column: one assuming the V gene is in-frame, the other assuming the J gene is.

* ``autoDCR`` reports FASTQ quality scores along with the original reads. When provided with a FASTA file instead it simply generates a length-matched dummy quality score, consisting of spaces.
