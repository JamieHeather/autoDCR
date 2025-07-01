Generating reference data
=========================

In order to work, ``autoDCR`` requires a few different reference files (which differ slightly depending on the mode being being run). These files are automatically generated using the ``autoDCR refs`` command, which in turn uses the `stitchr <https://jamieheather.github.io/stitchr/>`_\ -style output of `IMGTgeneDL <https://github.com/JamieHeather/IMGTgeneDL>`_, a tool used to automatically download and parse germline TCR allele data from `IMGT/GENE-DB <https://imgt.org/genedb/doc>`_.

``autoDCR`` data directory
--------------------------

``autoDCR`` generates its own data directory during installation, and will download and manage its content there itself. If everything is working as it should, end users shouldn't need to manually interact with their data directory much, but if it is required the location of this database can be printed by running ``autoDCR dd``.


``autoDCR refs``
----------------

``refs`` is the subcommand that grabs and manages the reference data.

Basic referencing
~~~~~~~~~~~~~~~~~

``autoDCR`` can theoretically be used to analyse any of the four common TCR loci (alpha, beta, gamma, delta) for which IMGT has sufficient information (see `this page of the stitchr documentation <https://jamieheather.github.io/stitchr/inputdata.html>`_ for more details). The majority of TCRseq research (including my own) deals with human alpha/beta data, which is reflected in most of the ``autoDCR`` default parameters.

``autoDCR`` identifies TCR genes used in a rearrangement by scanning reads with an Aho-Corasick trie built from short 'tag' sequences, themselves produced by taking overlapping subsequences from all allele gene sequences in the provided reference. In order to use ``autoDCR`` in its simplest form - i.e. annotating V, J, and CDR3 regions in TCRseq reads - one need only run:

.. code:: bash

    # for human data
    autoDCR refs

    # for other species
    autoDCR refs -s [species common name]
    autoDCR refs -s mouse

When you run that, ``autoDCR`` grabs and processes the necessary data for its style of annotation.

* By default this is done for all V and J genes of TRA and TRB genes together, so that the same process can be applied to all A/B human reads regardless of locus.
* Running ``autoDCR refs`` under its default conditions will therefore generate a 'HUMAN' folder in the data directory.
* That will contain several kinds of files, including:
    * General files relating to the establishment of the species-specific data directory:
        * ``imgt-data.fasta``, containing all of the unfiltered reads from the IMGT reference for the requested species (humans, by default).
        * Locus-specific filtered versions of that (``TRA.fasta, TRB.fasta`` etc).
        * Some automatically-produced motif files produced by ``IMGTgeneDL`` (``C-`` and ``J-regionmotifs.tsv``)
        * ``data-production-date.tsv``, which contains information about the date and release of IMGT/GENE-DB used (which should be retained for reporting when publising any results generated with a given reference).
    * Several files prefixed ``AB_JV`` (representing 'alpha/beta V and J genes'):
        * A ``*.fasta`` file, containing the subset of genes used for a specific trie.
        * A ``*.tags`` file, containing the specific tags used to populate the trie, and the specific alleles they are found in. These contain several fields:
            * tag sequences (by default 20-mers), used in the actual tag searching
            * tag jump values (comma-separated where applicable), which indicates the position in the gene of that tag
            * comma-separated lists of the alleles covered by each tag
        * A ``*.translate`` file, detailing automatically-inferred conserved V and J gene motifs used to help translate TCRs and identify CDR3 junctions, specifically with the fields:
            * TCR allele
            * relative position
            * amino acid found at the conserved residue position
            * Note that in both file types, positions are counted forwards from the 5' of the V or backwards from the 3' end of the J.
            * Not only is it conceptually easier to count backwards for the Js (due to the start of the J being obscured by V(D)J recombination), but the conserved FGXG motifs tend to fall at specific positions relative to the end of the Js, at least in functional genes, presumably due to evolutionary constraints.

        * A ``*.log`` file, which is not used in subsequent TCR annotation efforts but which contains metadata produced during the automatic processing of a given species' genes (see Notes section below).

* You can specify exactly which loci ``autoDCR`` files are generated for using the ``--loci / -L`` flag, e.g.:

.. code:: bash

    # Gamma/delta TCR references
    autoDCR refs -L GD

    # Single chain A/B references (to process individually)
    autoDCR refs -L A
    autoDCR refs -L B -sd

* Note the use of the ``skip_download / --sd`` flag
    * If not used, each use of ``autoDCR refs`` will over-write the previous data.
    * As some more advanced features of ``autoDCR`` require multiple iterations of ``refs`` to be applied to the same dataset, users must take care to ensure the correct order of commands using ``-sd`` is observed.


Novel allele supplementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As part of my interest in the impact of human TCR gene polymorphism, I maintain a `resource tracking putative novel human TCR alleles (not featured in IMGT/GENE-DB) identified in the literature <https://github.com/JamieHeather/novel-tcr-alleles/>`_. If users wish to include these alleles in their database, they can run the following command:

.. code:: bash

    autoDCR refs -nv

Note that this uses a function adapted from ``stitchr`` (version 0.3.0), which `is described in slightly more detail in the stitchr docs <https://jamieheather.github.io/stitchr/inputdata.html#using-novel-inferred-human-tcr-alleles>`_.


Additional region referencing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to typical V/J/CDR3 annotation, ``autoDCR`` also has some modes which permit additional annotation of the leader sequence and constant regions of a rearranged TCR.

In order to extract the necessary information for these modes, users must first generate the corresponding reference information, achieved by running ``autoDCR`` refs with the appropriate ``--regions / -r`` flag with the ``skip_download / -sd`` flag, after running regular reference production (+/- novel alleles) for that species. E.g.:

.. code:: bash

    # For human AB TCRs, first you must have run
    autoDCR refs
    # Or, if including novel alleles
    autoDCR refs -nv

    # Then to add the constant and leader files
    autoDCR refs -sd -r CL

TCR protein sequence referencing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``autoDCR`` is also capable of applying its TCR annotation functions to translated polypeptide sequences, as might be found in say repurposed structural data. In order to generate the corresponding files, ``autoDCR refs`` needs to be run with the boolean ``--protein / -aa`` flag. Users also need to supply a decreased sliding window length (which dictates the length of the tags) via the ``--sliding_window / -sw`` flag as below (using a length of 10, which seems to work reasonably well):

.. code:: bash

    # For human AB TCRs, first you must have acquired standard nt V/J sequences

    # Then generate amino acid versions from that
    autodcr refs -aa -sw 10

Use of the ``--protein / -aa`` flag automatically applies the ``skip_download`` flag, as it requires the corresponding nucleotide data to have been downloaded already.

Also note that this is one of the (`even more`) experimental features. As such it currently can only be used for standard V/J/CDR3 annotation, and may behave oddly with some less-common parameters.

Quick full human referencing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following commands will establish all of the necessary files for the current breadth of ``autoDCR`` functionality, at least as relates to human alpha/beta TCR analysis:

.. code:: bash

    autoDCR refs -nv
    autoDCR refs -sd -r CL
    autodcr refs -aa -sw 10

Things to note
--------------

* If you use ``pip`` to uninstall ``autoDCR`` will likely result in the deletion of this folder, so if you are likely to need the contents of the directory (e.g. if you have used a particular configuration to analyse data for publication) it is suggested that you make a prior of this directory before uninstalling.

* As all TRDV genes can be found recombined with TRAJ (even those genes not necessarily labelled as TRAVx/DVx, at least in humans), ``autoDCR`` automatically considers all TRDV genes when looking for alpha chain recombinations.

* When generating tags for a particular species for the first time, be sure to check the ``.log`` file produced in the relevant data directory,
    * Conserved C/FGXG motifs are detected using regexes manually produced by generating positional weight matrices of all four human TCR loci, which allows it to detect even non-canonical residues at the conserved positions, in allele sequences that may be incomplete.
    * However these motifs may be less conserved between species, and so if the log file shows that there are many predicted-functional genes not being found with high-confidence motifs users may wish to inspect and correct the `regexes` dictionary, stored in the ``regexes.json`` file located in the data directory.

* If for some reason multiple data directories for a given species are required to be maintained simultaneously, you can navigate to the data directory and rename the species folder (e.g. you could have 'HUMAN' and 'HUMAN2'), which could then be referred to specifically using the ``--species / -s`` of the different ``autoDCR`` commands.
    * Each directory needs to be generated with a recognised common name however, as this is what ``IMGTgeneDL`` uses as a reference to pull out the right reads from IMGT/GENE-DB.
