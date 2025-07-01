.. autoDCR documentation master file, created by
   sphinx-quickstart on Mon Apr 21 10:42:35 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

autoDCR
=======

A modified version of Decombinator for flexible TCR annotation
--------------------------------------------------------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   referencedata
   annotate
   tcrtranscripts
   cli
   tcrproteins
   importing


``autoDCR`` (short for **auto**\ matic **D**\ e\ **c**\ ombinato\ **r**\ ) is a python script to perform T cell receptor (TCR) gene annotation of DNA sequencing reads. This is inspired by and in part built upon the core functionality of `Decombinator <https://github.com/innate2adaptive/Decombinator>`_, the TCR analysis software developed by the Chain lab at UCL (in part by me). It uses a similar conceptual framework of using fast Aho-Corasick tries to search for the presence of 'tag' sequences in DNA reads, and use these to identify V and J TCR genes.

The first big difference between the scripts is what tags are looked for. Standard ``Decombinator`` (V1 to V4) aims to find single CDR3-proximal tags (i.e. in the 3' of V genes or 5' of J genes) which identify individual TCR genes - ideally covering all of the alleles of that gene - producing an extremely concise trie which can very rapidly search large datasets. Instead, ``autoDCR`` generates overlapping tags across the *entirety* of **all alleles** for all genes for both TRA and TRB, producing a much larger trie. Tags will then potentially cover multiple genes, but genes are called based on the intersection between genes/alleles covered in all tags which they matched.

This sacrifices the speed of ``Decombinator``, but allows for greater resolution of TCR genes and alleles, while using and retaining sequence information from a greater portion of the input TCR read. This tag approach also makes the automation of tag generation much simpler, making it much easier to either update ``autoDCR`` when new TCR alleles are discovered, or even to apply it to situations like non-natural repertoires, or for novel species for which TCR information has recently become available.

  * `Original Decombinator paper DOI: 10.1093/bioinformatics/btt004 <https://doi.org/10.1093/bioinformatics/btt004>`_
  * `V4 update paper DOI: 10.1093/bioinformatics/btaa758 <https://doi.org/10.1093/bioinformatics/btaa758>`_

**Note** that ``autoDCR`` is in active development and should be considered highly experimental, and does not offer as full a suite of functions as does ``Decombinator`` (nor seek to). Most especially it currently lacks the UMI based error-correction functionality of ``Decombinator``. Instead it is intended to supply a flexible framework of adaptable TCR annotation-related functions that can be tweaked for dataset-specific applications, including a few which I don't feel are necessarily well catered to with existing tools.
