# autoDCR

### v0.1.0
#### Jamie Heather, MGH, 2021

### Introduction

`autoDCR` (short for **auto**matic **D**e**c**ombinato**r**) is a python script to perform T cell receptor (TCR) gene annotation of DNA sequencing reads. This is inspired by and in part built upon the core functionality of [Decombinator](https://github.com/innate2adaptive/Decombinator), the TCR analysis software developed by the Chain lab at UCL (in part by me). It uses a similar conceptual framework of using fast Aho-Corasick tries to search for the presence of 'tag' sequences in DNA reads, and use these to identify V and J TCR genes. 

The first big difference between the scripts is what tags are looked for. `Decombinator` aims to find single CDR3-proximal tags (i.e. in the 3' of V genes or 5' of J genes) which aim to identify individual TCR genes - ideally covering all of the alleles of that gene - producing an extremely concise trie which can very rapidly search large datasets. Instead, `autoDCR` generates overlapping tags across the entirety of all alleles for all genes for both TRA and TRB, producing a much larger trie. Tags will then potentially cover multiple genes, but genes are calls based on the intersection between genes/alleles covered in all tags which they matched. 

This sacrifices the speed of `Decombinator`, but allows for greater resolution of TCR genes and alleles, while using and retaining sequence information from a greater portion of the input TCR read. This tag approach also makes the automation of tag generation much simpler, making it much easier to either update `autoDCR` when new TCR alleles are discovered, or even to apply it to novel species for which TCR information has become available.

  * [Original `Decombinator` paper DOI: 10.1093/bioinformatics/btt004](https://doi.org/10.1093/bioinformatics/btt004)
  * [Recent V4 update paper DOI: 10.1093/bioinformatics/btaa758](https://doi.org/10.1093/bioinformatics/btaa758)

**Note** that `autoDCR` is in development, and is not a full suite of functions as is `Decombinator`. Most especially it currently lacks the UMI based error-correction functionality of `Decombinator`.

### Generating `autoDCR` tag files

In order to use `autoDCR` you must first make the .tag and .translate files for the species you care about:

* Download all V and J genes for your given species from IMGT-GENE/DB.
    * Note that currently `autoDCR` is designed with TRA and TRB genes, and it's recommended that if you plan to regularly look for both you just generate files with both TRA and TRB combined. 
      * However if you know you're only interested in one specific locus you can download just the genes for that.
      * Note that in humans all TRDV genes can be found recombined with TRAJ (even those genes not necesessarily labelled as TRAVx/DVx), so I recommend including all TRDV genes for alpha/beta recombination searches.
    * Save as a FASTA file, retaining the '|' delimited header information.
    * Try to use a one word species name (without any special characters, including '.' characters) to name this file.
    * E.g. 'human.fasta' or 'mouse.fasta'.
* Generate appropriate tag and translation files using `generate-tag-files.py`, e.g.:
```
python generate-tag-files.py -in human.fasta
```
* This will generate two tab-separated files:
  * \[x\].tags, describing the tags used for searching TCRs, with the fields:
    * tag sequence (by default a 20-mer), used in the actual tag searching.
    * tag jump value (comma-separated where applicable), which indicates the position in the gene of that tag. Unlike regular `Decombinator`, `autoDCR` no longer uses these values in the TCR annotation process, but the values are retained to aid in manually verifying tag generation.
    * Comma-separated list of the alleles covered by this tag
  * \[x\].translate, used to identify the inferred location of the CDR3 junction-ending conserved residues, with the fields:
    * TCR allele
    * relative position
    * amino acid found at the conserved residue position
* Note that in both file types, positions are counted forwards from the 5' of the V or backwards from the 3' end of the J.
  * Not only is it conceptually easier to count backwards for the Js (due to the start of the J being obscured by V(D)J recombination), but the conserved FGXG motifs tend to fall at specific positions relative to the end of the Js, at least in functional genes, presumably due to evolutionary constraints. 
* If this this is the first time generating these files for a particular species, be sure to check the .log file.
    * Be particularly on the lookout for genes that are flagged with 'Manually verify': these are genes for which the automatic conserved motif detection was unsuccessful.

An example of this using the human TRA/TRB loci (downloaded from IMGT/GENE-DB on 2021-08-14) has been included in this repo.

### Running `autoDCR`

At its simplest, `autoDCR` can be run simply by pointing it to a particular FASTQ file that contains rearranged TCRs. (Note that it also works on FASTA files: in that case it will simply generate a fake quality score consisting of spaces.)

```bash
python autoDCR.py -fq FILE.fastq
```

`autoDCR` then produces a tab separated file detailing the rearrangements it detected in the [AIRR Community standard schema](https://docs.airr-community.org) as proposed in [Vander Heiden *et al.*, 2018](https://doi.org/10.3389/fimmu.2018.02206), s
                   'v_jump', 'j_jump' # This is the non-standard modificationupplemented with some additional custom fields. 

There are a few optional command line arguments that users can supply to override the default parameters:

* `-o / --out_path`: specify a directory or file name to output the results to. Note that the default is the name of the supplied fastq file (up to the first '.' character) + '.tsv', written to whatever directory `autoDCR` is being run in.
* `-or / --orientation`: specify the orientation of the DNA to search for TCRs in, choosing from forward/reverse/both (default = both). Can speed up the assignation process if a particular orientation is expected due to the library preparation.
* `-sp / --species`: name of the species (or at least the name used in the `generate-tag-files` process), used to identify the appropriate .fasta/.tags/.translate files. Default = huma.fasta/.tags/.translate n.
* `-dd / --data_dir`: specify the directory containing the .fasta/.tags/.translate files, if not located in the current working directory.
* `-bc / --barcoding`: toggles on the 'barcoding' option, which will split the read (and quality scores) into two strings, presuming that the 5' portion corresponds to non-TCR UMI/barcoding sequences (which could feed in to a downstream error-correction process).
* `-bl / --bclength`: the barcode length - i.e. the number of bases from the start of each read which are separated out to their own field for downstream analysis - to be used in conjunction with the -bc flag. Default = 42.
* `-dl / --deletion_limit`: upper limit of deletions to allow in a functional recombination. Default = 30.
* `-dt / --dont_tranlate`: toggle off automatic translation and CDR3 detection.
* `-dz / --dont_gzip`: toggle off automatic gzip compression of produced results.
