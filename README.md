# autoDCR

### v0.2.7
#### Jamie Heather, MGH, 2024

### Introduction

`autoDCR` (short for **auto**matic **D**e**c**ombinato**r**) is a python script to perform T cell receptor (TCR) gene annotation of DNA sequencing reads. This is inspired by and in part built upon the core functionality of [Decombinator](https://github.com/innate2adaptive/Decombinator), the TCR analysis software developed by the Chain lab at UCL (in part by me). It uses a similar conceptual framework of using fast Aho-Corasick tries to search for the presence of 'tag' sequences in DNA reads, and use these to identify V and J TCR genes. 

The first big difference between the scripts is what tags are looked for. `Decombinator` aims to find single CDR3-proximal tags (i.e. in the 3' of V genes or 5' of J genes) which aim to identify individual TCR genes - ideally covering all of the alleles of that gene - producing an extremely concise trie which can very rapidly search large datasets. Instead, `autoDCR` generates overlapping tags across the entirety of all alleles for all genes for both TRA and TRB, producing a much larger trie. Tags will then potentially cover multiple genes, but genes are called based on the intersection between genes/alleles covered in all tags which they matched. 

This sacrifices the speed of `Decombinator`, but allows for greater resolution of TCR genes and alleles, while using and retaining sequence information from a greater portion of the input TCR read. This tag approach also makes the automation of tag generation much simpler, making it much easier to either update `autoDCR` when new TCR alleles are discovered, or even to apply it to situations like non-natural repertoires, or for novel species for which TCR information has recently become available.

  * [Original `Decombinator` paper DOI: 10.1093/bioinformatics/btt004](https://doi.org/10.1093/bioinformatics/btt004)
  * [V4 update paper DOI: 10.1093/bioinformatics/btaa758](https://doi.org/10.1093/bioinformatics/btaa758)

**Note** that `autoDCR` is in development, and does not offer as full a suite of functions as does `Decombinator`. Most especially it currently lacks the UMI based error-correction functionality of `Decombinator`.

### Generating `autoDCR` tag files

In order to use `autoDCR` you must first make the .tag and .translate files for the species you care about:

* Download all V and J genes for your given species from IMGT-GENE/DB.
    * Note that currently `autoDCR` is designed with TRA and TRB genes, and it's recommended that if you plan to regularly look for both you just generate files with both TRA and TRB combined. 
      * However if you know you're only interested in one specific locus you can download just the genes for that.
      * Note that in humans all TRDV genes can be found recombined with TRAJ (even those genes not necessarily labelled as TRAVx/DVx), so I recommend including all TRDV genes for alpha/beta recombination searches.
    * Save as a FASTA file, retaining the pipe- ('`|`') delimited IMGT header information.
    * Try to include date and IMGT release number (and any other pertitent information) in the file name, as the output files will be renamed to the species in question (E.g. 'human.fasta' or 'mouse.fasta'). Input file name information will be retained in the log file.
  * The recommended way to achieve all this is by using the [IMGTgeneDL](https://github.com/JamieHeather/IMGTgeneDL) package:

```bash
pip install IMGTgeneDL
# Example download of all human V and J genes for TRA and TRB
IMGTgeneDL -s Homo+sapiens -L AB -v -j
```

* Generate appropriate tag and translation files using `generate-tag-files.py`, e.g. using the above `IMGTgeneDL` example:
```
python generate-tag-files.py -in 2024-03-13_202410-7_Homo+sapiens_AB_VJ.fasta 
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
* If this is the first time generating these files for a particular species, be sure to check the .log file.
    * Conserved C/FGXG motifs are detected using regexes manually produced by generating positional weight matrices of all four human TCR loci, which allows it to detect even non-canonical residues at the conserved positions, in allele sequences that may be incomplete. 
    * However these motifs may be less conserved between species, and so if the log file shows that there are many predicted-functional genes not being found with high-confidence motifs users may wish to inspect and correct the `regexes` dictionary in `generate-tag-files.py`.

An example of the products of this process using the human TRA/TRB loci (downloaded from IMGT/GENE-DB on 2024-03-13, release `202410-7`) has been included in this repo.

### Running `autoDCR`

At its simplest, `autoDCR` can be run simply by pointing it to a particular FASTQ file that contains rearranged TCRs. (Note that it also works on FASTA files: in that case it will simply generate a fake quality score consisting of spaces.)

```bash
python autoDCR.py -fq FILE.fastq
```

`autoDCR` then produces a tab separated file detailing the rearrangements it detected in the [AIRR Community standard schema](https://docs.airr-community.org) as proposed in [Vander Heiden *et al.*, 2018](https://doi.org/10.3389/fimmu.2018.02206), supplemented with some additional custom fields. 

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
* `-jv / --jump_values`: toggle on outputting the V/J 'jump' values (positions of the outermost V/J tag matches). Defaults to on when using the `-ad` option.
* `-ad / --allele_discovery`: toggle on a mode that permits rudimentary downstream novel V/J allele inference (see below).

#### Note on non-productive rearrangements

TCRs are deemed to be 'non-productive' if they contain stop codons, don't have in-frame variable and constant domains, and/or lack certain conserved motifs. While these rearrangements might be deemed less interesting, they are useful for some applications, and so `autoDCR` will still try to determine their CDR3 junction. Ordinarily `autoDCR` determines correct frames and conserved motifs using both the V and J genes, but in cases where there are detected stop codons it will only do in reference to the conserved J region FGXG motif. If it finds one, it will then just look for the nearest upstream V region C residue. In these circumstances it will report two comma-separated full inferred translations: one assuming the V gene is in-frame, the other assuming the J gene is.

#### Novel allele inference with `inferTCR`

In its determination of what alleles are used in a given rearrangement, `autoDCR` categorises exactly which tag sequences derived from which genes are present along a read. This makes it possible to perform rudimentary allele detection by looking for TCRs within an individual that have a consistent break in a run of tag matches, as an undescribed single nucleotide variant will result in two contiguous/overlapping tags not matching. Users can run `autoDCR` with the `-ad` flag set, which will add additional fields to the output data categorising such mismatches.

Note that feature is highly experimental, requiring high-quality TCR sequencing reads that span the entirety of the variable domain, which have been produced from a single donor's T cells. It also should only be run using tags that have been generated with default parameters in `generate-tag-files` (20 nt length, overlapping 10 nt). Note that it will also only infer novel alleles that are relatively commonly used in the donor being analysed, and those with distinguishing variant sites that are a) not present in the first/last 10 nucleotides of a gene and b) contained within a 10 bp window. 

```bash
# Run autoDCR with allele discovery mode on
python3 autoDCR.py -fq TCRs.fasta -ad

# Then run inferTCR on the output file produced
python3 inferTCR.py -in TCRs_infer-alleles.tsv.gz 
```

This then produces a few files:

* *_infer-alleles_alleles.fasta  
  * FASTA sequences of the potential novel alleles detects in the input file
* *_infer-alleles_sourcedata.tsv
  * The actual raw rows from the *_infer-alleles.tsv.gz file produced by `-ad` mode `autoDCR`
* *_infer-alleles_summary.tsv
  * A summary file tabulating the properties of the alleles which related to their being inferred as potentially novel 
