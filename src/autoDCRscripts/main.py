# -*- coding: utf-8 -*-

import os
import sys
import typer
from typing_extensions import Annotated
from . import referencedata as refdat
from . import autoDCRfunctions as fxn
from . import annotate as tcr
from . import commandline as cmd
from . import genotype as geno
from . import discover as disc
from . import __version__


# Ensure correct importlib-resources function imported
if sys.version_info < (3, 9):
    import importlib_resources                              # PyPI
else:
    import importlib.resources as importlib_resources       # importlib.resources


app = typer.Typer(context_settings={"help_option_names": ["-h", "--help"]})


data_files = importlib_resources.files("autoDCRdata")
data_dir = os.path.dirname(str(data_files / '__init__.py'))


@app.callback()
def callback():
    """
    autoDCR: a modified version of Decombinator to explore different experimental functions\n
    https://jamieheather.github.io/autoDCR/
    """


@app.command()
def refs(species: Annotated[str, typer.Option("--species", "-s",
                  help="The name of the species to generate data for (upper case common name)")] = 'HUMAN',
         sliding_window: Annotated[int, typer.Option("--sliding_window", "-sw",
                  help="Length of the sliding window for tag generation")] = 20,
         loci: Annotated[str, typer.Option("--loci", "-L",
                  help="Loci to generate autoDCR files for")] = "AB",
         regions: Annotated[str, typer.Option("--regions", "-r",
                  # TODO does this need to be specifiable? can have a sep option for 'mode' which could populate
                  help="Regions to generate autoDCR files for")] = "VJ",
         novel: Annotated[bool, typer.Option("--novel", "-nv",
                  help="Optional flag to also download putative novel alleles (human only)")] = False,
         protein: Annotated[bool, typer.Option("--protein", "-aa",
                  help="Generate amino acid tags for analysing TCR protein sequences")] = False,
         skip_download: Annotated[bool, typer.Option("--skip_download", "-sd",
                  help="Optional flag to skip downloading new reference data, just remaking files")] = False,
         source: Annotated[str, typer.Option("--source", "-src",
                  help="Germline reference source to use for annotation. IMGT, OGRDB, or path to a FASTA")] = "IMGT",
         name: Annotated[str, typer.Option("--name", "-n",
                  help="Optional string to rename the reference data folder (mandatory with the '-src' flag.")] = ""):
    """
    Automatically downloads and/or generates the necessary reference FASTA,
    tag, and translation files for autoDCR TCR sequence analysis
    """

    typer.echo(f"BUILDING REFERENCE DATA for species '{species}', with sliding window length {sliding_window}")
    refdat.get_reference_data(species, sliding_window, loci, regions, skip_download,
                              novel, protein, source, name, data_dir)


@app.command()
def annotate(in_file: Annotated[str, typer.Option('-in', '--in_file',
                help="FASTA or FASTQ DNA file containing TCR reads.")],
             mode: Annotated[str, typer.Option('-m', '--mode',
                help="AutoDCR mode: vjcdr3 (regular), full (+leaders/constants")] = 'vjcdr3',  # TODO add more modes?
             out_path: Annotated[str, typer.Option('-o', '--out_path',
                help="'Optionally specify path to output file.'")] = '[input-file-name].tsv',
             species: Annotated[str, typer.Option('-s', '--species',
                help="The name of the species to search for TCRs for (upper case common name)")] = 'HUMAN',
             loci: Annotated[str, typer.Option('-L', '--loci',
                help="TCR loci to search for rearrangements of.")] = 'AB',
             orientation: Annotated[str, typer.Option('-or', '--orientation', autocompletion=fxn.orientation_options,
                help="Specify the orientation to search in (forward/reverse/both).")] = 'both',
             protein: Annotated[bool, typer.Option("--protein", "-aa",
                help="Analyse TCR protein sequences.")] = False,
             barcoding: Annotated[bool, typer.Option('-bc', '--barcoding',
                # TODO needs adding, placeholder currently
                # TODO will need barcode length/specification parameters
                help="Flag for barcoded libraries.")] = False,
             deletion_limit: Annotated[int, typer.Option('-dl', '--deletion_limit',
                help="Upper limit of allowable deletions for each of V/J in a recombination.")] = 30,
             cdr3_limit: Annotated[int, typer.Option('-cl', '--cdr3_limit',
                help="Upper limit of allowable length of translated CDR3 junctions. Set to 0 for no limit.")] = 30,
             tag_coverage: Annotated[bool, typer.Option('-tc', '--tag_coverage',
                help="Include additional fields relating to position of outermost tag matches")] = False,
             dont_translate: Annotated[bool, typer.Option('-dt', '--dont_translate',
                help="Stop the automatic translation of TCRs.")] = False,
             dont_gzip: Annotated[bool, typer.Option('-dz', '--dont_gzip',
                help="Stop the output FASTQ files automatically being compressed with gzip")] = False):
    """
    Perform standard V/J/CDR3 TCR annotation on a specified FASTA or FASTQ file, for a set locus or loci.
    """
    typer.echo(f"Looking for TCRs in {in_file}")
    tcr.vjcdr3_annotate(mode, in_file, out_path, species, loci, orientation, protein, barcoding,
                             deletion_limit, cdr3_limit, tag_coverage, dont_translate, dont_gzip, data_dir)


@app.command()
def discover(in_file: Annotated[str, typer.Option('-in', '--in_file',
                help="FASTA or FASTQ DNA file containing TCR reads.")],
             out_path: Annotated[str, typer.Option('-o', '--out_path',
                help="'Optionally specify path to output file.'")] = '[input-file-name].allele',
             species: Annotated[str, typer.Option('-s', '--species',
                help="The name of the species to search for TCRs for (upper case common name)")] = 'HUMAN',
             loci: Annotated[str, typer.Option('-L', '--loci',
                help="TCR loci to search for rearrangements of.")] = 'AB',
             orientation: Annotated[str, typer.Option('-or', '--orientation', autocompletion=fxn.orientation_options,
                help="Specify the orientation to search in (forward/reverse/both).")] = 'both',
             barcoding: Annotated[bool, typer.Option('-bc', '--barcoding',
                # TODO needs adding, placeholder currently
                help="Flag for barcoded libraries.")] = False,
             dont_gzip: Annotated[bool, typer.Option('-dz', '--dont_gzip',
                help="Stop the output FASTQ files automatically being compressed with gzip")] = False):
    """
    Perform standard V/J/CDR3 TCR annotation on a specified FASTA or FASTQ file, for a set locus or loci.
    """
    typer.echo(f"Looking for TCRs in {in_file}")
    disc.discover(in_file, out_path, species, loci, orientation, barcoding, dont_gzip, data_dir)


@app.command()
def genotype(in_file: Annotated[str, typer.Option('-in', '--in_file',
                help="TSV file of rearrangements produced by 'autoDCR annotate.")],
             ref: Annotated[str, typer.Option('-rf', '--ref',
                help="Name of original TCR reference used during annotation (e.g. HUMAN). Required.")],
             mode: Annotated[str, typer.Option('-m', '--mode',
                help="AutoDCR mode: vjcdr3 (regular), full (+leaders/constants")] = 'vjcdr3',
             out_path: Annotated[str, typer.Option('-o', '--out_path',
                help="Optionally specify path to output file.")] = fxn.genotype_default,
             subset_loci: Annotated[str, typer.Option('-sL', '--subset_loci',
                help="TCR loci to consider if making a genotype-specific reference subset.")] = 'AB',
             allow_partial: Annotated[bool, typer.Option('-p', '--allow_partial',
               help="Allow partial transcripts to be included in the genotyping (not recommended).")] = False,
             read_threshold: Annotated[int, typer.Option('-rt', '--read_threshold',
                help="Minimum number of reads an allele must be detected in to be considered during genotyping.")] = 1
             ):
    """
    Aims to infer the genotype of a donor, based on the occurrence of gene-specific tags in an `autoDCR annotate` run.
    """
    geno.genotype(in_file, mode, out_path, ref, subset_loci, allow_partial, read_threshold, data_dir)


@app.command()
def cli(tcr: str,
        dcr_mode: Annotated[str, typer.Option('-m', '--dcr_mode',
            help="AutoDCR mode: vjcdr3 (regular), full (+leaders/constants")] = 'full',
        output_mode: Annotated[str, typer.Option('-om', '--output_mode',
            help="Output mode: stdout (terminal), gb (Genbank), json")] = 'stdout',
        tcr_name: Annotated[str, typer.Option('-n', '--tcr_name',
            help="Name of TCR for output purposes")] = 'autoDCR-TCR',
        species: Annotated[str, typer.Option('-s', '--species',
            help="The name of the species to search for TCRs for (upper case common name)")] = 'HUMAN',
        loci: Annotated[str, typer.Option('-L', '--loci',
            help="TCR loci to search for rearrangements of.")] = 'AB',
        orientation: Annotated[str, typer.Option('-or', '--orientation', autocompletion=fxn.orientation_options,
            help="Specify the orientation to search in (forward/reverse/both).")] = 'both',
        barcoding: Annotated[bool, typer.Option('-bc', '--barcoding',
            # TODO needs adding, placeholder currently
            help="Flag for barcoded libraries. (Placeholder, not currently active)")] = False,
        genbank_mode: Annotated[str, typer.Option('-gm', '--genbank_mode',
            help="GenBank output optional mode: read, inferred, or tags.")] = 'read',
        protein: Annotated[bool, typer.Option("--protein", "-aa",
            help="Analyse TCR protein sequences")] = False,
        deletion_limit: Annotated[int, typer.Option('-dl', '--deletion_limit',
            help="Upper limit of allowable deletions for each of V/J in a recombination.")] = 30,
        cdr3_limit: Annotated[int, typer.Option('-cl', '--cdr3_limit',
            help="Upper limit of allowable length of translated CDR3 junctions. Set to 0 for no limit.")] = 30):
    """
    Annotate the rearrangement of a single TCR directly in the command line.
    """

    cmd.cli_annotate(tcr, dcr_mode, output_mode, tcr_name, species, loci, orientation, barcoding,
                     protein, genbank_mode, deletion_limit, cdr3_limit, data_dir)


@app.command()
def dd(ls: Annotated[bool, typer.Option('-ls', '--list',
       help="List the available reference directories in the autodcr data directory.")] = False,
):
    """
    Print the location of the autodcr data directory
    """
    # TODO Add options to either archive or delete reference folders

    if not ls:
        print(data_dir)
    else:
        ref_dirs = [x for x in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, x))]
        ref_dirs = [x for x in ref_dirs if x[0] not in ['_', '.']]
        ref_dirs.sort()
        for ref in ref_dirs:
            print(ref)
    return data_dir


@app.command()
def version():
    """
    Print the autoDCR package version number.
    """
    print(__version__)
    return __version__
