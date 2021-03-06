import os
from typing import Optional
from ..lowlevel import call
from ..fasta import FastaParser


def simple_phylogeny(
        fasta: str,
        output: str,
        clustering: str = 'NJ',
        log: Optional[str] = None):
    """
    Make a tree file (newick/phylip format) from the aligned <fasta>
        using the command "clustalw -tree"

    Args:
        fasta: path-like
            The input fasta file, sequences must be aligned

        output: path-like
            The output tree file

        clustering:
            clustering method, 'NJ' or 'UPGMA'

        log: path-like
            The log file, default None -> <output>.log
    """
    with FastaParser(fasta) as parser:
        head, seq = parser.next()
        if set(seq.upper()) == set('ACGT-'):
            type_ = 'DNA'
        else:
            type_ = 'PROTEIN'

    if log is None:
        log = output + '.log'

    call(f'clustalw -tree -infile={fasta} -type={type_} -CLUSTERING={clustering} > {log}')

    # Rename the .ph file automatically generated by clustalw
    ph_file = fasta[:fasta.rfind('.')] + '.ph'
    os.rename(ph_file, output)