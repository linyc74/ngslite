from .fasta import FastaParser, FastaWriter


import os
import subprocess
from functools import partial
printf = partial(print, flush=True)


def __call(cmd):
    printf('CMD: ' + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        printf(inst)


def __rename_contig_id(genbank, contig_dict):
    """
    Args:
        genbank: str, path-like

        contig_dict: dict
            For example: {'contig_1_52671bp': 'assembly=control_contigs;id=NODE_400_length_52671_cov_411.055515;len=52671'}
                           (short header)        (original header)
    """
    with open(genbank) as reader:
        with open('temp.gb', 'w') as writer:
            for line in reader:
                if line.startswith('LOCUS       '):
                    key = line.split()[1]
                    writer.write(
                        line.replace(key, contig_dict[key])
                    )
                else:
                    writer.write(line)
    os.remove(genbank)
    os.rename('temp.gb', genbank)


def prokka(fasta, outdir, kingdom, locus_tag='LOCUS', proteins='',
           metagenome=True, evalue=1e-6, threads=4, log=None, keep=True):
    """
    A wrapper function for Prokka 1.12 (rapid bacterial genome annotation)

    This wrapper function takes a fasta and outputs a genbank file with the name <outdir>.gb

    Args
        fasta: str, path-like
            The input fasta to be annotated

        outdir: str, path-like
            The output directory and the genbank file name

        kingdom: str
            'Archaea', 'Bacteria', 'Mitochondria' or 'Viruses'

        proteins: str, path-like
            FASTA or GBK file to use as 1st priority

        metagenome: bool
            Improve gene predictions for highly fragmented genomes

        evalue: float or int
            Similarity e-value cut-off

        threads: int
            Number of CPUs to be used (0 = all)

        log: str, path-like
            The log file for stderr

        keep: bool
            Whether or not to keep all data in the <outdir>
    """
    # Create a temp fasta file to shorten the header
    #   because prokka does not take long headers
    # Also create a dictionary to index the shortened header
    contig_dict = dict()
    with FastaParser(fasta) as parser:
        with FastaWriter('temp.fa') as writer:
            i = 0
            for head, seq in parser:
                new_head = f"contig_{i}_{len(seq)}bp"
                writer.write(new_head, seq)
                contig_dict[new_head] = head
                i += 1

    # Prokka command line
    meta = ['', '--metagenome '][metagenome]
    evalue = format(evalue, 'f')
    if proteins:
        proteins = f"--proteins {proteins} "
    if log is not None:
        log = f" 2> {log}"
    else:
        log = f" 2> {outdir}.log"
    cmd = f"prokka --outdir {outdir} --prefix {outdir} --locustag {locus_tag} --kingdom {kingdom} --evalue {evalue} --cpus {threads} {proteins}{meta}temp.fa{log}"
    __call(cmd)
    os.remove('temp.fa')

    # In the outdir folder, write a tsv file to index shortened contig headers to the original headers
    with open(f"{outdir}/contig_index.tsv", 'w') as fh:
        for key, val in contig_dict.items():
            fh.write(f"{key}\t{val}\n")

    # Copy the genbank file out from the output directory
    __call(f"cp {outdir}/{outdir}.gbk {outdir}.gb")

    # Rename the contig header back to the original headers using <contig_dict>
    __rename_contig_id(f"{outdir}.gb", contig_dict)

    # Remove temporary files in the output dir
    if not keep: __call(f"rm -r {outdir}")
