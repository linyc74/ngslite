import os
from typing import Optional, Dict
from ..fasta import FastaParser, FastaWriter
from ..lowlevel import call
from ..filetools import get_temp_path


def prokka(
        fasta: str,
        outdir: str,
        kingdom: str,
        locus_tag: str = 'LOCUS',
        proteins: str = '',
        metagenome: bool = True,
        evalue: float = 1e-6,
        threads: int = 4,
        log: Optional[str] = None,
        keep: bool = True):
    """
    A wrapper function for Prokka 1.12 (rapid bacterial genome annotation)

    This wrapper function takes a fasta and outputs a genbank file with the name <outdir>.gb

    Args
        fasta:
            The nucleotide fasta to be annotated

        outdir:
            The output directory and the genbank file name

        kingdom:
            'Archaea', 'Bacteria', 'Mitochondria' or 'Viruses'

        locus_tag

        proteins: path-like
            FASTA or GBK file to use as 1st priority

        metagenome:
            Improve gene predictions for highly fragmented genomes

        evalue:
            Similarity e-value cut-off

        threads:
            Number of CPUs to be used (0 = all)

        log: path-like
            The log file for stderr

        keep:
            Whether or not to keep all data in the <outdir>
    """
    # Create a temp fasta file to shorten the header because prokka does not take long headers
    # Also create a dictionary to index the shortened header
    temp_fa = get_temp_path('temp', '.fa')  # for example, temp000.fa
    contig_dict = dict()
    with FastaParser(fasta) as parser:
        with FastaWriter(temp_fa) as writer:
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

    cmd = f"prokka --outdir {outdir} --prefix {outdir} --locustag {locus_tag} \
--kingdom {kingdom} --evalue {evalue} --cpus {threads} {proteins}{meta}{temp_fa}{log}"

    call(cmd)
    os.remove(temp_fa)

    # In the outdir folder, write a tsv file to index shortened contig headers to the original headers
    with open(f"{outdir}/contig_index.tsv", 'w') as fh:
        for key, val in contig_dict.items():
            fh.write(f"{key}\t{val}\n")

    # Copy the genbank file out from the output directory
    call(f"cp {outdir}/{outdir}.gbk {outdir}.gb")

    # Rename the contig header back to the original headers using <contig_dict>
    _rename_contig_id(f"{outdir}.gb", contig_dict)

    # Remove files in the output dir if keep == False
    if not keep:
        call(f"rm -r {outdir}")


def _rename_contig_id(
        genbank: str,
        contig_dict: Dict[str, str]):
    """
    Args:
        genbank: path-like

        contig_dict:
            {'contig_1_52671bp': 'assembly=control_contigs;id=NODE_400_length_52671_cov_411.055515;len=52671'}
               (short header)        (original header)
    """
    temp_gb = get_temp_path('temp', '.gb')  # for example, temp000.gb
    with open(genbank) as reader:
        with open(temp_gb, 'w') as writer:
            for line in reader:
                if line.startswith('LOCUS       '):
                    key = line.split()[1]
                    writer.write(
                        line.replace(key, contig_dict[key])
                    )
                else:
                    writer.write(line)
    os.remove(genbank)
    os.rename(temp_gb, genbank)
