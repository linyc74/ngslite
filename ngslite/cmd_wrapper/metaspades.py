from ..lowlevel import call, printf
from ..fasta import FastaParser, FastaWriter


def metaspades(
        fq1: str,
        fq2: str,
        output: str,
        min_contig_length: int = 1000,
        threads: int = 16,
        memory: int = 250):
    """
    Wrapper function of the command "spades.py --meta".

    Args:
        fq1:
            The read-1 fastq file

        fq2:
            The read-2 fastq file

        output:
            The output directory of metaspades
            The output fasta file containing assembled contigs

        min_contig_length:
            Minimum (inclusive) contig length in bp

        threads:
            Number of CPUs

        memory:
            Number of Gb of RAM
    """
    call(f'spades.py --meta -1 {fq1} -2 {fq2} -o {output} --threads {threads} --memory {memory} > {output}.log')

    printf(f'Retrieve contigs (>= {min_contig_length} bp) from {output}/contigs.fasta -> {output}.fa')

    with FastaParser(f'{output}/contigs.fasta') as parser:
        with FastaWriter(f"{output}.fa") as writer:
            for head, seq in parser:
                if len(seq) >= min_contig_length:
                    writer.write(head, seq)
