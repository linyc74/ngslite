from .fasta import *
from .lowlevel import __call
from functools import partial
printf = partial(print, flush=True)


def metaspades(fq1, fq2, output, min_contig_length=1000, threads=16, memory=250):
    """
    Wrapper function of the command "spades.py --meta".

    Args:
        fq1: str, path-like, the read-1 fastq file

        fq2: str, path-like, the read-2 fastq file

        output: str, path-like
            The output directory of metaspades, also the output fasta file containing assembled contigs

        min_contig_length: int
            Minimum (inclusive) contig length in bp

        threads: int, # of CPU cores

        memory: int, # of Gb of RAM
    """
    __call(f"spades.py --meta -1 {fq1} -2 {fq2} -o {output} --threads {threads} --memory {memory} > {output}.log")

    printf(f"Retrieve contigs (>= {min_contig_length} bp) from {output}/contigs.fasta -> {output}.fa",
          flush=True)

    with FastaParser(f"{output}/contigs.fasta") as parser:
        with FastaWriter(f"{output}.fa") as writer:
            for head, seq in parser:
                if len(seq) >= min_contig_length:
                    writer.write(head, seq)
