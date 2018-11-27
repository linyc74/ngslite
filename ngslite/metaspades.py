import subprocess


from .fasta import *


def __call(cmd):
    print('CMD: ' + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        print(inst)


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

    print(f"Retrieve contigs (>= {min_contig_length} bp) from {output}/contigs.fasta -> {output}.fa")
    parser = FastaParser(f"{output}/contigs.fasta")
    writer = FastaWriter(f"{output}.fa")
    while True:
        head, seq = parser.next()
        if head is None:
            break
        if len(seq) >= min_contig_length:
            writer.write(head, seq)
    parser.close()
    writer.close()

