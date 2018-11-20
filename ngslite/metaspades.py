def __call(cmd):
    print('CMD: ' + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        print(inst)


def metaspades(fq1, fq2, output, threads=16, memory=250):
    """
    Args:
        fq1: str, path-like, the read-1 fastq file

        fq2: str, path-like, the read-2 fastq file

        output: str, path-like
            The output directory of metaspades, also the output fasta file containing assembled contigs

        threads: int, # of CPU cores

        memory: int, # of Gb of RAM
    """
    __call(f"spades.py --meta -1 {fq1} -2 {fq2} -o {output} --threads {threads} --memory {memory}")

    # Copy the fasta file containing assembled contigs
    __call(f"cp {output}/contigs.fasta {output}.fa")

