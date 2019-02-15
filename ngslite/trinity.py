from functools import partial
printf = partial(print, flush=True)


def __call(cmd):
    printf('CMD: ' + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        printf(inst)


def trinity(fq1, fq2, output, threads=16, memory=128, min_contig_length=1000, normalize_reads=False, min_kmer_cov=1, kmer_size=25):
    """
    Args:
        fq1: str, path-like
            The read-1 fastq file

        fq2: str, path-like
            The read-2 fastq file

        output: str, path-like
            The output directory of Trinity

        threads: int
            # of CPU cores

        memory: int
            # of Gb of RAM

        min_contig_length: int
            Minimum (inclusive) contig length in bp

        normalize_reads: bool
            Normalize the read coverage by setting the max cov to 200 (default of Trinity),
                i.e. k-mers with counts > 200 will be arbitrarily set to 200

        min_kmer_cov: int
            Miminum k-mer coverage to be used for transcriptome assembly
            The default setting (1) of Trinity could overload the memory or its underlining Java implementation
            Thus I usually set it to 3

        kmer_size: int
            k-mer size (bp) used to construct the De Bruijn graph
    """
    no_normalize_reads = ['--no_normalize_reads', ''][normalize_reads]
    __call(f"Trinity --seqType fq --left {fq1} --right{fq2} --output {output} --CPU {threads} --max_memory {memory}G --min_contig_length {min_contig_length} {no_normalize_reads}--min_kmer_cov {min_kmer_cov} --KMER_SIZE {kmer_size} --no_version_check")

    # Copy the fasta file containing assembled contigs
    __call(f"cp {output}/Trinity.fasta {output}.fa")

