from .lowlevel import __call
from functools import partial
printf = partial(print, flush=True)


def bwa_mapping(ref, fq1, sam, fq2=None, threads=4, score=30):
    """
    Wrapper function of BWA mapping.

    Args:
        ref: str, path-like object
            The reference fasta

        fq1: str, path-like object
            The read-1 fastq

        fq2: str, path-like object
            The read-2 fastq. If none, use <fq1> for unpaired mapping.

        sam: str, path-like object
            The output SAM file

        threads: int,
            # of CPU threads

        score: int,
            Don’t output alignment with score lower than <score>
    """
    # Build index
    __call(f"bwa index {ref}")

    if fq2 is None:
        # Unpaired
        cmd = f"bwa mem -t {threads} -T {score} {ref} {fq1} > {sam}"
    else:
        # Paired-end
        cmd = f"bwa mem -t {threads} -T {score} {ref} {fq1} {fq2} > {sam}"
    __call(cmd)

    # Remove the index files
    __call('rm ref.*')

