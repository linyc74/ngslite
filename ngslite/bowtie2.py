from .lowlevel import __call
from functools import partial
printf = partial(print, flush=True)


def bowtie2_mapping(ref, fq1, sam, fq2=None):
    """
    Wrapper function of bowtie2 mapping.

    Args:
        ref: str, path-like object
            The reference fasta

        fq1: str, path-like object
            The read-1 fastq

        fq2: str, path-like object
            The read-2 fastq. If none, use <fq1> for unpaired mapping.

        sam: str, path-like object
            The output SAM file
    """
    # Build the .bt2 index files
    __call(f"bowtie2-build {ref} ref > bowtie2_build_{ref}.log")

    log = sam[:-4]+'.log'
    if fq2 is None:
        # Unpaired mapping
        cmd = f"bowtie2 -x ref -U {fq1} -S {sam} 2> {log}"
    else:
        # Paired-end mapping
        cmd = f"bowtie2 -x ref -1 {fq1} -2 {fq2} -S {sam} 2> {log}"
    __call(cmd)

    # Remove the .bt2 index files
    __call('rm *.bt2')

