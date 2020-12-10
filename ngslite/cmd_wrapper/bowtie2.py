from typing import Optional
from ..lowlevel import call


def bowtie2_mapping(
        ref: str,
        fq1: str,
        sam: str,
        fq2: Optional[str] = None):
    """
    Wrapper function of bowtie2 mapping

    Args:
        ref: path-like
            The reference fasta

        fq1: path-like
            The read-1 fastq

        fq2: path-like
            The read-2 fastq. If none, use <fq1> for unpaired mapping

        sam: path-like
            The output SAM file
    """
    # Build the .bt2 index files
    call(f"bowtie2-build {ref} ref > {ref}_bowtie2_build.log")

    log = sam[:-4]+'.log'
    if fq2 is None:
        # Unpaired mapping
        cmd = f"bowtie2 -x ref -U {fq1} -S {sam} 2> {log}"
    else:
        # Paired-end mapping
        cmd = f"bowtie2 -x ref -1 {fq1} -2 {fq2} -S {sam} 2> {log}"
    call(cmd)

    # Remove the .bt2 index files
    call('rm *.bt2')