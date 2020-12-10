from typing import Optional
from ..lowlevel import call


def bwa_mapping(
        ref: str,
        fq1: str,
        sam: str,
        fq2: Optional[str] = None,
        threads: int = 4,
        score: int = 30):
    """
    Wrapper function of BWA mapping

    Args:
        ref: path-like
            The reference fasta

        fq1: path-like
            The read-1 fastq

        fq2: path-like
            The read-2 fastq. If none, use <fq1> for unpaired mapping

        sam: path-like
            The output SAM file

        threads:
            Number of CPUs

        score:
            Don’t output alignment with score lower than <score>
    """
    # Build index
    call(f"bwa index {ref}")

    if fq2 is None:
        # Unpaired
        cmd = f"bwa mem -t {threads} -T {score} {ref} {fq1} > {sam}"
    else:
        # Paired-end
        cmd = f"bwa mem -t {threads} -T {score} {ref} {fq1} {fq2} > {sam}"
    call(cmd)

    # Remove the index files
    call('rm ref.*')