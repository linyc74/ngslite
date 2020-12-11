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
    Args:
        ref:
            The reference fasta

        fq1:
            The read-1 fastq

        fq2:
            The read-2 fastq. If none, use <fq1> for unpaired mapping

        sam:
            The output SAM file

        threads:
            Number of CPUs

        score:
            Don’t output alignment with score lower than <score>
    """
    call(f'bwa index {ref}')

    if fq2 is None:
        cmd = f'bwa mem -t {threads} -T {score} {ref} {fq1} > {sam}'  # unpaired
    else:
        cmd = f'bwa mem -t {threads} -T {score} {ref} {fq1} {fq2} > {sam}'  # paired
    call(cmd)

    call('rm ref.*')
