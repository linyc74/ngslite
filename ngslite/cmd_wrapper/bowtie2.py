from typing import Optional
from ..lowlevel import call


def bowtie2_mapping(
        ref: str,
        fq1: str,
        sam: str,
        fq2: Optional[str] = None):
    """
    Args:
        ref:
            The reference fasta

        fq1:
            The read-1 fastq

        fq2:
            The read-2 fastq
            If none, use <fq1> for unpaired mapping

        sam:
            The output SAM
    """
    call(f'bowtie2-build {ref} ref > {ref}_bowtie2_build.log')

    log = sam[:-4]+'.log'
    if fq2 is None:
        cmd = f'bowtie2 -x ref -U {fq1} -S {sam} 2> {log}'  # unpaired
    else:
        cmd = f'bowtie2 -x ref -1 {fq1} -2 {fq2} -S {sam} 2> {log}'  # paired
    call(cmd)

    call('rm *.bt2')
