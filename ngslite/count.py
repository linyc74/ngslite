from .filetools import gzip
import subprocess
import os
from functools import partial
printf = partial(print, flush=True)


def count_reads(file, mapped=True):
    """
    Count the reads in a file

    Supported format: fasta, fastq and sam/bam, compressed gz is also supported

    Args:
        file: str, path-like

        mapped: bool
            If True, only count reads that are mapped in a sam or bam file
            If False, count all reads

    Returns: int
        The number of reads or sequences contained in the input file
    """
    is_gz = False
    if file.endswith('.gz'):
        gzip(file, keep=True)
        file = file[:-3]
        is_gz = True

    ftypes = ['fq', 'fastq', 'fa', 'fasta', 'sam', 'bam']
    assert file.split('.')[-1] in ftypes, 'Invalid input file extention.'

    with open(file, 'r') as fh:
        # FASTQ: Use number of lines to count the reads in fastq file
        if file.endswith('.fq') or file.endswith('.fastq'):
            i = 0
            while fh.readline() != '':
                i += 1
            if i % 4 != 0:
                printf('Warning: Number of lines in the fastq file is not multiple of 4.')
            count = int(i / 4)

        # FASTA: Use '>' to count the reads in fasta file
        elif file.endswith('.fa') or file.endswith('.fasta'):
            i = 0
            while True:
                line = fh.readline()
                if line == '': break
                if line.startswith('>'):
                    i += 1
            count = i

        # SAM/BAM: Use 'samtools view -c' command to count reads
        else:
            mapped = ['', '-F 4 '][mapped]
            cmd = f"samtools view -c {mapped}{file}"
            count = int(subprocess.check_output(cmd, shell=True))

    if is_gz:
        os.remove(file)

    return count

