import os
import subprocess
from .lowlevel import printf
from .filetools import get_temp_path,gzip


def _unzip(file: str) -> str:

    file_type = file[:-3].split('.')[-1]

    temp = get_temp_path(
        prefix=f'{file[:-3]}_count',
        suffix=f'.{file_type}')

    file = gzip(file=file, keep=True, output=temp)

    return file


def count_reads(file: str, mapped: bool = True) -> int:
    """
    Count the reads in a file

    Supported format: fasta, fastq and sam/bam, compressed gz is also supported

    Args:
        file: path-like

        mapped: bool
            If True, only count reads that are mapped in a sam or bam file
            If False, count all reads

    Returns:
        The number of reads or sequences contained in the input file
    """

    is_gz = file.endswith('.gz')

    if is_gz:
        file = _unzip(file=file)

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
