from typing import List, Optional
from ..lowlevel import call


def sort_bam(file: str, keep: bool = False) -> str:
    """
    Args:
        file:
            BAM file

        keep:
            Keep the input file or not
    """
    file_out = f'{file[:-4]}_sorted.{file[-3:]}'
    call(f'samtools sort {file} > {file_out}')
    if keep:
        return file_out
    else:
        call(f'rm {file}')
        call(f'mv {file_out} {file}')
        return file


def index_bam(file: str) -> str:
    call(f'samtools index {file}')
    return f'{file}.bai'


def sam_to_indexed_bam(file: str, keep: bool = True) -> str:
    sam_to_bam(file, keep=keep)
    bam = file[:-4] + '.bam'
    sort_bam(bam, keep=False)
    index_bam(bam)
    return bam


def subset_bam_regions(
        file: str,
        regions: List[str],
        output: Optional[str] = None,
        replace_input: bool = False):
    """
    Args:
        file:
            The input BAM file

        regions:
            Each str is a region of the reference genome, e.g.
                chr1            chromosome 1
                chr3:1000-2000  chromosome 3 from 1000th (inclusive) to 2000th (inclusive) base

        output:
            The output BAM file
            If None, add subscript '_subset' to the input BAM file

        replace_input:
            If True, delete the input <file> and rename the output as the input <file>
            Overrides the <output> file name
    """
    # Convert regions list into a string
    # ['chr1', 'chr2:1001-2000'] -> ' chr1 chr2:1001-2000'
    regions = ' ' + ' '.join(regions)  # space-separated regions

    if output is None:
        output = file[:-len('.bam')] + '_subset.bam'

    # -b: output is a bam
    # -h: include header section
    call(f'samtools view -b -h {file}{regions} > {output}')

    if replace_input:
        call(f'rm {file}')
        call(f'mv {output} {file}')


def remove_unmapped(file: str, keep: bool = True):
    """
    Remove unmapped reads from an input SAM/BAM file
        by using the option '-F 4' of 'samtools view' command

    Args:
        file:
            The input SAM or BAM file

        keep:
            Keep the input file or not
    """
    b = ['', '-b '][file.endswith('.bam')]  # Output file is bam or not
    file_out = f'{file[:-4]}_remove_unmapped.{file[-3:]}'

    # -b: output is a bam
    # -h: include header section
    # -F 4: NOT including the flag 'read unmapped'
    call(f'samtools view -h -F 4 {b}{file} > {file_out}')

    if not keep:
        call(f'rm {file}')
        call(f'mv {file_out} {file}')


def sam_to_bam(
        file: str,
        keep: bool = True,
        output: Optional[str] = None) -> str:
    """
    Args:
        file:
            SAM file path

        keep:
            Keep the input file or not

        output:
            BAM file path
    """
    if output is None:
        output = f'{file[:-4]}.bam'

    # -S: input is a sam
    # -b: output is a bam
    # -h: include header section
    cmd = f'samtools view -S -b -h {file} > {output}'
    call(cmd=cmd)

    if not keep:
        call(f'rm {file}')

    return output


def bam_to_sam(
        file: str,
        keep: bool = True,
        output: Optional[str] = None) -> str:
    """
    Args:
        file:
            BAM file path

        keep:
            Keep the input file or not

        output:
            SAM file path
    """
    if output is None:
        output = f'{file[:-4]}.sam'

    # -h: include header section
    cmd = f'samtools view -h {file} > {output}'
    call(cmd=cmd)
    if not keep:
        call(f'rm {file}')

    return output
