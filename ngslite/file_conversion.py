from typing import Optional
from .lowlevel import call


def sam_to_bam(
        file: str,
        keep: bool = True,
        output: Optional[str] = None) -> str:
    """
    Use cmd "samtools view" to convert sam into bam

    Args:
        file: path-like

        keep:
            Keep the input file or not

        output: path-like
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
    Use cmd "samtools view" to convert bam into sam

    Args:
        file: path-like

        keep:
            Keep the input file or not

        output: path-like
    """
    if output is None:
        output = f'{file[:-4]}.sam'

    # -h: include header section
    cmd = f'samtools view -h {file} > {output}'
    call(cmd=cmd)
    if not keep:
        call(f'rm {file}')

    return output


def fq_to_fa(
        file: str,
        keep: bool = True,
        output: Optional[str] = None) -> str:
    """
    Use cmd "seqtk" to convert fastq into fasta

    Args:
        file: path-like

        keep:
            Keep the input file or not

        output: path-like
    """
    if output is None:
        if file.endswith('.fastq'):
            output = file[:-6] + '.fa'
        elif file.endswith('.fq'):
            output = file[:-3] + '.fa'
        else:
            output = file + '.fa'

    cmd = f'seqtk seq -A {file} > {output}'
    call(cmd=cmd)
    
    if not keep:
        call(f'rm {file}')
    
    return output    


def vcf_to_bcf(
        file: str,
        keep: bool = True,
        output: Optional[str] = None) -> str:
    """
    Use cmd "bcftools view" to convert VCF into BCF

    Args:
        file: path-like

        keep:
            Keep the input file or not
            
        output: path-like
    """
    if output is None:
        output = f'{file[:-4]}.bcf'
    
    # -O u: output uncompressed bcf
    # -o <file_out>
    cmd = f'bcftools view --no-version -O u -o {output} {file}'
    call(cmd=cmd)
    
    if not keep:
        call(f'rm {file}')
        
    return output


def bcf_to_vcf(
        file: str,
        keep: bool = True,
        output: Optional[str] = None) -> str:
    """
    Use cmd "bcftools view" to convert BCF into VCF

    Args:
        file: path-like

        keep:
            Keep the input file or not

        output: path-like
    """
    if output is None:
        output = f'{file[:-4]}.vcf'

    # -O v: output uncompressed vcf
    # -o <file_out>
    cmd = f'bcftools view --no-version -O v -o {output} {file}'
    call(cmd=cmd)

    if not keep:
        call(f'rm {file}')

    return output
