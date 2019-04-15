from .lowlevel import __call
from functools import partial
printf = partial(print, flush=True)


def sam_to_bam(file, keep=True):
    """
    Wrapper function of "samtools view" to convert sam into bam.

    Args:
        file: str, path-like

        keep: bool
            Keep the input file or not
    """
    # -S: input is a sam
    # -b: output is a bam
    # -h: include header section
    __call(f"samtools view -S -b -h {file} > {file[:-4]}.bam")
    if not keep:
        __call(f"rm {file}")


def bam_to_sam(file, keep=True):
    """
    Wrapper function of "samtools view" to convert bam into sam.

    Args:
        file: str, path-like

        keep: bool
            Keep the input file or not
    """
    # -h: include header section
    __call(f"samtools view -h {file} > {file[:-4]}.sam")
    if not keep:
        __call(f"rm {file}")


def fq_to_fa(file, keep=True):
    """
    Wrapper function of "seqtk" to convert fastq into fasta.

    Args:
        file: str, path-like

        keep: bool
            Keep the input file or not
    """
    if file.endswith('.fastq'):
        output = file[:-6] + '.fa'
    elif file.endswith('.fq'):
        output = file[:-3] + '.fa'

    __call(f"seqtk seq -A {file} > {output}")
    if not keep:
        __call(f"rm {file}")


def vcf_to_bcf(file, keep=True):
    """
    Wrapper function of "bcftools view" to convert VCF into BCF.

    Args:
        file: str, path-like

        keep: bool
            Keep the input file or not
    """
    # -Ou: output uncompressed bcf
    # -o <file_out>
    __call(f"bcftools view -Ou -o {file[:-4]}.bcf {file}")
    if not keep:
        __call(f"rm {file}")


def bcf_to_vcf(file, keep=True):
    """
    Wrapper function of "bcftools view" to convert BCF into VCF.

    Args:
        file: str, path-like

        keep: bool
            Keep the input file or not
    """
    # -Ov: output uncompressed vcf
    # -o <file_out>
    __call(f"bcftools view -Ov -o {file[:-4]}.vcf {file}")
    if not keep:
        __call(f"rm {file}")

