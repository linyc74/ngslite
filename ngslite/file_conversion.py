from .lowlevel import call


def sam_to_bam(file: str, keep: bool = True):
    """
    Wrapper function of "samtools view" to convert sam into bam

    Args:
        file: path-like

        keep:
            Keep the input file or not
    """
    # -S: input is a sam
    # -b: output is a bam
    # -h: include header section
    call(f"samtools view -S -b -h {file} > {file[:-4]}.bam")
    if not keep:
        call(f"rm {file}")


def bam_to_sam(file: str, keep: bool = True):
    """
    Wrapper function of "samtools view" to convert bam into sam

    Args:
        file: path-like

        keep:
            Keep the input file or not
    """
    # -h: include header section
    call(f"samtools view -h {file} > {file[:-4]}.sam")
    if not keep:
        call(f"rm {file}")


def fq_to_fa(file: str, keep: bool = True):
    """
    Wrapper function of "seqtk" to convert fastq into fasta

    Args:
        file: path-like

        keep:
            Keep the input file or not
    """
    if file.endswith('.fastq'):
        output = file[:-6] + '.fa'
    elif file.endswith('.fq'):
        output = file[:-3] + '.fa'
    else:
        output = file + '.fa'

    call(f"seqtk seq -A {file} > {output}")
    if not keep:
        call(f"rm {file}")


def vcf_to_bcf(file: str, keep: bool = True):
    """
    Wrapper function of "bcftools view" to convert VCF into BCF

    Args:
        file: path-like

        keep:
            Keep the input file or not
    """
    # -Ou: output uncompressed bcf
    # -o <file_out>
    call(f"bcftools view -Ou -o {file[:-4]}.bcf {file}")
    if not keep:
        call(f"rm {file}")


def bcf_to_vcf(file: str, keep: bool = True):
    """
    Wrapper function of "bcftools view" to convert BCF into VCF

    Args:
        file: path-like

        keep:
            Keep the input file or not
    """
    # -Ov: output uncompressed vcf
    # -o <file_out>
    call(f"bcftools view -Ov -o {file[:-4]}.vcf {file}")
    if not keep:
        call(f"rm {file}")
