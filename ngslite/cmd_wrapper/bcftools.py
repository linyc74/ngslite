from typing import Optional, List
from ..lowlevel import call


def bcftools_variant_call(
        ref: str,
        bam: str,
        output: Optional[str] = None,
        max_depth: int = 250,
        haploid: bool = False):
    """
    Wrapper function of "bcftools mpileup" -> BCF (genotype likelihood) -> "bcftools call" -> BCF (variants)

    Args:
        ref: path-like
            The reference fasta file

        bam: path-like
            The alignment BAM file

        output: path-like
            The output BCF file

        max_depth:
            Max per-file sequencing depth; avoids excessive memory usage, default 250

        haploid:
            This is for the option "--ploidy" in "bcftools call"
            If True, set the option "--ploidy 1" for haploid (bacteria)
            If False, ignore the option "--ploidy", the default is diploid
            Command "bcftools call --ploidy ?" prints all available presets: GRCh37, GRCh38, X, Y, 1
    """
    # -Ou: output uncompressed bcf
    # -m: use the default calling method, i.e. alternative model for multiallelic and rare-variant calling
    # -v: output only variant sites
    # -o <file_out>
    # -f <ref_fasta>
    # -d <max_depth>
    # --ploidy <GRCh37|GRCh38|X|Y|1>
    if output is None:
        output = bam[:-4] + '.bcf'
    ploidy = ['', '--ploidy 1 '][haploid]
    call(f'bcftools mpileup -Ou -d {max_depth} -f {ref} {bam} | bcftools call -Ou -m -v {ploidy}-o {output}')


def vcf_to_bcf(
        file: str,
        keep: bool = True,
        output: Optional[str] = None) -> str:
    """
    Args:
        file:
            VCF file path

        keep:
            Keep the input file or not

        output:
            BCF file path
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
    Args:
        file:
            BCF file path

        keep:
            Keep the input file or not

        output:
            VCF file path

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


def sort_bcf(
        file: str,
        keep: bool = False):
    """
    Args:
        file: path-like

        keep:
            Keep the input file or not
    """
    file_out = f'{file[:-4]}_sorted.{file[-3:]}'
    # -Ou: output uncompressed bcf
    # -o <file_out>
    call(f'bcftools sort -Ou -o {file_out} {file}')
    if not keep:
        call(f'rm {file}')
        call(f'mv {file_out} {file}')


def subset_bcf_regions(
        file: str,
        regions: List[str],
        output: Optional[str] = None,
        keep: bool = True):
    """
    Args:
        file:
            Input BCF file path

        regions:
            Each str is a region of the reference genome, e.g.
                chr1            chromosome 1
                chr3:1000-2000  chromosome 3 from 1000th (inclusive) to 2000th (inclusive) base

        output:
            Output BCF file path
            If None, add subscript '_subset' to the input <file>

        keep: bool
            If False, delete the input <file> and rename the output as the input <file>
            Overrides the <output> file name
    """
    # Convert regions list into a string
    # ['chr1', 'chr2:1001-2000'] -> ' chr1,chr2:1001-2000'
    regions = ' ' + ','.join(regions)  # comma-separated regions

    if output is None:
        output = file[:-4] + '_subset.bcf'

    # -Ou: output uncompressed bcf
    # -o <file_out>
    call(f'bcftools view -Ou -o {output} {file}{regions}')

    if not keep:
        call(f'rm {file}')
        call(f'mv {output} {file}')
