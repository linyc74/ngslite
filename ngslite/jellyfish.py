from .lowlevel import __call, __temp
from .file_conversion import fq_to_fa
from functools import partial
printf = partial(print, flush=True)


def __gzip(file, keep=True):
    """
    Call the "gzip" command to zip or unzip files.

    Args:
        file: str, path-like

        keep: bool, keep the input file or not
    """
    keep = ['', '-k '][keep]
    decomp = ['', '-d '][file.endswith('.gz')]
    __call(f"gzip {decomp}{keep}{file}")


def jellyfish_count(file, k, output, min_count=1, hash_size='100M', threads=4, canonical=True):
    """
    Wrapper function for the jellyfish k-mer count commands.

    From a high-level it takes an input fastq or fastq <file>, and eventually
        exports an <output> fasta file of k-mer counts.

    All the intermediate files are deleted.

    Args:
        file: str, path-like
            The input fasta or fastq file, also accepts .gz compressed format

        k: int
            The size of the k-mer

        output: str, path-like
            The output fasta file, which cannot be the same as the input <file>

        min_count: uint
            K-mers with counts < min_count will be excluded

        hash_size: str,
            e.g. '100M' for 100 million; '1G' for 1 billion

        threads: int (uint32)

        canonical: bool
            If True, 'AAAC' will be the same as 'GTTT'
    """
    # The <output> fasta should NOT be identical to the input fasta <file>
    if output == file: return

    input_name = file  # The original input name

    # Unzip
    if file.endswith('.gz'):
        __gzip(file, keep=True)
        file = file[:-3]

    # Convert fastq to fasta
    if file.endswith('.fastq') or file.endswith('.fq'):
        fq_to_fa(file, keep=True)
        # If the input <file> is not equal to the original input name, then don't keep it
        if not file == input_name: __call(f"rm {file}")

    if file.endswith('.fastq'): file = file[:-6] + '.fa'
    elif file.endswith('.fq'): file = file[:-3] + '.fa'

    # Count k-mers in the input fastq or fasta file
    canonical = ['', '-C '][canonical]
    temp_jf = __temp('temp', '.jf')  # for example, temp000.jf
    __call(f"jellyfish count -m {k} -s {hash_size} -t {threads} -o {temp_jf} {canonical}{file}")

    # Remove the counted file if it's not the original input file
    if not file == input_name: __call(f"rm {file}")

    # Dump the .jf file into a .fa file
    __call(f"jellyfish dump -L {min_count} --output {output} {temp_jf}")

    # The .jf file is not needed, so remove it.
    __call(f"rm {temp_jf}")
