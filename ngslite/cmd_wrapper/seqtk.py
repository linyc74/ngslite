from typing import Optional
from ..lowlevel import call


def fq_to_fa(
        file: str,
        keep: bool = True,
        output: Optional[str] = None) -> str:
    """
    Args:
        file:
            Fastq file path

        keep:
            Keep the input file or not

        output:
            Fasta file path
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
