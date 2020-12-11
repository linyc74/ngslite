from typing import Optional
from ..lowlevel import call


def muscle(
        fasta: str,
        output: str,
        maxiters: int = 16,
        outformat: str = 'fasta',
        log: Optional[str] = None):
    """
    Wrapper function for MSUCLE sequence alignment

    Args:
        fasta: path-like
            The input fasta file

        output: path-like
            The output file

        maxiters:
            Maximum number of iterations (default 16)

        outformat:
            'fasta' (default)
            'html'
            'msf' (GCF MSF)
            'clw' (ClustalW)

        log: path-like
            The log file, default None -> <output>.log
    """
    outformat = {'fasta': '', 'html': '-html', 'msf': '-msf', 'clw': '-clw'}[outformat]

    if log is None:
        log = output + '.log'

    call(f'muscle -in {fasta} -out {output} -maxiters {maxiters} -log {log} {outformat}')