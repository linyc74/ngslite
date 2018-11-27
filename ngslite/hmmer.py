import subprocess


from .fasta import FastaParser
from .fasta import FastaWriter
from .dnatools import translate
from .dnatools import rev_comp


def __call(cmd):
    print('CMD: ' + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        print(inst)


def __gzip(file, keep=True):
    """
    Call the "gzip" command to zip or unzip files.

    Args:
        file: str, path-like

        keep: bool, keep the input file or not
    """
    decomp = ['', '-d '][file.endswith('.gz')]
    keep = ['', '-k '][keep]
    __call(f"gzip {decomp}{keep}{file}")


def __is_dna(fasta):
    """
    Returns 'True' if the first sequence of the <fasta> file is a DNA
    """
    parser = FastaParser(fasta)
    is_dna = True
    for char in set(parser.next()[1]):
        if not char in ('A', 'C', 'G', 'T'):
            is_dna = False
            break
    parser.close()
    return is_dna


def __translate_dna_database(fasta):
    """
    Translate the input fasta file in six frames.
    Append ';frame=<frame>' in the header line.
    Write a new '<fasta>_translated.fa' file.

    Args:
        fasta: str, path-like

    Returns: str,
        The written fasta name '<fasta>_translated.fa'
    """
    if fasta.endswith('.fa'): output = fasta[:-3] + '_translated.fa'
    elif fasta.endswith('.fasta'): output = fasta[:-6] + '_translated.fa'

    parser = FastaParser(fasta)
    writer = FastaWriter(output, 'w')
    while True:
        head, seq = parser.next()
        if head is None or seq is None:
            break

        for frame in (1, 2, 3):
            new_head = f"{head};frame={frame}"
            aa_seq = translate(seq[frame-1:])
            writer.write(new_head, aa_seq)

        for frame in (-1, -2, -3):
            new_head = f"{head};frame={frame}"
            rc_seq = rev_comp(seq)
            aa_seq = translate(rc_seq[-frame-1:])
            writer.write(new_head, aa_seq)

    return output


def hmmsearch(hmm, database, output):
    """
    Args:
        hmm: str, path-like
            The .hmm file built by the command "hmmbuild",
                e.g. "hmmbuild -o Pfam-A_summary.txt Pfam-A.hmm Pfam-A.seed"
            Also accepts .gz format

        database: str, path-like
            The fasta database to be searched against
            Also accepts .gz format

        output: str, path-like
            The output text file reported by the command "hmmsearch"
    """
    db_is_gz = False
    if database.endswith('.gz'):
        __gzip(database, keep=True)
        database = database[:-len('.gz')]
        db_is_gz = True

    hmm_is_gz = False
    if hmm.endswith('.gz'):
        __gzip(hmm, keep=True)
        hmm = hmm[:-len('.gz')]
        hmm_is_gz = True

    db_is_dna = False
    if __is_dna(database):
        database = __translate_dna_database(database)
        db_is_dna = True

    # Run hmmsearch
    __call(f"hmmsearch {hmm} {database} > {output}")

    if db_is_gz or db_is_dna:
        __call(f"rm {database}")
    if hmm_is_gz:
        __call(f"rm {hmm}")


def hmmbuild(seed, hmm):
    """
    Args:
        seed: str, path-like
            The input .seed file

        hmm: str, path-like
            The output .hmm file
    """
    if not hmm.endswith('.hmm'):
        hmm = hmm + '.hmm'
    summary = hmm[:-len('.hmm')] + '_summary.txt'

    __call(f"hmmbuild -o {summary} {hmm} {seed}")

