from .gtftools import GtfWriter
from .fasta import FastaParser, FastaWriter


import subprocess
from functools import partial
printf = partial(print, flush=True)


def __call(cmd):
    printf('CMD: ' + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        printf(inst)


def __remove_extension(file):
    """
    Args:
        file: str, path-like

    Returns: str,
        A path without the file extension
    """
    return '.'.join(file.split('.')[:-1])


def __parser_glimmer3_result(file, output):
    """
    Args:
        file: str, path-like
            The glimmer3 output prediction file, for example:

            >fasta_header
             orf00001      107      448  +2     3.51
             orf00002      687      800  +3     2.86
             orf00003      725      549  -3     5.63
                          (start)  (end)(frame)(score)

        output: str, path-like
            The output GTF file
    """
    with open(file) as parser:
        with GtfWriter(output) as writer:
            head = parser.readline().rstrip()[1:]
            for line in parser:
                if line.startswith('>'):
                    head = line.rstrip()[1:]
                else:
                    name, start, end, frame, score = line.rstrip().split()
                    if frame[0] == '-':
                        start, end = end, start
                    writer.write(
                        (head,              # seqname
                        '.',                # source
                        'CDS',              # feature
                        start,              # start
                        end,                # end
                        score,              # score
                        frame[0],           # strand
                        frame[1:],          # frame
                        f"name \"{name}\"") # attribute
                    )


def glimmer3(fasta, output, linear=True, max_overlap=50, min_length=110, threshold=30, entropy_distance=1.15, log=None):
    """
    Args:
        fasta: str, path-like
            The input genome sequences

        output: str, path-like
            The output GTF file containing predicted ORFs

        linear: bool
            Linear or circular DNA

        max_overlap: int
            Max overlap (bp) between ORFs

        min_length: int
            Minimum length of ORFs

        threshold: int

        entropy_distance: float
            Only genes with entropy distance score less than this value will be considered
            Default 1.15 in the command 'long-orfs'

        log: str, path-like
            The log file for stderr
    """
    f = fasta
    for ext in ['.fa', '.fna', '.fasta']:
        if fasta.endswith(ext):
            f = __remove_extension(fasta)

    if type(log) is str:
        log = f" 2>> {log}"
    else:
        log = f" 2>> {f}_glimmer3.log"


    with FastaParser(fasta) as parser:
        for head, seq in parser:

            with FastaWriter('temp.fa') as writer:
                writer.write(head, seq)

            # --- Find long ORFs --- #
            # cutoff: Only genes with entropy distance score less than <cutoff> will be considered. Default 1.15
            temp_fa = 'temp.fa'
            linear = ['--linear ', ''][linear]
            cmd = f"long-orfs --no_header --cutoff {entropy_distance} {linear}{temp_fa} longorfs{log}"
            __call(cmd)

            # --- Extract the DNA sequences of long ORFs --- #
            cmd = f"extract {temp_fa} longorfs > train"
            __call(cmd)

            # --- Build ICM --- #
            # r: Use the reverse of input strings to build the model
            cmd = f"build-icm -r icm < train"
            __call(cmd)

            # Prepare output name
            if output.endswith('.gtf'): out = __remove_extension(output)
            else: out = output

            # --- Predict genes --- #
            # max_olap: Set maximum overlap length to <n>. Overlaps this short or shorter are ignored
            # gene_len: Set minimum gene length to <n>
            # threshold: Set threshold score for calling as gene to n. If the in-frame score >= <n>,
            #   then the region is given a number and considered a potential gene
            cmd = f"glimmer3 --max_olap {max_overlap} --gene_len {min_length} --threshold {threshold} {temp_fa} icm {out}{log}"
            __call(cmd)

            __parser_glimmer3_result(file=f"{out}.predict", output=output)

            __call(f"rm longorfs train icm {out}.predict")

