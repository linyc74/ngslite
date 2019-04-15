from .gtftools import GtfWriter
from .fasta import FastaParser, FastaWriter
from .lowlevel import __call, __temp
import os
from functools import partial
printf = partial(print, flush=True)


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


def glimmer3(fasta, output, linear=True, max_overlap=50, min_length=110,
             threshold=30, entropy_distance=1.15, verbose=False, log=None):
    """
    This is a wrapper function of glimmer3, which predicts ORFs.

    Glimmer has to run separately on each contig.

    For each contig, this function does the following things:
        Find all possible ORFs of the contig    (long-orfs)
        Extract the DNA sequences of long ORFs  (extract)
        Build ICM                               (build-icm)
        Predict ORFs with the ICM               (glimmer3)

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

        verbose: bool,
            If True, print every single command line
            If False, only print the command line for the first sequence of the input <fasta>

        log: str, path-like
            The log file for stderr
    """
    f = fasta
    for ext in ['.fa', '.fna', '.fasta']:
        if fasta.endswith(ext):
            f = __remove_extension(fasta)

    if log is None: log = f"{f}_glimmer3.log"
    log = f" 2>> {log}"

    if os.path.isfile(output): os.remove(output)
    output_gtf = open(output, 'a')

    with FastaParser(fasta) as parser:
        i = -1
        for head, seq in parser:
            i += 1
            if verbose or i == 0:
                print_cmd = True
                printf(f"contig: {head}")
            else:
                print_cmd = False

            temp = __temp(prefix='temp', suffix='')  # for example, temp000
            with FastaWriter(f"{temp}.fa") as writer:
                writer.write(head, seq)

            # --- Find long ORFs --- #
            # cutoff: Only genes with entropy distance score less than <cutoff> will be considered. Default 1.15
            line = ['', '--linear '][linear]
            cmd = f"long-orfs --no_header --cutoff {entropy_distance} {line}{temp}.fa longorfs{log}"
            __call(cmd, print_cmd)

            # --- Extract the DNA sequences of long ORFs --- #
            cmd = f"extract {temp}.fa longorfs > train"
            __call(cmd, print_cmd)

            # --- Build ICM --- #
            # r: Use the reverse of input strings to build the model
            cmd = f"build-icm -r icm < train"
            __call(cmd, print_cmd)

            # --- Predict genes --- #
            # max_olap: Set maximum overlap length to <n>. Overlaps this short or shorter are ignored
            # gene_len: Set minimum gene length to <n>
            # threshold: Set threshold score for calling as gene to n. If the in-frame score >= <n>,
            #   then the region is given a number and considered a potential gene
            cmd = f"glimmer3 --max_olap {max_overlap} --gene_len {min_length} --threshold {threshold} {temp}.fa icm {temp}{log}"
            __call(cmd, print_cmd)
            # glimmer3 outputs two files: {temp}.detail and {temp}.predict

            __parser_glimmer3_result(file=f"{temp}.predict", output=f"{temp}.gtf")

            # Append the temp gtf to the output gtf
            with open(f"{temp}.gtf") as fh:
                output_gtf.write(fh.read())

            __call(f"rm {temp}.fa longorfs train icm {temp}.detail {temp}.predict {temp}.gtf", print_cmd)

        output_gtf.close()

