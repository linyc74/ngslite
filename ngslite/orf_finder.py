from .fasta import *
from .gtftools import *
from .lowlevel import __call, __temp
import os
from functools import partial
printf = partial(print, flush=True)


def orf_finder(fasta, output, min_length=75):
    """
    Use the command line tool NCBI ORFfinder to find Open Reading Frames
      in the input fasta file, and output a GTF file

    Args:
        fasta: str, path-like

        output: str, path-like

        min_length: int
            Minimal length of ORF. 75 is the default of ORFfinder
    """
    # ORFfinder can only take fasta headers <= 50 characters
    # Create a temporary fasta file with short headers
    #   and a dictionary matching the short and the original headers
    temp_fa = temp('temp', '.fa')  # for example, temp000.fa
    contig_dict = {}
    with FastaParser(fasta) as parser:
        with FastaWriter(temp_fa) as writer:
            for i, contig in enumerate(parser):
                head, seq = contig
                writer.write(str(i), seq)
                contig_dict[str(i)] = head

    # -outfmt 3: output format = feature table
    # -ml [int]: minimal length
    temp_table = __temp('temp', '.table')  # for example, temp000.table
    __call(f"ORFfinder -in {temp_fa} -out {temp_table} -outfmt 3 -ml {min_length}")

    os.remove(temp_fa)

    with GtfWriter(output) as writer:
        with open(temp_table) as fh:
            while True:
                line1 = fh.readline().rstrip()
                line2 = fh.readline().rstrip()
                line3 = fh.readline().rstrip()
                if line1 == '':
                    break

                # line1: Get contig_number and orf_id
                contig_number = line1.split('_')[1].split(':')[0]
                orf_id = line1.split('>Feature lcl|')[1].split('_')[0]

                # Get the original contig header in the input fasta file
                header = contig_dict[contig_number]

                # line2: Get start and end (zero-based)
                start, end = line2.split('\t')[0:2]
                start, end = int(start), int(end)

                if start < end:
                    strand = '+'
                else:
                    strand = '-'
                    start, end = end, start

                # The start and end of GTF is 1-based
                feature = \
                    (header              ,  # 1 seqname   str
                     '.'                 ,  # 2 source    str
                     'CDS'               ,  # 3 feature   str
                     start + 1           ,  # 4 start     int
                     end + 1             ,  # 5 end       int
                     '.'                 ,  # 6 score     float
                     strand              ,  # 7 strand    str ('+', '-')
                     0                   ,  # 8 frame     int (0, 1, 2)
                     f"name \"{orf_id}\"")  # 9 attribute str
                writer.write(feature)
        os.remove(temp_table)

