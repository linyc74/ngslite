from .lowlevel import __call
import random
import os
from functools import partial
printf = partial(print, flush=True)


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


def __subsample_fq(file, fraction, output):
    """
    Randomly subsample the input fastq file

    Args:
        file: str, path-like
            The input fastq file

        fraction: float, fraction of reads or sequences to be retrieved

        output: str, path-like
            The output file name
    """
    fh_in = open(file, 'r')
    fh_out = open(output, 'w')

    # Count the total number of lines in the fastq file
    i = 0
    while fh_in.readline() != '':
        i += 1

    # Go back to the beginning of the file, because it needs to go through the file again
    fh_in.seek(0)

    N = int(i / 4)  # The total count of reads
    n = int(N * fraction)  # The count of reads to be retrieved

    # A random array of int corresponding to which reads to be retrieved
    rand_array = random.sample(list(range(N)), n)
    rand_array.sort()

    # Prepare the iterating variable pos (current line position) and next_pos
    pos = -4
    # 4 lines for each read in fastq, so if the next read is 10th read,
    #     then it starts at 40th line
    next_pos = rand_array.pop(0) * 4
    while len(rand_array) > 0:
        # Read 4 lines at a time
        four_lines = [fh_in.readline() for _ in range(4)]
        # Every time read 4 lines, update the current position pos
        pos += 4
        if pos == next_pos:
            # If pos is at the next position, write four lines
            for line in four_lines:
                fh_out.write(line)
            # Update the next position
            next_pos = rand_array.pop(0) * 4
        # If pos is not at the next position, then just continue onto the next read

    fh_in.close()
    fh_out.close()


def __subsample_fa(file, fraction, output):
    """
    Randomly subsample the input fasta file

    Args:
        file: str, path-like
            The input fasta file

        fraction: float, fraction of reads or sequences to be retrieved

        output: str, path-like
            The output file name
    """
    fh_in = open(file, 'r')
    fh_out = open(output, 'w')

    # All line positions of the header line '>'
    all_entry_pos = []
    for i, line in enumerate(fh_in):
        if line.startswith('>'):
            all_entry_pos.append(i)

    # Go back to the beginning of the file, because it needs to go through the file again
    fh_in.seek(0)

    N = len(all_entry_pos)  # Total number of reads
    n = int(N * fraction)  # The number of reads to be retrieved
    # Randomly select a subset of entry positions
    rand_entry_pos = random.sample(all_entry_pos, n)
    rand_entry_pos.sort()

    # Prepare the iterating variable line, i and next_pos
    line = fh_in.readline()
    pos = 0
    next_pos = rand_entry_pos.pop(0)
    while True:
        if pos != next_pos:
            # If the current position is not the next position to look for,
            #     then just move on to the next line
            line = fh_in.readline()
            pos += 1
        else:
            # pos == next_pos, hit the header line to look for!
            fh_out.write(line)  # Write the header line

            # The following while loop writes lines until the next header line '>' is seen
            while True:
                line = fh_in.readline()
                pos += 1
                if line.startswith('>') or line == '':
                    # Starts with '>', i.e. a header line, or reach the end of the input file
                    # Break out this while loop for writing the sequence
                    break
                else:
                    # The line is still before the next header line '>'
                    # So it must be part of the DNA sequence to be written into the output file
                    fh_out.write(line)

            # The sequence has been written into the output file
            if len(rand_entry_pos) > 0:
                # If there's more entry position to look for, then update the next_pos
                next_pos = rand_entry_pos.pop(0)
            else:  # len(rand_entry_pos) == 0, no more entry position to look for
                break  # Break out from the main while loop and return

    fh_in.close()
    fh_out.close()


def __subsample_fq_pair(file1, file2, fraction, output1, output2):
    """
    Randomly subsample the input fastq file pair

    Args:
        file1: str, path-like
            The input fastq file 1

        file2: str, path-like
            The input fastq file 2

        fraction: float, fraction of reads or sequences to be retrieved

        output1: str, path-like
            The output fastq file 1

        output2: str, path-like
            The output fastq file 2
    """

    fh_in1 = open(file1, 'r')
    fh_in2 = open(file2, 'r')
    fh_out1 = open(output1, 'w')
    fh_out2 = open(output2, 'w')

    # Count the total number of lines in the fastq file
    i = 0
    while fh_in1.readline() != '':
        i += 1

    # Go back to the beginning of the file, because it needs to go through the file again
    fh_in1.seek(0)

    N = int(i / 4)  # The total count of reads
    n = int(N * fraction)  # The count of reads to be retrieved

    # A random array of int corresponding to which reads to be retrieved
    rand_array = random.sample(list(range(N)), n)
    rand_array.sort()

    # Prepare the iterating variable pos (current line position) and next_pos
    pos = -4
    # 4 lines for each read in fastq, so if the next read is 10th read,
    #     then it starts at 40th line
    next_pos = rand_array.pop(0) * 4
    while len(rand_array) > 0:
        # Read 4 lines at a time
        lines1 = [fh_in1.readline() for _ in range(4)]  # .1.fq
        lines2 = [fh_in2.readline() for _ in range(4)]  # .2.fq
        # Every time read 4 lines, update the current position pos
        pos += 4
        if pos == next_pos:
            # If pos is at the next position, write four lines
            for l in lines1:
                fh_out1.write(l)  # .1.fq
            for l in lines2:
                fh_out2.write(l)  # .2.fq
            # Update the next position
            next_pos = rand_array.pop(0) * 4
        # If pos is not at the next position, then just continue onto the next read

    fh_in1.close()
    fh_in2.close()
    fh_out1.close()
    fh_out2.close()


def __subsample_sam(file, fraction, output):
    """
    Randomly subsample the input SAM or BAM file

    Args:
        file: str, path-like
            The input SAM or BAM file

        fraction: float, fraction of reads or sequences to be retrieved

        output: str, path-like
            The output file name
    """

    # The output file format (SAM or BAM) depends on the input format
    if file.endswith('.sam'): output_bam = ''
    elif file.endswith('.bam'): output_bam = '-b '

    __call(f"samtools view -s {fraction} {output_bam}{file} > {output}")


def subsample(file, fraction, file2='', output='', output2=''):
    """
    Randomly subsample the reads in the input <file>, and write them into the <output> file

    Supported file formats:
        Fastq (single or paired end)
        Fasta
        SAM/BAM

    The input file could be from other folder,
        but the output file is always written in the current folder.

    Supports compressed input files (automatically detects .gz file name).

    Args:
        file: str, path-like
            The input file (fastq, fasta, sam, bam)

        fraction: float
            Fraction of reads or sequences to be randomly sampled

        file2: str, path-like
            The filename of the second fastq file for paired-end

        output: str, path-like
            The output file name
            If '', then add the input <file> with prefix 'subset_' as the <output> file name

        output2: str, path-like
            The output file name for <file2> of paired-end fastq
            If '', then add the input <file2> with prefix 'subset_' as the <output2> file name
    """
    # Unzip the input file if endswith '.gz'
    if file.endswith('.gz'):
        __gzip(file, keep=True)
        file = file[:-3]  # Remove '.gz' suffix
        is_gz = True
    else:
        is_gz = False

    if file2.endswith('.gz'):
        __gzip(file2, keep=True)
        file2 = file2[:-3]  # Remove '.gz' suffix
        is_gz2 = True
    else:
        is_gz2 = False

    # Determine file type
    if file.endswith('.fa'):
        ftype = 'fa'
    elif (file.endswith('.fq') or file.endswith('.fastq'))\
            and file2 == '':
        ftype = 'fq'
    elif (file.endswith('.fq') or file.endswith('.fastq'))\
            and (file2.endswith('.fq') or file2.endswith('.fastq')):
        ftype = 'fq_paired'
    elif file.endswith('.sam') or file.endswith('.bam'):
        ftype = 'sam'
    else:
        ftype = ''

    if output == '':
        output = 'subset_' + file.split('/')[-1]  # Remove preceding folder path
    if output2 == '' and ftype == 'fq_paired':
        output2 = 'subset_' + file2.split('/')[-1]

    if ftype == 'fa': __subsample_fa(file, fraction, output)
    elif ftype == 'fq': __subsample_fq(file, fraction, output)
    elif ftype == 'fq_paired': __subsample_fq_pair(file, file2, fraction, output, output2)
    elif ftype == 'sam': __subsample_sam(file, fraction, output)

    # Remove the temporary unzipped file
    if is_gz: os.remove(file)
    if is_gz2: os.remove(file2)

