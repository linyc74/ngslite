import subprocess
import random
import os


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


def __subset_fq(file, fraction):
    """
    A low-level method for the public method random_subset()

    Args:
        file: str, path-like

        fraction: float, fraction of reads or sequences to be retrieved
    """
    fh_in = open(file, 'r')
    fh_out = open('subset_'+file.split('/')[-1], 'w')

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


def __subset_fa(file, fraction):
    """
    Args:
        file: str, path-like
            The input fasta file

        fraction: float, fraction of reads or sequences to be retrieved
    """
    fh_in = open(file, 'r')
    fh_out = open('subset_'+file.split('/')[-1], 'w')

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


def __subset_fq_pair(file1, file2, fraction):
    """
    A low-level method for the public method random_subset()

    Args:
        file1: str, path-like

        file2: str, path-like

        fraction: float, fraction of reads or sequences to be retrieved
    """

    fh_in1 = open(file1, 'r')
    fh_in2 = open(file2, 'r')
    fh_out1 = open('subset_'+file1.split('/')[-1], 'w')
    fh_out2 = open('subset_'+file2.split('/')[-1], 'w')

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


def random_subset(file, fraction, file2=None):
    """
    Write a file containing a random subset of the input file(s).
    The output file adds a prefix 'subset_' to the input file(s).
    Supports paired-end reads.

    Args:
        file: str, path-like object
            The input fasta or fastq file

        fraction: float
            Fraction of reads or sequences to be randomly sampled

        file2: str, path-like object
            The filename of the second fastq file for paired-end
    """
    # Unzip the input file (keep original)
    is_file_gz = False
    if file.endswith('.gz'):
        __gzip(file)
        file = file[:-3]  # Remove '.gz' suffix
        is_file_gz = True

    # Unzip input file2 (keep original), if files is not None and endswith '.gz'
    is_file2_gz = False
    if not file2 is None:
        if file2.endswith('.gz'):
            __gzip(file2)
            file2 = file2[:-3]  # Remove '.gz' suffix
            is_file2_gz = True

    # Simply file extensions
    if file.endswith('.fasta'):
        file = file[:-6] + '.fa'
    elif file.endswith('.fastq'):
        file = file[:-6] + '.fq'
    if not file2 is None:
        if file2.endswith('.fastq'):
            file2 = file2[:-6] + '.fq'

    # There could only be three cases
    # 1) Single fasta
    if file.endswith('.fa'):
        __subset_fa(file, fraction)
    # 2) Single fastq
    elif file.endswith('.fq') and file2 is None:
        __subset_fq(file, fraction)
    # #) Paired fastq
    elif file.endswith('.fq') and file2.endswith('.fq'):
        __subset_fq_pair(file, file2, fraction)

    # Remove the temporary unzipped file
    if is_file_gz:
        os.remove(file)
    if is_file2_gz:
        os.remove(file2)

