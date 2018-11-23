def count_reads(file):
    """
    Count the reads in a fasta or fastq file.

    Args:
        file: str
            The file of the fasta or fastq file: ('.fq', '.fastq', '.fa', '.fasta')

    Returns:
        int: The number of reads or sequences contained in the input file
    """
    assert file.split('.')[-1] in ['.fq', '.fastq', '.fa', '.fasta'],\
        'File "{}" is neither a fasta nor a fastq.'.format(file)

    with open(file, 'r') as fh:

        # Use number of lines to count the reads in fastq file
        if file.endswith('.fq') or file.endswith('.fastq'):
            i = 0
            while fh.readline() != '':
                i += 1
            if i % 4 != 0:
                print('Warning: Number of lines in the fastq file is not multiple of 4.')
            count = int(i / 4)

        # Use '>' to count the reads in fasta file
        else:
            i = 0
            while True:
                line = fh.readline()
                if line == '':
                    break
                if line.startswith('>'):
                    i += 1
            count = i

    return count

