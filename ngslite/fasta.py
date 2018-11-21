class FastaParser:
    """
    A simple fasta parser that parses each read of a fasta file
    """
    def __init__(self, file):
        """
        Args:
            file: str, path-like object
        """
        self.__fasta = open(file, 'r')

    def next(self):
        """
        Returns: tuple
            The next read of the fasta file.
            If it reaches the end of the file, return None.
        """
        header = self.__fasta.readline().rstrip()[1:]
        if header == '':
            return (None, ) * 2

        seq = ''
        while True:
            pos = self.__fasta.tell()
            line = self.__fasta.readline().rstrip()
            if line.startswith('>'):
                self.__fasta.seek(pos)
                return header, seq
            if line == '':
                return header, seq
            seq = seq + line

    def close(self):
        self.__fasta.close()


class FastaWriter:
    """
    A simple fasta writer that writes a single read (header and sequence) each time
    """
    def __init__(self, file, mode='w'):
        """
        Args
            file: str, path-like object
            mode: str, 'w' for write or 'a' for append
        """
        self.__fasta = open(file, mode)

    def write(self, header, sequence):
        """
        Args:
            header: str
            sequence: str, i.e. DNA sequence
        """
        self.__fasta.write('>' + header + '\n' + sequence + '\n')

    def close(self):
        self.__fasta.close()


def read_fasta(file):
    """
    Args:
        file: str, path-like object
            The input fasta file

    Returns: list of tuples
        [(header_1, sequence_1), (header_2, sequence_2), ...]
    """
    with open(file, 'r') as fh:
        ret_list = []
        seq = ''
        # Get the first header
        head = fh.readline().rstrip()[1:]
        while True:
            line = fh.readline().rstrip()

            # If the line is a new header
            if line.startswith('>'):
                # Append the old header and the sequence accumulated so far
                ret_list.append((head, seq))
                # Reset the new header and clear the sequence (to be accumulated again)
                head = line[1:]
                seq = ''

            # If the end of the file
            elif line == '':
                # Append the last header and the sequence of the fasta
                ret_list.append((head, seq))
                break

            # Not a header, so just concatenate (i.e. accumulate) the sequence
            else:
                seq = seq + line
    return ret_list


def fasta_replace_blank_with(file, s, output):
    """
    Replace blank spaces with <s> in the headers of the input fasta file.

    Args:
        file: str, path-like
            The input fasta file

        s: str
            The string used to replace blank spaces in the headers

        output:
            The output fasta file
    """
    parser = FastaParser(file)
    writer = FastaWriter(output)

    while True:
        head, seq = parser.next()
        if head is None:
            break

        head = s.join(head.split(' '))
        writer.write(head, seq)

    parser.close()
    writer.close()

