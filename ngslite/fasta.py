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

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def __iter__(self):
        return self

    def __next__(self):
        r = self.next()
        if r:
            return r
        else:  # r is None
            raise StopIteration

    def next(self):
        """
        Returns: tuple
            The next read of the fasta file.
            If it reaches the end of the file, return None.
        """
        header = self.__fasta.readline().rstrip()[1:]
        if header == '':
            return None

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

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

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

        If no sequences from the fasta, return an empty list
    """
    with FastaParser(file) as parser:
        return [(head, seq) for head, seq in parser]


def subset_fasta(file, headers, output):
    """
    Args:
        file: str, path-like
            The input fasta file

        headers: str, or list of str
            Each str is a header to be included

        output: str, path-like
            The output fasta file
    """
    if isinstance(headers, str):
        headers = [headers]

    with FastaParser(file) as parser:
        with FastaWriter(output) as writer:
            for head, seq in parser:
                if head in headers:
                    writer.write(head, seq)

