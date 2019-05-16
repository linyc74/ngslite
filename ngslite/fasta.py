def _wrap(seq, width):
    """
    Args:
        seq: str
            A single line of DNA or protein sequence without '\n'

    Returns: str
    """
    if len(seq) <= width: return seq  # no need to wrap
    w = width
    L = []
    for i in range(int(len(seq)/w) + 1):
        L.append(seq[i*w:(i+1)*w])
    return '\n'.join(L)


class FastaParser(object):
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
        self.__fasta.seek(0)
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


class FastaWriter(object):
    """
    A simple fasta writer that writes a single read (header and sequence) each time
    """
    def __init__(self, file, mode='w'):
        """
        Args
            file: str, path-like

            mode: str
                'w' for write or 'a' for append
        """
        self.__fasta = open(file, mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def write(self, header, sequence, wrap=80):
        """
        Args:
            header: str

            sequence: str

            wrap: int
        """
        self.__fasta.write('>' + header + '\n' + _wrap(sequence, wrap) + '\n')

    def close(self):
        self.__fasta.close()


def read_fasta(file, as_dict=False):
    """
    Args:
        file: str, path-like object
            The input fasta file

        as_dict: bool
            If True, returns a dictionary

    Returns: list of tuples, or dict

        [(head_1, seq_1), (head_2, seq_2), ...]

        or

        {head_1: seq_1, head_2: seq_2, ...}

        If no sequences from the fasta, return an empty list or dict
    """
    with FastaParser(file) as parser:
        if as_dict:
            return {head: seq for head, seq in parser}
        else:
            return [(head, seq) for head, seq in parser]


def write_fasta(data, file):
    """
    Take the data in the format returned by read_fasta()
        and write it into a new fasta file

    Args:
        data: list of tuples, or dict

        file: str, path-like
            The output fasta file
    """
    with FastaWriter(file) as writer:
        if type(data) is dict:
            data = data.items()
        for head, seq in data:
            writer.write(head, seq)


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

