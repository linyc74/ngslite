from typing import Optional, Tuple, Union, Dict, List


class FastaParser:

    def __init__(self, file: str):
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

    def next(self) -> Optional[Tuple[str, str]]:
        """
        Returns the next read of the fasta file
        If it reaches the end of the file, return None
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

    def __init__(self, file: str, mode: str = 'w'):
        """
        Args:
            file: path-like

            mode: 'w' for write or 'a' for append
        """
        assert mode in ['w', 'a']

        self.__fasta = open(file, mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __wrap(self, seq: str, length: int) -> str:
        """
        Wraps a single line of string by a specified length

        Args:
            seq: a single line of DNA or protein sequence without '\n'

            length: length of each line
        """
        if len(seq) <= length:
            return seq  # no need to wrap

        w = length
        list_ = []
        for i in range(int(len(seq)/w) + 1):
            list_.append(seq[i*w:(i+1)*w])

        return '\n'.join(list_)

    def write(self, header: str, sequence: str, wrap: int = 80):
        """
        Args:
            header: Fasta header

            sequence: DNA or protein sequence without '\n'

            wrap: length of each wrapped line for the sequence
        """
        seq = self.__wrap(sequence, wrap)
        self.__fasta.write(f'>{header}\n{seq}\n')

    def close(self):
        self.__fasta.close()


def read_fasta(
        file: str,
        as_dict: bool = False,
        strip_header: bool = False) \
        -> Union[Dict[str, str], List[Tuple[str, str]]]:
    """
    Args:
        file: path-like

        as_dict:
            Return dictionary or not

        strip_header:
            Strip everything after the first space in the header

    Returns:
        {
            header: sequence, ...
            header: sequence, ...
        }

        [
            (header, sequence), (header, sequence), ...
        ]
    """

    parser = FastaParser(file)

    if as_dict:
        if strip_header:
            data = {h.split(' ')[0]: s for h, s in parser}
        else:
            data = {h: s for h, s in parser}
    else:
        if strip_header:
            data = [(h.split(' ')[0], s) for h, s in parser]
        else:
            data = [(h, s) for h, s in parser]

    parser.close()

    return data


def write_fasta(
        data: Union[Dict[str, str], List[Tuple[str, str]]],
        file: str):
    """
    Args:
        data:
            Fasta data return by the function 'read_fasta'

        file: path-like
    """

    if type(data) is list:
        iterator = data
    else:
        iterator = data.items()

    with FastaWriter(file) as writer:
        for head, seq in iterator:
            writer.write(header=head, sequence=seq)


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
