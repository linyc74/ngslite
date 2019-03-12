from functools import partial
printf = partial(print, flush=True)


class GtfParser:
    def __init__(self, file):
        """
        Args:
            file: str, path-like object
        """
        self.__gtf = open(file, 'r')

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
        Each line of the GTF file has 9 fields

        #   Col	Field       Type
        0   1   seqname     string
        1   2   source      string
        2   3   feature     string
        3   4   start       int
        4   5   end         int
        5   6   score       float
        6   7   strand      string ('+', '-')
        7   8   frame       int (0, 1, 2)
        8   9   attribute   string

        Returns: tuple of str or int
            9 fields of a feature (i.e. a line)
        """
        line = self.__gtf.readline().rstrip()
        if line:
            fields = line.split('\t')
            for i in (3, 4, 7):
                if fields[i] != '.':
                    fields[i] = int(fields[i])
            if fields[5] != '.':
                fields[5] = float(fields[5])
            return tuple(fields)
        else:  # line == ''
            return None

    def close(self):
        self.__gtf.close()


class GtfWriter:
    def __init__(self, file, mode='w'):
        """
        Args:
            file: str, path-like object
            mode: str, 'w' or 'a'
        """
        self.__gtf = open(file, mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def write(self, feature):
        """
        Args:
            feature: tuple of str or int
                Containing 9 fields of a line of GTF file
        """
        self.__gtf.write('\t'.join(map(str, feature)) + '\n')

    def close(self):
        self.__gtf.close()


def subset_gtf(file, seqname, output):
    """
    Args:
        file: str, path-like
            The input gtf file

        seqname: str, or list of str
            Each str is a seqname (chromosome name) to be included

        output: str, path-like
            The output gtf file
    """
    if isinstance(seqname, str):
        seqname = [seqname]

    with GtfParser(file) as parser:
        with GtfWriter(output) as writer:
            for feature in parser:
                if feature[0] in seqname:
                    writer.write(feature)


def print_gtf(feature=None):
    """
    Pretty print a feature (tuple) from GTF or GFF file

    Args:
        feature: tuple or list
            Containing 9 fields of a feature from GTF or GFF file
    """
    if feature is None:
        text = """\
#   Col	Field       Type
0   1   seqname     string
1   2   source      string
2   3   feature     string
3   4   start       int
4   5   end         int
5   6   score       float
6   7   strand      string ('+', '-')
7   8   frame       int (0, 1, 2)
8   9   attribute   string"""
        printf(text)

    elif type(feature) is tuple or type(feature) is list:
        fields = ['seqname  ', 'source   ', 'feature  ', 'start    ', 'end      ',
                  'score    ', 'strand   ', 'frame    ', 'attribute']
        for i in range(9):
            printf(f"{i}\t{fields[i]}\t{feature[i]}")

