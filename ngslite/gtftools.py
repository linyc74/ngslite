from functools import partial
from collections import namedtuple
printf = partial(print, flush=True)


GtfFeature = namedtuple('GtfFeature', 'seqname source feature start end score strand frame attribute')


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
        0   1   seqname     str
        1   2   source      str
        2   3   feature     str
        3   4   start       int
        4   5   end         int
        5   6   score       float
        6   7   strand      str ('+', '-')
        7   8   frame       int (0, 1, 2)
        8   9   attribute   str

        Returns: namedtuple of str or int
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
            return GtfFeature._make(fields)
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
            feature: namedtuple of str or int
                Containing 9 fields of a line of GTF file
        """
        self.__gtf.write('\t'.join(map(str, feature)) + '\n')

    def close(self):
        self.__gtf.close()


def subset_gtf(file, seqname, output):
    """
    Args:
        file: str, path-like
            The input GTF file

        seqname: str, or list of str
            Each str is a seqname (chromosome name) to be included

        output: str, path-like
            The output GTF file
    """
    if isinstance(seqname, str):
        seqname = [seqname]

    with GtfParser(file) as parser:
        with GtfWriter(output) as writer:
            for feature in parser:
                if feature[0] in seqname:
                    writer.write(feature)


def gtf_to_dict(file):
    """
    Args:
        file: str, path-like
            The input GTF file

    Returns: dict
        {
            seqname: [GtfFeature(), GtfFeature(), ...],
            seqname: [GtfFeature(), GtfFeature(), ...],
        }
    """
    D = dict()
    with GtfParser(file) as parser:
        for feature in parser:
            seqname = feature.seqname
            D.setdefault(seqname, [])
            D[seqname].append(feature)
    return D


def dict_to_gtf(dict_, output):
    """
    Args:
        dict_: dict
            The dictionary returned by gtf_to_dict()
            {
                seqname: [GtfFeature(), GtfFeature(), ...],
                seqname: [GtfFeature(), GtfFeature(), ...],
            }

        output: str, path-like
            The output GTF file
    """
    with GtfWriter(output) as writer:
        for features in dict_.values():
            for f in features:
                writer.write(f)


def attribute_str_to_dict(str_):
    """
    Args:
        str_: str
            name "2OG-FeII_Oxy_3  [M=96]";accession "PF13640.6";description "2OG-Fe(II) oxygenase superfamily";E_value "7e-19"

    Returns: dict
        {
            'name': '2OG-FeII_Oxy_3  [M=96]',
            'accession': 'PF13640.6'
            'description': '2OG-Fe(II) oxygenase superfamily'
            'E_value': '7e-19'
        }
    """
    D = dict()
    for a in str_.split(';'):
        key = a.split(' "')[0]
        val = a[len(key)+2:-1]
        D[key] = val
    return D


def attribute_dict_to_str(dict_):
    """
    Args:
        dict_: dict
            {
                'name': '2OG-FeII_Oxy_3  [M=96]',
                'accession': 'PF13640.6'
                'description': '2OG-Fe(II) oxygenase superfamily'
                'E_value': '7e-19'
            }

    Returns: str
        name "2OG-FeII_Oxy_3  [M=96]";accession "PF13640.6";description "2OG-Fe(II) oxygenase superfamily";E_value "7e-19"
    """
    s = ''
    for key, val in dict_.items():
        s = s + f"{key} \"{val}\";"
    return s[:-1]


def print_gtf(feature=None):
    """
    Pretty print a feature (namedtuple) from GTF or GFF file

    Args:
        feature: namedtuple, tuple or list
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

    else:
        fields = ['seqname  ', 'source   ', 'feature  ', 'start    ', 'end      ',
                  'score    ', 'strand   ', 'frame    ', 'attribute']
        for i in range(9):
            printf(f"{i}\t{fields[i]}\t{feature[i]}")
