from typing import List, Union, Optional, Dict
from .lowlevel import printf
from .dataclass import GenericFeature, GtfFeature
from .feature_conversion import generic_to_gtf_feature


class GtfParser:
    def __init__(self, file: str):
        """
        Args:
            file: path-like
        """
        self.__gtf = open(file, 'r')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def __iter__(self):
        self.__gtf.seek(0)
        return self

    def __next__(self):
        r = self.next()
        if r:
            return r
        else:  # r is None
            raise StopIteration

    def next(self) -> Optional[GtfFeature]:
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
            fields: List[Union[str, int, float]] = line.split('\t')
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
    def __init__(self, file: str, mode: str = 'w'):
        """
        Args:
            file: path-like

            mode:
                The file mode: 'w' or 'a'
        """
        self.__gtf = open(file, mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def write(self, feature: GtfFeature):
        """
        Args:
            feature: namedtuple 'GtfFeature'
                Containing 9 fields of a line of GTF file
        """
        self.__gtf.write('\t'.join(map(str, feature)) + '\n')

    def close(self):
        self.__gtf.close()


def read_gtf(
        file: str,
        as_dict: bool = False) \
        -> Union[List[GtfFeature], Dict[str, GtfFeature]]:
    """
    Args:
        file: path-like object
            The input GTF file

        as_dict: bool
            If True, returns a dictionary

    Returns: list of GtfFeature objects, or dict

        [GtfFeature_1, GtfFeature_2, ...]

        or

        {
            seqname_1: [GtfFeature_1, ...],
            seqname_2: [GtfFeature_2, ...], ...
        }
    """
    with GtfParser(file) as parser:
        features = [feature for feature in parser]
    if as_dict:
        D = {}
        for f in features:
            D.setdefault(f.seqname, [])
            D[f.seqname].append(f)
        return D
    return features


def write_gtf(
        data: Union[List[Union[GtfFeature, GenericFeature]],
                    Dict[str, Union[GtfFeature, GenericFeature]]],
        file: str):
    """
    Take the data in the format returned by read_gtf()
        and write it into a new GTF file

    Accepts GtfFeature and GenericFeature objects

    Args:
        data: list of GtfFeature objects, or dict

        file: path-like
            The output GTF file
    """
    with GtfWriter(file) as writer:
        if type(data) is dict:
            for feature_arr in data.values():
                for feature in feature_arr:
                    if type(feature) is GenericFeature:
                        feature = generic_to_gtf_feature(feature)
                    writer.write(feature)
        elif type(data) is list:
            for feature in data:
                if type(feature) is GenericFeature:
                    feature = generic_to_gtf_feature(feature)
                writer.write(feature)


def subset_gtf(
        file: str,
        seqname: Union[str, List[str]],
        output: str):
    """
    Args:
        file: path-like
            The input GTF file

        seqname:
            Each str is a seqname (chromosome name) to be included

        output: path-like
            The output GTF file
    """
    if isinstance(seqname, str):
        seqname = [seqname]

    with GtfParser(file) as parser:
        with GtfWriter(output) as writer:
            for feature in parser:
                if feature[0] in seqname:
                    writer.write(feature)


def print_gtf(feature: Optional[GtfFeature] = None):
    """
    Pretty print a feature (namedtuple) from GTF or GFF file

    Args:
        feature:
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
        fields = [
            'seqname  ',
            'source   ',
            'feature  ',
            'start    ',
            'end      ',
            'score    ',
            'strand   ',
            'frame    ',
            'attribute',
        ]
        for i in range(9):
            printf(f"{i}\t{fields[i]}\t{feature[i]}")
