from typing import List, Optional, Dict, Union
from .lowlevel import printf
from .dataclass import GenericFeature, GffFeature
from .feature_conversion import generic_to_gff_feature


class GffParser:
    """
    File parser for GFF3 format
    """
    def __init__(self, file: str):
        """
        Args:
            file: path-like
        """
        self.__gff = open(file, 'r')
        line1 = self.__gff.readline().rstrip()
        assert line1 == '##gff-version 3', 'The first line is not "##gff-version 3"'

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def __iter__(self):
        self.__gff.seek(0)  # 0 = first character
        self.__gff.readline()
        return self

    def __next__(self):
        r = self.next()
        if r:
            return r
        else:  # r is None
            raise StopIteration

    def next(self) -> Optional[GffFeature]:
        """
        Each line of the GFF3 file has 9 fields

        #   Col	Field        Type
        0   1   seqid        str
        1   2   source       str
        2   3   type         str
        3   4   start        int
        4   5   end          int
        5   6   score        float
        6   7   strand       str ('+', '-')
        7   8   phase        int (0, 1, 2)
        8   9   attributes   str

        Returns: namedtuple of str or int
            9 fields of a feature (i.e. a line)
        """
        line = self.__gff.readline().rstrip()
        if line:
            fields: List[Union[str, int, float]] = line.split('\t')
            for i in (3, 4, 7):
                if fields[i] != '.':
                    fields[i] = int(fields[i])
            if fields[5] != '.':
                fields[5] = float(fields[5])
            return GffFeature._make(fields)
        else:  # line == ''
            return None

    def close(self):
        self.__gff.close()


class GffWriter:
    def __init__(self, file: str, mode: str = 'w'):
        """
        Args:
            file: path-like

            mode:
                The file mode: 'w' or 'a'
        """
        self.__gff = open(file, mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def write(self, feature: GffFeature):
        """
        Args:
            feature:
                Containing 9 fields of a line of GFF3 file
        """
        self.__gff.write('\t'.join(map(str, feature)) + '\n')

    def close(self):
        self.__gff.close()


def read_gff(
        file: str,
        as_dict: bool = False) \
        -> Union[List[GffFeature], Dict[str, GffFeature]]:
    """
    Args:
        file: path-like
            The input GFF file

        as_dict:
            If True, returns a dictionary

    Returns: list of GffFeature objects, or dict

        [GffFeature_1, GffFeature_2, ...]

        or

        {
            seqid_1: [GffFeature_1, ...],
            seqid_2: [GffFeature_2, ...], ...
        }
    """
    with GffParser(file) as parser:
        features = [feature for feature in parser]
    if as_dict:
        D = {}
        for f in features:
            D.setdefault(f.seqid, [])
            D[f.seqid].append(f)
        return D
    return features


def write_gff(
        data: Union[List[Union[GffFeature, GenericFeature]],
                    Dict[str, Union[GffFeature, GenericFeature]]],
        file: str):
    """
    Take the data in the format returned by read_gff()
        and write it into a new GFF file

    Accepts GffFeature and GenericFeature objects

    Args:
        data: list of GffFeature objects, or dict

        file: path-like
            The output GFF file
    """
    with GffWriter(file) as writer:
        if type(data) is dict:
            for feature_arr in data.values():
                for feature in feature_arr:
                    if type(feature) is GenericFeature:
                        feature = generic_to_gff_feature(feature)
                    writer.write(feature)
        elif type(data) is list:
            for feature in data:
                if type(feature) is GenericFeature:
                    feature = generic_to_gff_feature(feature)
                writer.write(feature)


def subset_gff(
        file: str,
        seqid: Union[str, List[str]],
        output: str):
    """
    Args:
        file: path-like
            The input GFF file

        seqid:
            Each str is a seqid to be included

        output: path-like
            The output GFF file
    """
    if isinstance(seqid, str):
        seqid = [seqid]

    with GffParser(file) as parser:
        with GffWriter(output) as writer:
            for feature in parser:
                if feature[0] in seqid:
                    writer.write(feature)


def print_gff(
        feature: Optional[GffFeature] = None):
    """
    Pretty print a feature (namedtuple) from GFF or GFF file

    Args:
        feature:
            Containing 9 fields of a feature from GFF or GFF file
    """
    if feature is None:
        text = """\
#   Col	Field       Type
0   1   seqid     string
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
            'seqid    ',
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
