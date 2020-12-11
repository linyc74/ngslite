from typing import Optional, List, Tuple, Union, Dict
from .lowlevel import printf


class VcfParser:
    """
    A VCF file parser for VCFv4.3
    """
    def __init__(self, file: str):
        """
        The header section is parsed upon instantiation

        Args:
            file: path-like
        """
        self.__vcf = open(file, 'r')
        header = ''
        while True:
            pos = self.__vcf.tell()
            line = self.__vcf.readline()
            if line.startswith('#'):
                header = header + line
            else:
                self.__vcf.seek(pos)
                self.pos_0 = pos
                break
        self.header = header.rstrip()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def __iter__(self):
        self.__vcf.seek(self.pos_0)
        return self

    def __next__(self):
        r = self.next()
        if r:
            return r
        else:  # r is None
            self.__vcf.close()
            raise StopIteration

    def next(self) -> Optional[Tuple[Union[str, int]]]:
        """
        Each line of VCF file has 8 fixed, mandatory fields

        #   Col	Field   Type    Description
        0   1   CHROM   String  The name of the sequence (typically a chromosome)
        1   2   POS     Int     The 1-based position of the variation on the given sequence
        2   3   ID      String  The identifier of the variation
        3   4   REF     String  The reference base (or bases in the case of an indel) at the given position on the given reference sequence
        4   5   ALT     String  The list (separated by ,) of alternative alleles at this position
        5   6   QUAL    Int     A quality score associated with the inference of the given alleles.
        6   7   FILTER  String  A flag indicating which of a given set of filters the variation has passed
        7   8   INFO    String  An extensible list of key-value pairs (fields) describing the variation, for example, NS=2;DP=10;AF=0.333,0.667;AA=T;DB

        Returns: tuple of str or int
            8 fields of a variant (i.e. a line) plus optional fields
        """
        line = self.__vcf.readline().rstrip()
        if line:
            fields: List[Union[str, int]] = line.split('\t')
            for i in [1, 5]:
                fields[i] = int(fields[i])
            return tuple(fields)
        else:  # line == ''
            return None

    def close(self):
        self.__vcf.close()


class VcfWriter:
    def __init__(self, file: str, header: str, mode: str = 'w'):
        """
        The header section should be written in upon instantiation

        Args:
            file: path-like

            header:
                The header section

            mode:
                File mode: 'w' or 'a'
        """
        self.__vcf = open(file, mode)
        if mode == 'w':
            if not header == '':
                self.__vcf.write(header.rstrip() + '\n')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def write(self, variant: Tuple[Union[str, int]]):
        """
        Args:
            variant: tuple
                Containing at least 8 fields of a line of VCF file
        """
        self.__vcf.write('\t'.join(map(str, variant)) + '\n')

    def close(self):
        self.__vcf.close()


def print_vcf(var: Tuple[Union[str, int]] = None):
    """
    Pretty print a variant (tuple) from vcf file

    Args:
        var: tuple or list of str or int
            Containing (at least) 8 fields of a variant from vcf file
    """
    if var is None:
        text = ''''\
#   Col	Field   Type    Description
0   1   CHROM   String  The name of the sequence (typically a chromosome)
1   2   POS     Int     The 1-based position of the variation on the given sequence
2   3   ID      String  The identifier of the variation
3   4   REF     String  The reference base (or bases in the case of an indel) at the given position on the given reference sequence
4   5   ALT     String  The list (separated by ,) of alternative alleles at this position
5   6   QUAL    Int     A quality score associated with the inference of the given alleles.
6   7   FILTER  String  A flag indicating which of a given set of filters the variation has passed
7   8   INFO    String  An extensible list of key-value pairs (fields) describing the variation, for example, NS=2;DP=10;AF=0.333,0.667;AA=T;DB
8   9                   Optional fields...'''
        printf(text)

    elif type(var) is tuple or type(var) is list:
        fields = ['CHROM ', 'POS   ', 'ID    ', 'REF   ', 'ALT   ',
                  'QUAL  ', 'FILTER', 'INFO  ']
        for i in range(8):
            printf(f"{i}\t{fields[i]}\t{var[i]}")
        if len(var) > 8:
            for i in range(8, len(var)):
                printf(f"{i}\t     \t{var[i]}")


def unpack_vcf_info(var: Tuple[Union[str, int]]) -> Dict[str, str]:
    """
    Args:
        var: tuple or list of str or int
            Containing (at least) 8 fields of a variant from vcf file

    Return: dict
        Key, value pairs from the INFO (8th column) of a variant
        Values still str, not converted to int or float
    """
    def func(x: str) -> List[str]:
        return x.split('=') if '=' in x else (x, None)
    return {key: val for key, val in map(func, var[7].split(';'))}
