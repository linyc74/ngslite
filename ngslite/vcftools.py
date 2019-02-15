from functools import partial
printf = partial(print, flush=True)


class VcfParser:
    """
    A VCF file parser for VCFv4.3
    """
    def __init__(self, file):
        """
        The header section is parsed upon instantiation.

        Args:
            file: str, path-like object
        """
        self.__vcf = open(file, 'r')
        header = ''
        while True:
            # Get the current position
            pos = self.__vcf.tell()
            # Readline and move on to the next position
            line = self.__vcf.readline()
            if line.startswith('#'):
                header = header + line
            else:
                # If reaching the alignment section, that is,
                #   the header section has been parsed completely,
                #   then go back one line and break out the loop
                self.__vcf.seek(pos)
                break
        # Store the header string
        self.__header = header

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
            self.__vcf.close()
            raise StopIteration

    def next(self):
        """
        Each line of the SAM file has at least 11 fields

        #   Col	Field   Type    Description
        0   1   CHROM   String  The name of the sequence (typically a chromosome)
        1   2   POS     Int     The 1-based position of the variation on the given sequence
        2   3   ID      String  The identifier of the variation
        3   4   REF     String  The reference base (or bases in the case of an indel) at the given position on the given reference sequence
        4   5   ALT     String  The list (separated by ,) of alternative alleles at this position
        5   6   QUAL    Int     A quality score associated with the inference of the given alleles.
        6   7   FILTER  String  A flag indicating which of a given set of filters the variation has passed
        7   8   INFO    String  An extensible list of key-value pairs (fields) describing the variation, for example, NS=2;DP=10;AF=0.333,0.667;AA=T;DB
        8   9                   Optional fields...

        Returns: tuple of str or int
            The length of tuple might be > 8,
              depending of the presence of additional optional fields
        """
        line = self.__vcf.readline().rstrip()
        if line:
            fields = line.split('\t')
            for i in (1, 5):
                fields[i] = int(fields[i])
            return tuple(fields)
        else:  # line == ''
            return None

    def get_header(self):
        """Returns the header section (str) of SAM file"""
        return self.__header

    def close(self):
        self.__vcf.close()


def print_vcf(var=None):
    """
    Pretty print a variant (tuple) from vcf file

    Args:
        var: tuple or list of str or int
            Containing (at least) 8 fields of a variant from vcf file
    """
    if var is None:
        text = """\
#   Col	Field   Type    Description
0   1   CHROM   String  The name of the sequence (typically a chromosome)
1   2   POS     Int     The 1-based position of the variation on the given sequence
2   3   ID      String  The identifier of the variation
3   4   REF     String  The reference base (or bases in the case of an indel) at the given position on the given reference sequence
4   5   ALT     String  The list (separated by ,) of alternative alleles at this position
5   6   QUAL    Int     A quality score associated with the inference of the given alleles.
6   7   FILTER  String  A flag indicating which of a given set of filters the variation has passed
7   8   INFO    String  An extensible list of key-value pairs (fields) describing the variation, for example, NS=2;DP=10;AF=0.333,0.667;AA=T;DB
8   9                   Optional fields..."""
        printf(text)

    elif type(var) is tuple or type(var) is list:
        fields = ['CHROM ', 'POS   ', 'ID    ', 'REF   ', 'ALT   ',
                  'QUAL  ', 'FILTER', 'INFO  ']
        for i in range(8):
            printf(f"{i}\t{fields[i]}\t{var[i]}")
        if len(var) > 8:
            for i in range(8, len(var)):
                printf(f"{i}\t     \t{var[i]}")


def get_info(var):
    """
    Args:
        var: tuple or list of str or int
            Containing (at least) 8 fields of a variant from vcf file

    Return: dict
        Key, value pairs from the INFO (column 7) of a variant
        Values still str, not converted to int or float
    """
    func = lambda x: x.split('=')
    return {key: val for key, val in map(func, var[7].split(';'))}

