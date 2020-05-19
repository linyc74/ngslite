from typing import List, Union, Optional, Tuple, Dict
from .lowlevel import call, printf
from .file_conversion import sam_to_bam


def sort_bam(file: str, keep: bool = False):
    """
    Args:
        file: path-like

        keep:
            Keep the input file or not
    """
    file_out = f"{file[:-4]}_sorted.{file[-3:]}"
    call(f"samtools sort {file} > {file_out}")
    if not keep:
        call(f"rm {file}")
        call(f"mv {file_out} {file}")


def index_bam(file: str):
    """
    Args:
        file: path-like
    """
    cmd = f"samtools index {file}"
    call(cmd)


def sam_to_indexed_bam(file: str, keep: bool = True):
    """
    Args:
        file: path-like
            The input sam file

        keep:
            Keep the input file or not
    """
    sam_to_bam(file, keep=keep)
    bam = file[:-4] + '.bam'
    sort_bam(bam, keep=False)
    index_bam(bam)


def subset_bam_regions(
        file: str,
        regions: List[str],
        output: Optional[str] = None,
        keep: bool = True):
    """
    Args:
        file: path-like
            The input bam file

        regions:
            Each str is a region of the reference genome, e.g.
                chr1            chromosome 1
                chr3:1000-2000  chromosome 3 from 1000th (inclusive) to 2000th (inclusive) base

        output: path-like
            The output bam file
            If None, add subscript '_subset' to the input bam file

        keep:
            If False, delete the input <file> and rename the output as the input <file>
            Overrides the <output> file name
    """
    # Convert regions list into a string
    # ['chr1', 'chr2:1001-2000'] -> ' chr1 chr2:1001-2000'
    regions = ' ' + ' '.join(regions)  # space-separated regions

    if output is None:
        output = file[:-len('.bam')] + '_subset.bam'

    # -b: output is a bam
    # -h: include header section
    call(f"samtools view -b -h {file}{regions} > {output}")

    if not keep:
        call(f"rm {file}")
        call(f"mv {output} {file}")


def remove_unmapped(file: str, keep: bool = True):
    """
    Remove unmapped reads from an input SAM/BAM file
        by using the option '-F 4' of 'samtools view' command

    Args:
        file: path-like
            The input sam/bam file

        keep:
            Keep the input file or not
    """
    b = ['', '-b '][file.endswith('.bam')]  # Output file is bam or not
    file_out = f"{file[:-4]}_remove_unmapped.{file[-3:]}"

    # -b: output is a bam
    # -h: include header section
    # -F 4: NOT including the flag 'read unmapped'
    call(f"samtools view -h -F 4 {b}{file} > {file_out}")

    if not keep:
        call(f"rm {file}")
        call(f"mv {file_out} {file}")


class SamParser:
    """
    A SAM file parser
    """
    def __init__(self, file: str):
        """
        The header section is parsed upon instantiation.

        Args:
            file: path-like
        """
        self.__sam = open(file, 'r')
        header = ''
        while True:
            # Get the current position
            pos = self.__sam.tell()
            # Readline and move on to the next position
            line = self.__sam.readline()
            if line.startswith('@'):
                header = header + line
            else:
                # If reaching the alignment section, that is,
                #   the header section has been parsed completely,
                #   then go back one line and break out the loop
                self.__sam.seek(pos)
                # Remember where the first line of data is
                self.pos_0 = pos
                break
        # Store the header string
        self.header = header

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def __iter__(self):
        self.__sam.seek(self.pos_0)
        return self

    def __next__(self):
        r = self.next()
        if r:
            return r
        else:  # r is None
            raise StopIteration

    def next(self) -> Optional[Tuple[Union[str, int]]]:
        """
        Each line of the SAM file has at least 11 fields

        #   Col	Field   Type    Description
        0   1   QNAME   String  Query template NAME
        1   2   FLAG    Int     bitwise FLAG
        2   3   RNAME   String  References sequence NAME
        3   4   POS     Int     1-based leftmost mapping POSition
        4   5   MAPQ    Int     MAPping Quality
        5   6   CIGAR   String  CIGAR String
        6   7   RNEXT   String  Ref. name of the mate/NEXT read
        7   8   PNEXT   Int     Position of the mate/NEXT read
        8   9   TLEN    Int     observed Template LENgth
        9   10  SEQ     String  segment SEQuence
        10  11  QUAL    String  ASCII of Phred-scaled base QUALity+33
        11  12  optional fields...

        Returns:
            The length of tuple might be > 11,
              depending of the presence of additional optional fields
        """
        line = self.__sam.readline().rstrip()
        if line:
            fields: List[Union[str, int]] = line.split('\t')
            for i in (1, 3, 4, 7, 8):
                fields[i] = int(fields[i])
            return tuple(fields)
        else:  # line == ''
            return None

    def close(self):
        self.__sam.close()


class SamWriter:
    def __init__(self, file: str, header: str, mode: str = 'w'):
        """
        The header section should be written in upon instantiation.

        Args:
            file: path-like

            header:
                The header section

            mode:
                The file mode: 'w' or 'a'
        """
        self.__sam = open(file, mode)
        if not header == '':
            self.__sam.write(header.rstrip() + '\n')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def write(self, read: Tuple[Union[str, int]]):
        """
        Args:
            read:
                Containing 11 (at least) fields of a line of SAM file
        """
        self.__sam.write('\t'.join(map(str, read)) + '\n')

    def close(self):
        self.__sam.close()


FLAGS_PROPERTIES = (
    'read paired',  # 0
    'read mapped in proper pair',  # 1
    'read unmapped',  # 2
    'mate unmapped',  # 3
    'read reverse strand',  # 4
    'mate reverse strand',  # 5
    'first in pair',  # 6
    'second in pair',  # 7
    'not primary alignment',  # 8
    'read fails platform/vendor quality checks',  # 9
    'read is PCR or optical duplicate',  # 10
    'supplementary alignment'  # 11
)


def decode_flag(flag: Union[int, str], return_tuple: bool = False):
    """
    The flag in SAM files is a decimal integer of 12 bits for 12 properties.
    The range of the integer is 0 - 4095 (2^12 -1)
    This function decodes the flag and returns an easily readable format.

    000000000000
    ||||||||||||
    |||||||||||read paired
    ||||||||||read mapped in proper pair
    |||||||||read unmapped
    ||||||||mate unmapped
    |||||||read reverse strand
    ||||||mate reverse strand
    |||||first in pair
    ||||second in pair
    |||not primary alignment
    ||read fails platform/vendor quality checks
    |read is PCR or optical duplicate
    supplementary alignment

    Args:
        flag:
            decimal flag, ranging from 0 - 4095

        return_tuple:
            If True, returns a tuple of True/False, to speed up the function

    Returns:
        A dictionary of the properties of the flag, with
            each key as the property name
            each value as True/False
    """
    # Format decimal to binary strings
    bin_str = '{0:b}'.format(int(flag))
    # Pad zeros up to 12 bits
    bin_str = '0'*(12-len(bin_str)) + bin_str
    # Starting from the right-most bit, make a list of True/False
    bin_tuple = tuple(bool(int(b)) for b in bin_str[::-1])
    if return_tuple:
        return bin_tuple
    # Make {PROPERTY: True/False} pairs and return the dictionary
    return {key: val for key, val in zip(FLAGS_PROPERTIES, bin_tuple)}


def encode_flag(flag_dict: Dict[str, bool], return_str: bool = False):
    """
    Encode a dictionary of flags into an int, which represents the 12-bit binary string

    Args:
        flag_dict:
            For example:

            {'first in pair': True,
             'read reverse strand': False},

            Flags not specified are defaulted False

        return_str: bool
            If True, returns a 12-bit str

    Returns: int or str
    """

    bits = ['0'] * 12
    for key, val in flag_dict.items():
        # This flag is True, set the corresponding bit to '1'
        if val:
            which = FLAGS_PROPERTIES.index(key)
            bits[which] = '1'
    # First flag is on the right (last element), so reverse the list <bits>
    bits = bits[::-1]
    bits = ''.join(bits)
    if return_str:
        return bits
    return int(bits, 2)  # Binary string -> int


def filter_sam_by_flag(
        file_in: str,
        file_out: str,
        flag_sets: List[Dict[str, bool]]):
    """
    Args:
        file_in: path-like
            The input SAM file

        file_out: path-like
            The output SAM file

        flag_sets:
            Each dictionary is a set of flags that have to be completely satisfied.
            There could be more than one set of flags, therefore many dictionaries.
            If a read satisfies any one of the dictionary, then it is included, so other
                there's no need to look into other dictionaries in the list.

            For example:

            [
                {'first in pair': True,
                 'read reverse strand': False},

                {'first in pair': False,
                 'read reverse strand': True},
            ]

            The two dictionaries represent two sets of flags:
                "read 1 mapped to the forward strand"
                "read 2 mapped to the reverse strand"

            These two sets of flags indicate the DNA insert is from the forward strand.
    """
    parser = SamParser(file_in)
    writer = SamWriter(file_out, header=parser.header)

    for A in parser:
        decoded_flag = decode_flag(A[1])  # FLAG is the second field of each line
        for each_set in flag_sets:
            # each_set is a dict

            passed = True
            for key in each_set.keys():
                if each_set[key] != decoded_flag[key]:
                    # If any flag among the set is not satisfied,
                    #   break from the loop of the dictionary <each_set>
                    passed = False
                    break

            if passed:
                # If all of the flags among the set is satisfied,
                #   write the read into the output sam file
                writer.write(A)
                # Break from the loop of all sets,
                # A read can only be included once, even if it satisfies other set of flags.
                break

    parser.close()
    writer.close()


def print_flag(flag: Optional[Union[int, Dict[str, bool]]] = None):
    """
    Decode (int -> dict) or encode (dict -> int) a flag and then print it

    Args: None or int or dict
        If None:
            Print an example dictionary of flags

        If int:
            Print decoded (int -> dict) flag

        If dict:
            print encode (dict -> int) flag
    """
    if type(flag) == dict:
        printf(encode_flag(flag))
        return

    elif flag is None:
        D = {key: True for key in FLAGS_PROPERTIES}

    elif type(flag) == int:
        D = decode_flag(flag)

    t = ''
    for key, val in D.items():
        t = t + f"'{key}': {val},\n "
    t = '{' + t[:-3] + '}'
    printf(t)


def print_sam(read: Tuple[Union[str, int]] = None):
    """
    Pretty print a read (tuple) from sam file

    Args:
        read:
            Containing (at least) 11 fields of a read from sam file
    """
    if read is None:
        text = """\
#   Col	Field   Type    Description
0   1   QNAME   String  Query template NAME
1   2   FLAG    Int     bitwise FLAG
2   3   RNAME   String  References sequence NAME
3   4   POS     Int     1-based leftmost mapping POSition
4   5   MAPQ    Int     MAPping Quality
5   6   CIGAR   String  CIGAR String
6   7   RNEXT   String  Ref. name of the mate/NEXT read
7   8   PNEXT   Int     Position of the mate/NEXT read
8   9   TLEN    Int     observed Template LENgth
9   10  SEQ     String  segment SEQuence
10  11  QUAL    String  ASCII of Phred-scaled base QUALity+33"""
        printf(text)

    elif type(read) is tuple or type(read) is list:
        fields = ['QNAME', 'FLAG ', 'RNAME', 'POS  ', 'MAPQ ',
                  'CIGAR', 'RNEXT', 'PNEXT', 'TLEN ', 'SEQ  ', 'QUAL ']
        for i in range(11):
            printf(f"{i}\t{fields[i]}\t{read[i]}")
        if len(read) > 11:
            for i in range(11, len(read)):
                printf(f"{i}\t     \t{read[i]}")
