import subprocess


from .file_conversion import *


def __call(cmd):
    print('CMD: ' + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as inst:
        print(inst)


def sort_bam(file, keep=False):
    """
    Args:
        file: str, path-like

        keep: bool
            Keep the input file or not
    """
    file_out = '{}_sorted.{}'.format(file[:-4], file[-3:])
    cmd = 'samtools sort {} > {}'.format(file, file_out)
    __call(cmd)
    if not keep:
        __call('rm {}'.format(file))
        __call('mv {} {}'.format(file_out, file))


def index_bam(file):
    """
    Args:
        file: str, path-like
    """
    cmd = 'samtools index {}'.format(file)
    __call(cmd)


def sam_to_indexed_bam(file, keep=True):
    """
    Args:
        file: str, path-like
            The input sam file

        keep: bool
            Keep the input file or not
    """
    sam_to_bam(file, keep=keep)
    bam = file[:-4] + '.bam'
    sort_bam(bam, keep=False)
    index_bam(bam)


def subset_bam_regions(file, regions, output=None, keep=True):
    """
    Args:
        file: str, path-like
            The input bam file

        regions: list of str
            Each str is a region of the reference genome, e.g.
                chr1            chromosome 1
                chr3:1000-2000  chromosome 3 from 1000th (inclusive) to 2000th (inclusive) base

        output: str, path-like
            The output bam file
            If None, add subscript '_subset' to the input bam file

        keep: bool
            If False, delete the input <file> and rename the output as the input <file>
    """
    # Convert regions list into a string
    # ['chr1', 'chr2'] -> '"chr1" "chr2"'
    add_quote = lambda x: '"{}"'.format(x)
    regions = ' '.join(map(add_quote, regions))

    if output is None:
        output = file[:-len('.bam')] + '_subset.bam'

    # -b: output is a bam
    # -h: include header section
    cmd = 'samtools view -b -h {} {} > {}'.format(file, regions, output)
    __call(cmd)

    if not keep:
        __call('rm {}'.format(file))
        __call('mv {} {}'.format(output, file))


class SamParser:
    """
    A SAM file parser
    """
    def __init__(self, file):
        """
        The header section is parsed upon instantiation.

        Args:
            file: str, path-like object
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
                break
        # Store the header string
        self.__header = header

    def next(self):
        """
        Each line of the SAM file has at least 11 fields

            Col	Field   Type    Brief Description
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

        Returns: tuple of str or int
            The length of tuple might be > 11,
              depending of the presence of additional optional fields
        """
        line = self.__sam.readline().rstrip()
        if line == '':
            return (None, ) * 11
        else:
            fields = line.split('\t')
            for i in (1, 3, 4, 7, 8):
                fields[i] = int(fields[i])
            return tuple(fields)

    def get_header(self):
        """Returns the header section (str) of SAM file"""
        return self.__header

    def close(self):
        self.__sam.close()


class SamWriter:
    def __init__(self, file, header, mode='w'):
        """
        The header section should be written in upon instantiation.

        Args:
            file: str, path-like object

            header: str
                The header section

            mode: str, 'w' or 'a'
        """
        self.__sam = open(file, mode)
        if not header == '':
            self.__sam.write(header.rstrip() + '\n')

    def write(self, alignment):
        """
        Args:
            alignment: tuple of str or int
                Containing 11 (at least) fields of a line of SAM file
        """
        self.__sam.write('\t'.join(map(str, alignment)) + '\n')

    def close(self):
        self.__sam.close()


# global variable for decode_flag()
FLAGS_PROPERTIES = ('read paired',  # 0
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
                    'supplementary alignment')  # 11
def decode_flag(flag, return_tuple=False):
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
        flag: int or str
            decimal flag, ranging from 0 - 4095

        return_tuple: bool
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


def filter_sam_by_flag(file_in, file_out, flag_sets):
    """
    Args:
        file_in: str, path-like object
            The input SAM file

        file_out: str, path-like object
            The output SAM file

        flag_sets: list of dictionary
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
    writer = SamWriter(file_out, header=parser.get_header())

    while True:
        A = parser.next()
        if A[0] is None:  # The end of input file, so break
            break

        flag = decode_flag(A[1])  # FLAG is the second field of each line
        for each_set in flag_sets:
            # each_set is a dict

            passed = True
            for key in each_set.keys():
                if not flag[key] == each_set[key]:
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

