from collections import namedtuple
from typing import List, Tuple, Union, Dict, Optional
from .gtftools import GtfWriter
from .gfftools import GffWriter
from .fasta import FastaWriter
from .dataclass import GenericFeature, FeatureArray, Chromosome
from .feature_conversion import generic_to_gtf_feature, generic_to_gff_feature


GenbankText = namedtuple('GenbankText', 'locus_text features_text origin_text')


def get_seqname(locus_text: str) -> str:
    """
    Args:
        locus_text:
            The LOCUS (i.e. header section) of genbank file

    Returns:
        The seqname of chromosome
    """
    line1 = locus_text[:locus_text.find('\n')]
    return line1[len('LOCUS'):].lstrip().split(' '*3)[0]


def get_circular(locus_text: str) -> bool:
    """
    Args:
        locus_text:
            The LOCUS (i.e. header section) of genbank file
    """
    line1 = locus_text[:locus_text.find('\n')]
    if ' circular ' in line1: return True
    return False


def split_features_text(features_text: str) -> List[str]:
    """
    Args:
        features_text:
            The FEATURES (i.e. annotation) section of genbank file, for example:

            'FEATURES             Location/Qualifiers
            '     source          1..168903
            '                     /organism="Enterobacteria phage T4"
            '                     /mol_type="genomic DNA"
            '                     /db_xref="taxon:10665"
            '     CDS             complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

    Returns:
        Break it into a list of single feature text, for example

        feature_text_1:
            '     source          1..168903
            '                     /organism="Enterobacteria phage T4"
            '                     /mol_type="genomic DNA"
            '                     /db_xref="taxon:10665"

        feature_text_2:
            '     CDS             complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"
    """
    list_ = []
    for line in features_text.splitlines()[1:]:  # first line is useless
        if line.strip() == '':  # skip empty line
            continue
        if line.startswith(' '*21):
            list_[-1] += line + '\n'
        else:
            list_.append(line + '\n')
    return list_


def get_feature_type(feature_text: str) -> str:
    """
    Args:
        feature_text: endswith '\n'
            For example:
            '     CDS             complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

    Returns:
        The example would be 'CDS'
    """
    line1 = feature_text.splitlines()[0]
    return line1.strip().split()[0]


def get_location_string(feature_text: str) -> str:
    """
    Args:
        feature_text: endswith '\n'
            For example:
            '     CDS             complement(join(<360626..360849,360919..360948,
            '                     361067..361220,361292..361470,361523..>361555))
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

    Returns:
        The example would be
        'complement(join(<360626..360849,360919..360948,361067..361220,361292..361470,361523..>361555))'
    """
    # The location string is before ' '*21 + '/', sometimes could be multiple lines
    pos = feature_text.find(' ' * 21 + '/')
    ret = feature_text[:pos]

    # Multiple -> single line
    ret = ret.replace('\n' + ' ' * 21, '')

    # Remove the feature type such as 'CDS'
    ret = ret.strip().split()[1]

    return ret


def get_feature_location(feature_text: str) -> Tuple[int, int, str, List, bool, bool]:
    """
    Args:
        feature_text: endswith '\n'
            For example:
            '     CDS             complement(join(<360626..360849,360919..360948,
            '                     361067..361220,361292..361470,361523..>361555))
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

    Returns:
        start: int
            The example would be 12

        end: int
            The example would be 2189

        strand: str, '+' or '-'
            The example would be '-'

        regions: list of tuple (int, int, str)
            Indicating start, end, strand of each region (intron)
            The example would be [(12, 2189, '-')]

        partial_start: bool

        partial_end: bool
    """
    locstr = get_location_string(feature_text)

    if locstr.startswith('complement('):
        # Remove 'complement(' and ')'
        locstr = locstr[len('complement('):-1]
        all_complement = True
    else:
        all_complement = False

    if locstr.startswith('join('):
        # Remove 'join(' and ')'
        locstr = locstr[len('join('):-1]

    # loclist = list of strings
    # e.g. ["complement(2853..2990)", "complement(2458..2802))"]
    loclist = locstr.split(',')

    partial_start, partial_end = False, False
    regions = []  # e.g. [(100, 200, '-'), (<300, 400, '+'), (500, >600, '+')]
    for i, s in enumerate(loclist):
        # Tell the strand
        if s.startswith('complement('):
            # Remove 'complement(' and ')'
            s = s[len('complement('):-1]
            c = '-'           # c = strand
        elif all_complement:
            c = '-'
        else:
            c = '+'

        a, b = s.split('..') if ('..' in s) else (s, s)  # a is start, b is end

        # First start has '<' --> partial start
        if i == 0 and a.startswith('<'):
            a = a[1:]
            partial_start = True

        # Last end has '>' --> partial end
        if i == len(loclist) - 1 and b.startswith('>'):
            b = b[1:]
            partial_end = True

        a, b = int(a), int(b)

        if a > b:
            a, b = b, a  # a must be < b

        regions.append((a, b, c))

    regions = sorted(regions, key=lambda x: x[0])
    start, end, strand = regions[0][0], regions[-1][1], regions[0][2]

    return start, end, strand, regions, partial_start, partial_end


def get_feature_attributes(feature_text: str) -> List[Tuple[str, str]]:
    """
    Args:
        feature_text: endswith '\n'
            For example:
            '     CDS             complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

    Returns: list of tuples of (key, value) pairs
        The example would be
        [
            ('gene', 'rIIA'),
            ('locus_tag', 'T4p001'),
            ('db_xref', 'GeneID:1258593')
        ]
    """
    # If no attribute at all, return empty list
    if (' '*21 + '/') not in feature_text:
        return []

    # Remove the location information before ' '*21 + '/'
    pos = feature_text.find(' '*21 + '/')
    feature_text = feature_text[pos:]

    # Parse the feature_text into a list of str, where each str is 'key="text"' or 'key=value'
    attr_list = list()
    for line in feature_text.splitlines():
        if line.startswith(' '*21+'/') and '=' in line:
            attr_list.append(line[22:])
        else:
            # Add a blank space ' ' for changing line
            attr_list[-1] += ' ' + line.lstrip()

    # Unpack key="value" or key=value into a list of tuples (key, value)
    for i, attr in enumerate(attr_list):

        # With quote, value is text
        if '="' in attr:
            pos = attr.find('="')
            key = attr[:pos]
            val = attr[pos+2:-1]

        # Without quote
        elif '=' in attr:
            key, val = attr.split('=')
            if val.isdigit():
                val = int(val)
        else:
            continue

        # Protein sequence should not contain blank space
        if key == 'translation':
            val = val.replace(' ', '')

        attr_list[i] = (key, val)

    return attr_list


def construct_generic_feature(feature_text: str, seqname: str) -> GenericFeature:
    """
    Construct a GenericFeature object from a single feature text
    Introns (joined regions) are supported

    Args:
        feature_text: endswith '\n'
            For example:
            '     CDS             complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

        seqname: str
    """

    start, end, strand, regions, partial_start, partial_end = get_feature_location(feature_text)

    return GenericFeature(
        seqname=seqname,
        type_=get_feature_type(feature_text),
        start=start,
        end=end,
        strand=strand,
        regions=regions,
        frame=1,
        attributes=get_feature_attributes(feature_text),
        partial_start=partial_start,
        partial_end=partial_end
    )


def contruct_feature_array(
        features_text: str, seqname: str, genome_size: int, circular: bool) -> FeatureArray:
    """
    Construct a FeatureArray from the complete FEATURES section of genbank file

    Args:
        features_text: str
            The FEATURES (i.e. annotation) section of genbank file

            For example:
            'FEATURES             Location/Qualifiers
            '     source          1..168903
            '                     /organism="Enterobacteria phage T4"
            '                     /mol_type="genomic DNA"
            '                     /db_xref="taxon:10665"
            '     CDS             complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

        seqname: str

        genome_size: int

        circular: bool
    """
    features = []

    for feature_text in split_features_text(features_text):

        features.append(
            construct_generic_feature(feature_text, seqname)
        )

    return FeatureArray(
        seqname=seqname, genome_size=genome_size,
        features=features, circular=circular)


def get_sequence(origin_text: str) -> str:
    """
    Args:
        origin_text: str
            The ORIGIN (i.e. sequence section) of genbank file

    Returns: str
        The genome sequence
    """
    lines = []

    for line in origin_text.splitlines():
        lines.append(line[9:].replace(' ', ''))

    return ''.join(lines)


def is_valid_first_line_of_feature(line: str) -> bool:

    if len(line) <= 21:
        return False

    if ' ' in line[:21].strip():
        return False

    if not line[21].isdigit() and \
            not line[21:].startswith('complement') and \
            not line[21:].startswith('join') and \
            not line[21] == '<':
        return False

    return True


def get_valid_features_text(genbank_text: GenbankText) -> str:

    lines = []

    for line in genbank_text.features_text.splitlines():

        if line == 'FEATURES             Location/Qualifiers':
            lines.append(line)

        elif line.startswith(' ' * 21):
            lines.append(line)

        elif line.startswith(' ' * 5):
            if is_valid_first_line_of_feature(line=line):
                lines.append(line)

    return '\n'.join(lines)


def construct_chromosome(genbank_text: GenbankText) -> Chromosome:

    locus_text = genbank_text.locus_text
    features_text = get_valid_features_text(genbank_text=genbank_text)
    origin_text = genbank_text.origin_text

    seqname = get_seqname(locus_text=locus_text)
    sequence = get_sequence(origin_text=origin_text)
    circular = get_circular(locus_text=locus_text)

    features = contruct_feature_array(
        features_text=features_text,
        seqname=seqname,
        genome_size=len(sequence),
        circular=circular)

    return Chromosome(
        seqname=seqname,
        sequence=sequence,
        features=features,
        circular=circular,
        genbank_locus_text=locus_text)


class GenbankTextParser:

    def __init__(self, file: str):
        """
        Args:
            file: path-like
        """
        self.__gbk = open(file)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def __iter__(self):
        self.reset()
        return self

    def __next__(self):
        r = self.next()
        if r:
            return r
        else:  # r is None
            raise StopIteration

    def reset(self):
        self.__gbk.seek(0)

    def next(self) -> Optional[GenbankText]:
        pos = self.__gbk.tell()
        if not self.__gbk.readline().startswith('LOCUS'):
            return None

        self.__gbk.seek(pos)
        lines = []
        while True:
            line = self.__gbk.readline()
            if line.startswith('//'):
                break
            line = line[:-1]  # remove '\n'
            lines.append(line)
        text = '\n'.join(lines)

        locus_text = text[0:text.find('FEATURES')]
        features_text = text[text.find('FEATURES'):text.find('ORIGIN')]
        origin_text = text[text.find('ORIGIN'):]

        return GenbankText(locus_text, features_text, origin_text)

    def close(self):
        self.__gbk.close()


class GenbankParser:

    def __init__(self, file):

        self.__text_parser = GenbankTextParser(file=file)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def __iter__(self):
        self.reset()
        return self

    def __next__(self):
        r = self.next()
        if r:
            return r
        else:  # r is None
            raise StopIteration

    def reset(self):
        self.__text_parser.reset()

    def next(self) -> Optional[Chromosome]:
        genbank_text = self.__text_parser.next()
        if genbank_text is None:
            return None
        else:
            return construct_chromosome(genbank_text)

    def close(self):
        self.__text_parser.close()


def read_genbank(
        file: str,
        as_dict: bool = False) -> \
        Union[List[Chromosome], Dict[str, Chromosome]]:
    """
    Read a genbank file into Chromosome ojbects, i.e. annotated genome

    Args:
        file: path-like
            The input genbank file

        as_dict:
            Return as dictionary instead of list

    Returns: list of Chromosome objects, or dict

        [Chromosome_1, Chromosome_2, ...]

        or

        {
            Chromosome_1.seqname: Chromosome_1,
            Chromosome_2.seqname: Chromosome_2, ...
        }
    """
    with GenbankParser(file) as parser:
        chromosomes = [chrom for chrom in parser]

    if as_dict:
        chromosomes = {chrom.seqname: chrom for chrom in chromosomes}

    return chromosomes


def genbank_to_fasta(file: str, output: str):
    """
    Args:
        file: path-like
            Input genbank file

        output: path-like
            Output fna file
    """
    with FastaWriter(output) as writer:
        with GenbankParser(file) as parser:
            for chromosome in parser:
                writer.write(chromosome.seqname, chromosome.sequence)


DEFAULT_SKIP_TYPES = ['source']
DEFAULT_SKIP_ATTRIBUTES = ['translation', 'codon_start', 'transl_table']


def genbank_to_gtf(
        file: str, output: str,
        skip_types: Optional[Union[str, List[str]]] = None,
        skip_attributes: Optional[Union[str, List[str]]] = None):
    """
    Extract features in the genbank file and write them into a GTF file

    Args:
        file: path-like
            Input genbank file

        output: path-like
            Output GTF file

        skip_types:
            Feature of types not to be included

        skip_attributes:
            Attributes not to be included
    """
    if skip_types is None:
        skip_types = DEFAULT_SKIP_TYPES
    elif type(skip_types) is str:
        skip_types = [skip_types]

    if skip_attributes is None:
        skip_attributes = DEFAULT_SKIP_ATTRIBUTES
    elif type(skip_attributes) is str:
        skip_attributes = [skip_attributes]

    with GtfWriter(output) as writer:
        with GenbankParser(file) as parser:
            for chromosome in parser:
                for feature in chromosome.features:
                    if feature.type in skip_types:
                        continue
                    for key in skip_attributes:
                        feature.remove_attribute(key)
                    gtf_feature = generic_to_gtf_feature(feature)
                    writer.write(gtf_feature)


def genbank_to_gff(
        file: str, output: str,
        skip_types: Optional[Union[str, List[str]]] = None,
        skip_attributes: Optional[Union[str, List[str]]] = None):
    """
    Extract features in the genbank file and write them into a GFF3 file

    Args:
        file: path-like
            Input genbank file

        output: path-like
            Output GFF file

        skip_types:
            Feature of types not to be included

        skip_attributes:
            Attributes not to be included
    """
    if skip_types is None:
        skip_types = DEFAULT_SKIP_TYPES
    elif type(skip_types) is str:
        skip_types = [skip_types]

    if skip_attributes is None:
        skip_attributes = DEFAULT_SKIP_ATTRIBUTES
    elif type(skip_attributes) is str:
        skip_attributes = [skip_attributes]

    with GffWriter(output) as writer:
        with GenbankParser(file) as parser:
            for chromosome in parser:
                for feature in chromosome.features:
                    if feature.type in skip_types:
                        continue
                    for key in skip_attributes:
                        feature.remove_attribute(key)
                    gff_feature = generic_to_gff_feature(feature)
                    writer.write(gff_feature)
