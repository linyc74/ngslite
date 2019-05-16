from collections import namedtuple
from .gtftools import GtfWriter
from .fasta import FastaWriter
from .data_class import GenericFeature, FeatureArray, Chromosome


GenbankItem = namedtuple('GenbankItem', 'LOCUS FEATURES ORIGIN')


class GenbankParser(object):
    def __init__(self, file):
        """
        Args:
            file: str, path-like object
        """
        self.__gbk = open(file, 'r')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def __iter__(self):
        self.__gbk.seek(0)
        return self

    def __next__(self):
        r = self.next()
        if r:
            return r
        else:  # r is None
            raise StopIteration

    def next(self):
        pos = self.__gbk.tell()
        if not self.__gbk.readline().startswith('LOCUS'):
            return None

        self.__gbk.seek(pos)
        lines = []
        while True:
            line = self.__gbk.readline()[:-1]  # Remove trailing '\n'
            if line.startswith('//'): break
            lines.append(line)
        text = '\n'.join(lines)

        LOCUS = text[0:text.find('FEATURES')]
        FEATURES = text[text.find('FEATURES'):text.find('ORIGIN')]
        ORIGIN = text[text.find('ORIGIN'):]

        return GenbankItem(LOCUS, FEATURES, ORIGIN)

    def close(self):
        self.__gbk.close()


def _get_seqname(LOCUS):
    """
    Args:
        LOCUS: str
            The LOCUS (i.e. header section) of genbank file

    Returns: str
        The unique seqname of chromosome
    """
    line1 = LOCUS[:LOCUS.find('\n')]
    return line1[len('LOCUS'):].lstrip().split(' '*3)[0]


def _is_circular(LOCUS):
    """
    Args:
        LOCUS: str
            The LOCUS (i.e. header section) of genbank file

    Returns: bool
    """
    line1 = LOCUS[:LOCUS.find('\n')]
    if ' circular ' in line1: return True
    return False


def _split_feature_text(FEATURES):
    """
    Args:
        FEATURES: str
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

    Returns: list of str
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
    L = []
    for line in FEATURES.splitlines()[1:]:  # first line is useless
        if line.strip() == '':  # skip empty line
            continue
        if line.startswith(' '*21):
            L[-1] += line + '\n'
        else:
            L.append(line + '\n')
    return L


def _get_feature_type(feature_text):
    """
    Args:
        feature_text: str

            '     CDS             complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

    Returns: str
        The example would be 'CDS'
    """
    line1 = feature_text.splitlines()[0]
    return line1.strip().split()[0]


def _get_feature_location(feature_text):
    """
    Args:
        feature_text: str,

            '     CDS             complement(12..2189)
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
    """
    line1 = feature_text.splitlines()[0]
    region = line1.strip().split()[1]

    if region.startswith('complement('):
        strand = '-'
        region = region[len('complement('):-1]
    else:
        strand = '+'

    start = int(region.split('..')[0])
    if '..' in region:
        end = int(region.split('..')[1])
    else:
        end = start

    return start, end, strand


def _get_feature_attributes(feature_text):
    """
    Args:
        feature_text: str,

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
    # Remove the first line
    feature_text = '\n'.join(feature_text.splitlines()[1:])

    # Parse as a list of str, each str is 'key="text"' or 'key=value'
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


def _construct_feature(feature_text, seqname):
    """
    Construct a GenericFeature object from a single feature text

    Introns (joined regions) are not supported

    Args:
        feature_text: str

            '     CDS             complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

        seqname: str
    """
    if 'join(' in feature_text.split('\n')[0]:
        return None

    start, end, strand = _get_feature_location(feature_text)
    return GenericFeature(
        seqname=seqname,
        type_=_get_feature_type(feature_text),
        start=start,
        end=end,
        strand=strand,
        attributes=_get_feature_attributes(feature_text),
        frame=1
    )


def _contruct_feature_array(FEATURES, seqname, genome_size, circular):
    """
    Construct a FeatureArray from the complete FEATURE section of genbank file

    Args:
        FEATURES: str
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

        seqname: str

        genome_size: int

        circular: bool
    """
    features = []
    for feature_text in _split_feature_text(FEATURES):
        # f is GenericFeature
        f = _construct_feature(feature_text, seqname)
        if f is not None:
            features.append(f)
    return FeatureArray(seqname, genome_size, features, circular)


def _get_sequence(ORIGIN):
    """
    Args:
        ORIGIN: str
            The ORIGIN (i.e. sequence section) of genbank file

    Returns: str
        The genome sequence
    """
    lines = []
    for line in ORIGIN.splitlines():
        lines.append(line[9:].replace(' ', ''))
    return ''.join(lines)


def _pack_attributes(feature, skip_attributes):
    """
    Pack feature.attributes into a sinlge line of str,
        which is in the format of the GTF attribute field, for example:

        gene "AAA  [M=132]";accession "PF00004.29";...

    Args:
        feature: GenericFeature() object

        skip_attributes: str, or list of str
            Attributes not to be included, for example ['translation', 'codon_start']

    Returns: str
    """
    if len(feature.attributes) == 0:
        return '.'
    s = ''
    for key, val in feature.attributes:
        if not key in skip_attributes:
            s += f"{key} \"{val}\";"
    return s[:-1]  # Remove trailing ";"


def construct_chromosome(genbank_item):
    """
    Construct a Chromosome object from a GenbankItem object

    Args:
        genbank_item: GenbankItem (namedtuple)
            Containing 3 text sections: LOCUS, FEATURES, ORIGIN

    Returns: Chromosome object
    """
    item = genbank_item

    seqname = _get_seqname(LOCUS=item.LOCUS)
    sequence = _get_sequence(ORIGIN=item.ORIGIN)
    circular = _is_circular(LOCUS=item.LOCUS)

    feature_array = _contruct_feature_array(
        FEATURES=item.FEATURES,
        seqname=seqname,
        genome_size=len(sequence),
        circular=circular
    )

    return Chromosome(
        seqname=seqname,
        sequence=sequence,
        feature_array=feature_array,
        circular=circular,
        genbank_LOCUS=item.LOCUS
    )


def read_genbank(file, as_dict=False):
    """
    Read a genbank file into Chromosome ojbects, i.e. annotated genome

    Args:
        file: str, path-like
            The input genbank file

        as_dict: bool

    Returns: list of Chromosome objects, or dict

        [Chromosome_1, Chromosome_2, ...]

        or

        {
            Chromosome_1.seqname: Chromosome_1,
            Chromosome_2.seqname: Chromosome_2, ...
        }
    """
    with GenbankParser(file) as parser:
        chromosomes = [construct_chromosome(item) for item in parser]
    if as_dict:
        return {chrom.seqname: chrom for chrom in chromosomes}
    return chromosomes


def genbank_to_fasta(file, output):
    """
    Args:
        file: str, path-like
            The input genbank file

        output: str, path-like
            The output fasta file
    """
    with FastaWriter(output) as writer:
        with GenbankParser(file) as parser:
            for item in parser:
                chrom = construct_chromosome(item)
                writer.write(chrom.seqname, chrom.sequence)


# Default skip types and attributes
default_skip_types = ['source']
default_skip_attributes = ['translation', 'codon_start', 'transl_table']

def genbank_to_gtf(file, output, skip_types=None, skip_attributes=None):
    """
    Extract features in the genbank file and write them into a GTF file

    Args:
        file: str, path-like
            The input genbank

        output: str, path-like
            The output GTF

        skip_types: str, or list of str
            Feature of types not to be included
            By default 'source' is skipped because it's just from the start to the end of genome

        skip_attributes: str, or list of str
            Attributes not to be included
            By default skip "translation", 'codon_start' and 'transl_table'
    """
    if skip_types is None:
        skip_types = default_skip_types
    elif type(skip_types) is str:
        skip_types = [skip_types]

    if skip_attributes is None:
        skip_attributes = default_skip_attributes
    elif type(skip_attributes) is str:
        skip_attributes = [skip_attributes]

    with GtfWriter(output) as writer:
        with GenbankParser(file) as parser:
            for item in parser:
                chrom = construct_chromosome(item)
                for feature in chrom.feature_array:
                    if feature.type in skip_types:
                        continue
                    for key in skip_attributes:
                        feature.remove_attribute(key)
                    gtf_feature = feature.to_gtf_feature()
                    writer.write(gtf_feature)
