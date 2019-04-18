from .gtftools import GtfWriter
from .fasta import FastaWriter


class Locus(object):
    """
    A data class to store the information of a feature (or locus, or interval)
    """
    def __init__(self, type_, start, end, strand, attributes):
        """
        Args:
            type_: str,
                Types defined in genbank file, for example 'gene', 'CDS', 'misc_feature'

            start: int,
                1-based, inclusive

            end: int,
                1-based, inclusive

            strand: str
                '+', '-'

            attributes: dict
                for example
                {
                    'gene': 'AAA  [M=132]',
                    'accession': 'PF00004.29',
                    'description': 'ATPase family associated with variouscellular activities (AAA)'
                }
        """
        self.type = type_
        self.start = start
        self.end = end
        self.strand = strand
        self.attributes = attributes


def __get_contig_id(contig_text):
    """
    Args:
        contig_text: str
            The complete text block of a contig (or genome) in a genbank file,
                starting from 'LOCUS' and ending at '//'

    Returns: str,
        The contig_id (or accession number) in the very first line
    """
    line1 = contig_text.splitlines()[0]
    # After 'LOCUS' and before several space
    contig_id = line1[len('LOCUS'):].lstrip().split(' '*3)[0]
    return contig_id


def __split_features(all_feature_text):
    """
    Split the complete feature text block into smaller blocks, in which
        each small block is a feature

    Args:
        all_feature_text: str
            for example

            'FEATURES             Location/Qualifiers
            '     source          1..168903
            '                     /organism="Enterobacteria phage T4"
            '                     /mol_type="genomic DNA"
            '                     /db_xref="taxon:10665"
            '     gene            complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

    Returns: list of str
        for example [str1, str2], where

        str1 =
            '     source          1..168903
            '                     /organism="Enterobacteria phage T4"
            '                     /mol_type="genomic DNA"
            '                     /db_xref="taxon:10665"

        str2 =
            '     gene            complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"
    """
    L = []
    for line in all_feature_text.splitlines()[1:]:
        if line.strip() == '':  # empty line with only blank spaces
            continue
        if line.startswith(' '*21):
            L[-1] = L[-1] + line + '\n'
        else:
            L.append(line + '\n')

    return L


def __get_feature_type(feature_text):
    """
    Args:
        feature_text: str,
            for example

            '     CDS             complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

    Returns: str
        The example would be 'CDS'
    """
    line1 = feature_text.splitlines()[0]
    return line1.strip().split()[0]


def __get_feature_location(feature_text):
    """
    Args:
        feature_text: str,
            for example

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
    end = int(region.split('..')[1])

    return start, end, strand


def __get_feature_attributes(feature_text):
    """
    Args:
        feature_text: str,
            for example

            '     CDS             complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

    Returns: dict
        The example would be
        {
            'gene': 'rIIA',
            'locus_tag': 'T4p001'
            'db_xref': 'GeneID:1258593'
        }
    """
    # Remove the first line
    feature_text = '\n'.join(feature_text.splitlines()[1:])

    # Parse as a list of str, each str is 'key="text"' or 'key=value'
    attr_list = list()
    for line in feature_text.splitlines():
        if line.startswith(' '*21+'/') and '=' in line:
            attr_list.append(line[22:])
        else:
            attr_list[-1] = attr_list[-1] + line.lstrip()

    # Unpack kay="text" or key=value into a dictionary
    attr_dict = dict()
    for a in attr_list:
        if '="' in a:
            key, val = a.split('="')
            val = val[:-1]
        elif '=' in a:
            key, val = a.split('=')
        else:
            continue
        attr_dict[key] = val

    return attr_dict


def __parse_feature(feature_text):
    """
    Parse a single feature text block and returns a Locus()
        that stores the information of the feature

    Args:
        feature_text: str,
            for example

            '     gene            complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

    Returns: Locus() object
    """
    type_ = __get_feature_type(feature_text)
    start, end, strand = __get_feature_location(feature_text)
    attributes = __get_feature_attributes(feature_text)
    return Locus(type_, start, end, strand, attributes)


def __extract_features(all_feature_text):
    """
    Take the complete feature text block of one contig,
        parse the features into a list of Locus() objects

    Args:
        all_feature_text: str, for example

            'FEATURES             Location/Qualifiers
            '     source          1..168903
            '                     /organism="Enterobacteria phage T4"
            '                     /mol_type="genomic DNA"
            '                     /db_xref="taxon:10665"
            '     gene            complement(12..2189)
            '                     /gene="rIIA"
            '                     /locus_tag="T4p001"
            '                     /db_xref="GeneID:1258593"

    Returns: list of Locus() object
    """
    L = []  # list of loci

    for feature_text in __split_features(all_feature_text):
        locus = __parse_feature(feature_text)
        if locus.start <= locus.end:
            L.append(locus)

    return L


def extract_genbank_features_to_dict(file):
    """
    Args:
        file: str, path-like
            The input genbank file

    Returns: dict()
        {
            contig_id_1: [locus1, locus2, ...],
            contig_id_2: [locus3, locus4, ...]
        }

        in which,
            contig_id is str
            locus is Locus() object
    """
    D = dict()

    with open(file) as fh:
        for contig_text in fh.read().strip().split('\n//\n'):

            contig_id = __get_contig_id(contig_text)

            # The <all_feature_text> is after 'FEATURES' and before 'ORIGIN'
            all_feature_text = contig_text.split('\nFEATURES')[1].split('\nORIGIN')[0]
            all_feature_text = 'FEATURES' + all_feature_text

            # Extract <all_feature_text> -> a list of Locus() objects
            feature_list = __extract_features(all_feature_text)

            D[contig_id] = feature_list

    return D


def genbank_to_gtf(file, output, skip_attributes):
    """
    Extract features in the genbank file into GTF file

    Args:
        file: str, path-like
            The input genbank

        output: str, path-like
            The output GTF

        skip_attributes: str, or list of str
            Attributes not to be included, for example ['translation', 'codon_start']
    """
    if type(skip_attributes) is str:
        skip_attributes = [skip_attributes]

    with GtfWriter(output) as writer:

        D = extract_genbank_features_to_dict(file)
        for contig_id, features in D.items():
            for locus in features:

                # Pack locus.attributes (dict) into a single line of text, e.g.
                #   gene "AAA  [M=132]";accession "PF00004.29"
                attributes = ''
                for key, val in locus.attributes.items():
                    if key in skip_attributes:
                        continue
                    attributes = attributes + f"{key} \"{val}\";"
                attributes = attributes[:-1]  # Remove trailing ;

                writer.write(
                    (contig_id,    # seqname
                     '.',          # source
                     locus.type,   # feature type
                     locus.start,  # start
                     locus.end,    # end
                     '.',          # score
                     locus.strand, # strand
                     0,            # frame
                     attributes)   # attribute
                )


def genbank_to_fasta(file, output):
    """
    Args:
        file: str, path-like
            The input genbank

        output: str, path-like
            The output fasta
    """
    with open(file) as fh:
        with FastaWriter(output) as writer:
            for contig_text in fh.read().strip().split('\n//\n'):

                contig_id = __get_contig_id(contig_text)

                # The formatted sequence text is after 'ORIGIN'
                seq_text = contig_text.strip().split('\nORIGIN\n')[1]
                seq = ''
                for line in seq_text.splitlines():
                    seq = seq + line[9:].replace(' ', '')

                writer.write(contig_id, seq)
