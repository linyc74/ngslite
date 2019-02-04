from .gtftools import GtfParser
from .fasta import FastaParser
from .dnatools import translate, rev_comp


from datetime import date


def __today():
    """
    Returns: str,
        Today, e.g. '20-JUN-2006'
    """
    today = date.today()
    day = today.strftime('%d')
    month = today.strftime('%b').upper()
    year = today.strftime('%Y')
    return f"{day}-{month}-{year}"


def __wrap_paragraph(paragraph, length, indent, by_char=False):
    """
    Wraps a paragraph (contains '\n'). The returned paragraph ends with '\n'.

    Args:
        paragraph: str
        length: int
        indent: int
        by_char: bool

    Returns: str
    """
    T = ''
    # Wrap line by line
    for line in paragraph.rstrip().splitlines():
        if by_char:
            T = T + __wrap_line_by_char(line, length, indent) + '\n'
        else:  # wrap by word
            T = T + __wrap_line_by_word(line, length, indent) + '\n'
    return T


def __wrap_line_by_word(line, length, indent):
    """
    Wraps a line (usually does not contain '\n') by word.
    The returned line does not end with '\n'.

    Args:
        line: str
        length: int
        indent: int

    Returns: str
    """
    L = ''
    for word in line.split(' '):
        last_line = L.split('\n')[-1]
        if len(last_line) + 1 + len(word) <= length:
            L = L + ' ' + word
        else:
            L = L + '\n' + ' '*indent + word
    return L[1:]


def __wrap_line_by_char(line, length, indent):
    """
    Wraps a line (usually does not contain '\n') by character.
    The returned line does not end with '\n'.

    Args:
        line: str
        length: int
        indent: int

    Returns: str
    """
    L = ''
    for char in line:
        last_line = L.split('\n')[-1]
        if len(last_line) + 1 <= length:
            L = L + char
        else:
            L = L + '\n' + ' '*indent + char
    return L


def __read_gtf_attribute(feature):
    """
    From a feature of GTF file,
        read the attributes (in the last field) as a dict.

    Args:
        feature: tuple or list

    Returns: dict
    """
    D = {}
    for attr in feature[-1].split(';'):
        attr = attr.strip()
        # Use the first empty space to separate key and "value"
        pos = attr.find(' ')
        key, val = attr[:pos], attr[pos+1:]
        # Remove the open and close quotes
        val = val[1:-1]
        # Add the key value pair to the dictionary
        D[key] = val
    return D


# HEADER
def __make_header(molecule, length, shape, ACCESSION, DEFINITION, KEYWORDS,
                     SOURCE, ORGANISM, division='ENV'):
    """
    Create a header section (str) of a genbank file.

    Args:
        molecule: str, e.g. 'DNA', 'mRNA'
        length: int
        shape: str, 'linear' or 'circular'
        ACCESSION: str
        DEFINITION: str
        KEYWORDS: str
        SOURCE: str
        ORGANISM: str
        division: str,
            The three-letter tag for one of the 18 divisions in the GenBank database.
            Default 'ENV' for environmental sampling sequences

    Returns: str
        A text (paragraph) of genbank header
    """

    # Replace all '\n' with ' ' to become single line
    ACCESSION  = ACCESSION.replace('\n', ' ')
    DEFINITION = DEFINITION.replace('\n', ' ')
    KEYWORDS   = KEYWORDS.replace('\n', ' ')
    SOURCE     = SOURCE.replace('\n', ' ')

    # The 1st line of ORGANISM should be the species name
    ORGANISM_1 = ORGANISM.split('\n')[0]
    # Starting from the 2nd line, it should be the taxanomy
    # Make everything from the 2nd line a single line
    ORGANISM_2 = ' '.join(ORGANISM.split('\n')[1:])

    text = f"""\
LOCUS       {ACCESSION}        {length} bp    {molecule}    {shape}   {division} {__today()}
DEFINITION  {DEFINITION}
ACCESSION   {ACCESSION}
KEYWORDS    {KEYWORDS}
SOURCE      {SOURCE}
  ORGANISM  {ORGANISM_1}
            {ORGANISM_2}\n"""

    lines = text.splitlines()
    # block = the text block from 2nd to the last line
    block = '\n'.join(lines[1:])
    # Do not wrap the 1st line; wrap everything starting from the 2nd line
    block = __wrap_paragraph(block, length=79, indent=12)
    # Put the first line and the block together -> genbank header section
    header = lines[0] + '\n' + block

    return header


# FEATURE
def __init_feature(length, organism, mol_type='genomic DNA'):
    """
    Initialize a feature section of genbank file.
    The feature section always start with the mandatory 'source' feature

    Args:
        length: int
            Length of the genome

        organism: str

        mol_type: str
            Usually 'genomic DNA'
    """
    text = f"""\
FEATURES             Location/Qualifiers
     source          1..{length}
                     /mol_type="{mol_type}"
                     /organism="{organism}"\n"""
    return text


def __add_feature(key, start, end, strand, **kwargs):
    """
    Create a feature of genbank file.
    The returned str can be added to an existing feature section of a genbank file.

    Example:

    '     CDS             complement(205..390)
    '                     /locus_tag="XPXV15_gp01"
    '                     /codon_start=1
    '                     /transl_table=11
    '                     /product="hypothetical protein"
    '                     /translation="MSRKELRTSESEWRRTKHTAYMPKGTHVTRTWRSQAIKVATVIM
    '                     IGTIITAFWYEFLKGAA"
     1234567890123456789012

    Args:
        key: str
            e.g. 'CDS', 'gene'

        start: int

        end: int

        strand: str, '+' or '-'

        kwargs:
            qualifier=value pairs
            For example, if qualifier='locus_tag' and value='ABC001', the following will be added:
                /locus_tag="ABC001"

    Returns: str
    """
    # First line, e.g.: '     CDS             complement(205..390)'
    text = ' '*5 + key + ' '*(21-5-len(key))
    if strand == '-':
        text = text + f"complement({start}..{end})\n"
    else:  # strand == '+'
        text = text + f"{start}..{end}\n"

    # From the 2nd line
    for qualifier, value in kwargs.items():
        if type(value) is int or type(value) is float:
            # If value is a number, add it without quotes
            line = ' '*21 + f"/{qualifier}={value}"
        else:  # type(value) is str
            # If value is str, add it with quotes
            line = ' ' * 21 + f"/{qualifier}=\"{value}\""

        if qualifier == 'translation':
            # If it's a protein sequence, there's no space, wrap the line by character
            line = __wrap_line_by_char(line, length=79, indent=21)
        else:
            line = __wrap_line_by_word(line, length=79, indent=21)

        # Add the /qualifier="value" line to the text block
        text = text + line + '\n'

    return text


# ORIGIN
def __format_ref_seq(seq):
    """
    Take a DNA sequence (str) and format it to the ORIGIN section of genbank file.

    Args:
        seq: str

    Returns: str
    """
    # Remove all ' ' and '\n'
    seq = seq.replace(' ', '').replace('\n', '')

    text = 'ORIGIN\n'
    # Each line contains 60 nucleotides
    num_lines = int(len(seq) / 60) + 1
    for i in range(num_lines):
        pos = i * 60 + 1
        text = text + ' '*(9 - len(str(pos))) + str(pos)
        for j in range(6):
            a = pos + j * 10 - 1
            b = a + 10
            text = text + ' ' + seq[a:b]
        text = text + '\n'
    return text


# The master public function
def make_genbank(fasta, gtf, output, shape='linear', DEFINITION='.',
                 KEYWORDS='.', SOURCE='.', ORGANISM='.', genbank_division='ENV'):
    """
    Merge a fasta file and a GTF file into a genbank file, in which
        the fasta file provides the genome sequence and
        the GTF file provides annotation (e.g. CDS).

    The fasta file can contain more than one sequence. Headers in the fasta file
        is used as the ACCESSION in the genbank file.

    GFF and GTF formats are both supported.

    Args:
        fasta: str, path-like
            The input fasta file

        gtf: str, path-like
            The input GTF (or GFF) file

        output: str, path-like
            The output genbank file

        shape: str
            The shape of DNA or RNA, either 'linear' or 'circular'

        DEFINITION: str
            The "DEFINITION" field in the genbank file

        KEYWORDS: str
            The "KEYWORDS" field in the genbank file

        SOURCE: str
            The "SOURCE" field in the genbank file

        ORGANISM: str
            The "ORGANISM" field in the genbank file

        genbank_division: str,
            The three-letter tag for one of the 18 divisions in the GenBank database.
            Default 'ENV' for environmental sampling sequences
    """

    # Group the features in the GTF file in a dictionary
    # {'seqname': [list of features]}
    #                      where each feature is a tuple of 9 fields
    feature_dict = {}
    for ftr in GtfParser(gtf):
        seqname = ftr[0]
        feature_dict.setdefault(seqname, [])
        feature_dict[seqname].append(ftr)

    gbk = open(output, 'w')
    
    for genome_id, genome_seq in FastaParser(fasta):
        # --- HEADER --- #
        text = __make_header(molecule = 'DNA',
                             length = len(genome_seq),
                             shape = shape,
                             ACCESSION = genome_id,
                             DEFINITION = DEFINITION,
                             KEYWORDS = KEYWORDS,
                             SOURCE = SOURCE,
                             ORGANISM = ORGANISM,
                             division=genbank_division)
        gbk.write(text)

        # --- FEATURE --- #
        text = __init_feature(length=len(genome_seq),
                              organism=ORGANISM.split('\n')[0])
        gbk.write(text)

        # Get the features belonging to the current genome
        features = feature_dict.get(genome_id, [])
        for ftr in features:
            ftr_type, start, end, strand, frame = ftr[2], ftr[3], ftr[4], ftr[6], ftr[7]

            # Read the attributes of the feature as a dictionary
            attr_dict = __read_gtf_attribute(ftr)

            # If the feature is a CDS, we need to add:
            #   protein sequence (tranlation)
            #   frame (codon start)
            if ftr_type == 'CDS':
                if strand == '+':
                    cds = genome_seq[start - 1:end]
                else:  # strand == '-'
                    cds = rev_comp(genome_seq[start - 1:end])
                attr_dict['translation'] = translate(cds)
                attr_dict['codon_start'] = frame + 1  # codon start is 1-based

            text = __add_feature(key=ftr_type,
                                 start=start,
                                 end=end,
                                 strand=strand,
                                 **attr_dict)
            gbk.write(text)

        # --- ORIGIN --- #
        text = __format_ref_seq(genome_seq)
        gbk.write(text)

        # Finally, write '//' to separate genomes
        gbk.write('//\n')

    gbk.close()

