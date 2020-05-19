from datetime import date
from .fasta_gtf import read_fasta_gtf
from .dnatools import translate, rev_comp
from .genbank_parse import GenbankText


class GenbankTextWriter:
    def __init__(self, file: str, mode: str = 'w'):
        """
        Args
            file: path-like

            mode:
                'w' for write or 'a' for append
        """
        self.__gbk = open(file, mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def write(self, genbank_text: GenbankText):
        """
        Args:
            genbank_text: namedtuple GenbankItem
        """
        for i in range(3):
            self.__gbk.write(genbank_text[i].rstrip() + '\n')
        self.__gbk.write('//\n')

    def close(self):
        self.__gbk.close()


def get_today() -> str:
    """
    Returns: str,
        Today, e.g. '20-JUN-2006'
    """
    today = date.today()
    day = today.strftime('%d')
    month = today.strftime('%b').upper()
    year = today.strftime('%Y')
    return f"{day}-{month}-{year}"


def wrap_line_by_word(line, length, indent, sep=' ', keep_sep=False) -> str:
    """
    Wraps a line (not containing any '\n') by word
    The returned line does not end with '\n'

    Args:
        line: str

        length: int

        indent: int

        sep: str
            The separator of words, usually ' '

        keep_sep: bool
            When wrapping, whether or not to keep the separator attached to the word

    Returns: str
    """
    words = line.split(sep)
    if keep_sep:
        # separator always attached to the word (except the last word)
        for i, word in enumerate(words[:-1]):
            words[i] = word + sep

    indent = ' ' * indent
    text = indent + words[0]

    for word in words[1:]:
        last_line = text.split('\n')[-1]

        if keep_sep:
            # Separator is always attached to the word, no need to consider the separator's length
            if len(last_line) + len(word) <= length:
                text += word
            else:
                text += '\n' + indent + word

        else:
            # Need to consider the separator's length
            if len(last_line) + 1 + len(word) <= length:
                text += sep + word
            else:
                text += '\n' + indent + word

    return text


def wrap_line_by_char(line, length, indent) -> str:
    """
    Wraps a line (not containing any '\n') by character
    The returned line does not end with '\n'

    Args:
        line: str

        length: int

        indent: int

    Returns: str
    """
    text = ' ' * indent
    for char in line:
        last_line = text.split('\n')[-1]
        if len(last_line) + 1 <= length:
            text += char
        else:
            text += '\n' + ' '*indent + char
    return text


def generic_feature_to_genbank_text(feature):
    """
    Use the information in a GenericFeature object to write a genbank feature
        text section, for example:

    `     CDS             complement(join(<484868..485129,485184..485520,
    `                     485576..>486058))
    `                     /locus_tag="GN000052.00_006979"
    `                     /codon_start=3
    `                     /product="hypothetical protein"
    `                     /protein_id="ncbi:GN000052.00_006979-T1"
    `                     /translation="PAPQAPAPQAPAPQAPAPQAPAPHAPAPHAPAPHAPAQHAPAQQ
    `                     ASLQGAPAHQVQMPLGVPPPTPSSLRHIQSAIIREEPVYPEPVHRTQSPSPRTHGTKK
    `                     RANKEYTADSESGPSWGFDDVYGDGFKLADEDKTDLVELTEAQAKNVKTLGLAAGGKL
    `                     IQDIYKDPFPAAIWNHAGARILHVHMIDPESCERVTHIVPQPPPMAVEEYVKSGGQFF
    `                     VVEEKVDERLDGGDFDNVKSVSQMDQHLGITTEPEFDPKKPKMCTTCERRLCDCIIRP
    `                     CNHQFCNVCIKRLDEAGDTEQSVQQARHWKCPTCNSPVSHVAGFSAPMNLPGEERLLR
    `                     TKVPVHVLKIEDGRMRLSSMQTSRV"

    Args:
        feature: GenericFeature object

    Returns: str
    """
    f = feature

    lines = []

    # First line: regions
    regions = [f'{start}..{end}' for start, end, _ in f.regions]
    regions = ','.join(regions)
    if f.partial_start:
        regions = '<' + regions
    if f.partial_end:
        # Insert '>' after the last '..'
        pos = regions.rfind('..') + 2
        regions = regions[:pos] + '>' + regions[pos:]
    if len(f.regions) > 1:
        regions = f'join({regions})'
    if f.strand == '-':
        regions = f'complement({regions})'

    lines.append(regions)

    # Second to the last line
    for key, val in f.attributes:
        if type(val) is int or type(val) is float:
            lines.append(f'/{key}={val}')
        else:
            lines.append(f'/{key}="{val}"')

    text = ''
    for i, line in enumerate(lines):
        if i == 0:
            # Wrap the first line (regions) by ',' and keep ','
            text += wrap_line_by_word(line, length=79, indent=21, sep=',', keep_sep=True)
        elif line.startswith('/translation'):
            text += wrap_line_by_char(line, length=79, indent=21)
        else:
            text += wrap_line_by_word(line, length=79, indent=21, sep=' ', keep_sep=False)
        text += '\n'

    # Fill in the feature type, e.g. 'CDS'
    head = ' '*5 + f.type
    text = head + text[len(head):]

    return text


def make_header(
        molecule,
        length,
        shape,
        ACCESSION,
        DEFINITION,
        KEYWORDS,
        SOURCE,
        ORGANISM,
        division='ENV') -> str:
    """
    Create a header section of a genbank file

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
        Text of genbank header
    """

    # Replace all '\n' with ' ' to become single line
    ACCESSION  = ACCESSION.replace('\n', ' ')
    DEFINITION = DEFINITION.replace('\n', ' ')
    KEYWORDS   = KEYWORDS.replace('\n', ' ')
    SOURCE     = SOURCE.replace('\n', ' ')

    # Wrap and then remove the indent of the first line
    ACCESSION  = wrap_line_by_word(ACCESSION, length=79, indent=12).lstrip()
    DEFINITION = wrap_line_by_word(DEFINITION, length=79, indent=12).lstrip()
    KEYWORDS   = wrap_line_by_word(KEYWORDS, length=79, indent=12).lstrip()
    SOURCE     = wrap_line_by_word(SOURCE, length=79, indent=12).lstrip()

    # If there are multiple lines in ORGANISM
    if '\n' in ORGANISM:
        lines = ORGANISM.split('\n')
        # The 1st line of ORGANISM should be the species name
        # Starting from the 2nd line, it should be the taxanomy
        # Make everything from the 2nd line a single line without wrapping
        ORGANISM = lines[0] + '\n' + ' ' * 12 + ''.join(lines[1:])

    header = f"""\
LOCUS       {ACCESSION}        {length} bp    {molecule}    {shape}   {division} {get_today()}
DEFINITION  {DEFINITION}
ACCESSION   {ACCESSION}
KEYWORDS    {KEYWORDS}
SOURCE      {SOURCE}
  ORGANISM  {ORGANISM}\n"""

    return header


def init_feature(length, organism, mol_type='genomic DNA') -> str:
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
    text = f'''\
FEATURES             Location/Qualifiers
     source          1..{length}
                     /mol_type="{mol_type}"
                     /organism="{organism}"\n'''
    return text


def translate_feature(feature, sequence) -> str:
    """
    Args:
        feature: GenericFeature object

        sequence: str
            DNA sequence to which the feature belongs
    """
    # Build CDS by concatenating introns
    cds = ''
    for start, end, _ in feature.regions:
        cds += sequence[(start - 1):end]

    if feature.strand == '-':
        cds = rev_comp(cds)

    codon_start = feature.get_attribute('codon_start')
    if codon_start is not None:
        cds = cds[(codon_start-1):]

    translation = translate(cds)

    if translation.endswith('*'):
        translation = translation[:-1]  # remove stop codon

    if feature.partial_start and feature.strand == '+':
        return translation  # not full protein --> no need to start with 'M'
    elif feature.partial_end and feature.strand == '-':
        return translation  # not full protein --> no need to start with 'M'
    else:
        return 'M' + translation[1:]  # full protein --> start with 'M'


def format_ref_seq(seq) -> str:
    """
    Take a DNA sequence (str) and format it to the ORIGIN section of genbank file.

    Args:
        seq: str

    Returns: str
    """
    # Remove all ' ' and '\n'
    seq = seq.replace(' ', '').replace('\n', '')

    lines = []
    num_lines = int(len(seq) / 60) + 1  # 60 nucleotides each line
    for i in range(num_lines):
        new = ''  # new line
        pos = i * 60 + 1
        new += ' '*(9 - len(str(pos))) + str(pos)
        for j in range(6):
            a = pos + j * 10 - 1
            b = a + 10
            new += ' ' + seq[a:b]
        lines.append(new)

    return 'ORIGIN\n' + '\n'.join(lines) + '\n'


def write_genbank(
        data,
        file,
        DEFINITION='.',
        KEYWORDS='.',
        SOURCE='.',
        ORGANISM='.',
        genbank_division='ENV',
        use_locus_text=True):
    """
    Write <data> into a new genbank <file>

    Args:
        data: list of Chromosome objects, or dict

            [Chromosome_1, Chromosome_2, ...]

            or

            {
                Chromosome_1.seqname: Chromosome_1,
                Chromosome_2.seqname: Chromosome_2, ...
            }

        file: str, path-like
            The output genbank file

        DEFINITION: str
            The "DEFINITION" field in the genbank file

        KEYWORDS: str
            The "KEYWORDS" field in the genbank file

        SOURCE: str
            The "SOURCE" field in the genbank file

        ORGANISM: str
            The "ORGANISM" field in the genbank file

        genbank_division: str
            The three-letter tag for one of the 18 divisions in the GenBank database
            Default 'ENV' for environmental sampling sequences

        use_locus_text: bool
            Use the genbank_locus_text from the original input gbk file for the LOCUS section
    """
    if type(data) is dict:
        data = data.values()

    with GenbankTextWriter(file) as writer:
        for chromosome in data:
            sequence = chromosome.sequence

            # --- LOCUS text section --- #
            if use_locus_text:
                locus_text = chromosome.genbank_locus_text
            else:
                locus_text = make_header(
                    molecule='DNA',
                    length=len(sequence),
                    shape=['linear', 'circular'][chromosome.circular],
                    ACCESSION=chromosome.seqname,
                    DEFINITION=DEFINITION,
                    KEYWORDS=KEYWORDS,
                    SOURCE=SOURCE,
                    ORGANISM=ORGANISM,
                    division=genbank_division
                )

            # --- FEATURES text section --- #
            features_text = init_feature(
                length=len(sequence),
                organism=ORGANISM.split('\n')[0]
            )

            for f in chromosome.feature_array:
                # 'source' feature was made in the _init_feature()
                if f.type == 'source':
                    continue

                # For CDS, make 'translation', 'codon_start' if absent
                if f.type == 'CDS':
                    if f.get_attribute('translation') is None:
                        aa = translate_feature(feature=f, sequence=sequence)
                        f.add_attribute('translation', aa)
                    if f.get_attribute('codon_start') is None:
                        f.set_attribute('codon_start', f.frame)

                features_text += generic_feature_to_genbank_text(feature=f)

            # --- ORIGIN text sections --- #
            origin_text = format_ref_seq(sequence)

            # --- Write into file --- #
            writer.write(
                GenbankText(
                    locus_text=locus_text,
                    features_text=features_text,
                    origin_text=origin_text
                )
            )


def make_genbank(
        fasta,
        gtf,
        output,
        shape='linear',
        DEFINITION='.',
        KEYWORDS='.',
        SOURCE='.',
        ORGANISM='.',
        genbank_division='ENV'):
    """
    Merge a fasta file and a GTF file into a genbank file, in which
        the fasta file provides the genome sequence and
        the GTF file provides annotation (e.g. CDS).

    The fasta file can contain more than one sequence. Headers in the fasta file
        is used as the LOCUS name (i.e. ACCESSION) in the genbank file.

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
    data = read_fasta_gtf(
        fasta=fasta,
        gtf=gtf,
        circular=(shape == 'circular')
    )

    write_genbank(
        data=data,
        file=output,
        DEFINITION=DEFINITION,
        KEYWORDS=KEYWORDS,
        SOURCE=SOURCE,
        ORGANISM=ORGANISM,
        genbank_division=genbank_division
    )
