from .fasta_gtf import read_fasta_gtf
from .dnatools import translate, rev_comp
from .genbank_parse import GenbankItem
from datetime import date


class GenbankWriter(object):
    def __init__(self, file, mode='w'):
        """
        Args
            file: str, path-like object
            mode: str, 'w' for write or 'a' for append
        """
        self.__gbk = open(file, mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def write(self, genbank_item):
        """
        Args:
            genbank_item: namedtuple GenbankItem
        """
        for i in range(3):
            self.__gbk.write(genbank_item[i].rstrip() + '\n')
        self.__gbk.write('//\n')

    def close(self):
        self.__gbk.close()


def _today():
    """
    Returns: str,
        Today, e.g. '20-JUN-2006'
    """
    today = date.today()
    day = today.strftime('%d')
    month = today.strftime('%b').upper()
    year = today.strftime('%Y')
    return f"{day}-{month}-{year}"


def _wrap_paragraph(paragraph, length, indent, by_char=False):
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
            T = T + _wrap_line_by_char(line, length, indent) + '\n'
        else:  # wrap by word
            T = T + _wrap_line_by_word(line, length, indent) + '\n'
    return T


def _wrap_line_by_word(line, length, indent):
    """
    Wraps a line (not containing any '\n') by word.
    The returned line does not end with '\n'.

    Args:
        line: str

        length: int

        indent: int

    Returns: str
    """
    words = line.split(' ')
    L = ' ' * indent + words[0]
    for word in words[1:]:
        last_line = L[L.rfind('\n')+1:]
        if len(last_line) + 1 + len(word) <= length:
            L += ' ' + word
        else:
            L += '\n' + ' '*indent + word
    return L


def _wrap_line_by_char(line, length, indent):
    """
    Wraps a line (not containing any '\n') by character.
    The returned line does not end with '\n'.

    Args:
        line: str

        length: int

        indent: int

    Returns: str
    """
    L = ' ' * indent
    for char in line:
        last_line = L[L.rfind('\n') + 1:]
        if len(last_line) + 1 <= length:
            L = L + char
        else:
            L = L + '\n' + ' '*indent + char
    return L


def _generic_feature_to_genbank_text(generic_feature):
    """
    Use the information in <generic_feature> object to write
        a genbank feature text section, for example:

        '     CDS             complement(12..2189)
        '                     /gene="rIIA"
        '                     /locus_tag="T4p001"
        '                     /db_xref="GeneID:1258593"

    Args:
        generic_feature: GenericFeature object

    Returns: str
    """
    f = generic_feature

    text = ''

    # --- First line --- #
    first_line = f"{' '*5}{f.type}{' '*(16 - len(f.type))}"

    ps = ['', '<'][f.partial_start]
    pe = ['', '>'][f.partial_end]

    # Contiguous segment, no intron
    if len(f.regions) == 1:
        if f.strand == '+':
            first_line += f"{ps}{f.start}..{f.end}{pe}"
        else:
            first_line += f"complement({ps}{f.start}..{f.end}{pe})"

    # Introns
    else:  # len(f.regions) >= 2
        s = ''
        for r in f.regions:
            if r[2] == '+':
                s += f"{r[0]}..{r[1]},"
            else:
                s += f"complement({r[0]}..{r[1]}),"
        first_line += f"join({s[:-1]})"  # Remove trailing ','

    text += first_line + '\n'

    # --- Second line to the end --- #
    for key, val in f.attributes:
        if type(val) is int or type(val) is float:
            newline = f"/{key}={val}"
        else:
            newline = f"/{key}=\"{val}\""

        if key == 'translation':
            text += _wrap_line_by_char(newline, length=79, indent=21) + '\n'
        else:
            text += _wrap_line_by_word(newline, length=79, indent=21) + '\n'

    return text


def _make_header(molecule, length, shape, ACCESSION, DEFINITION, KEYWORDS,
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

    # Wrap and then remove the indent of the first line
    ACCESSION  = _wrap_line_by_word(ACCESSION , length=79, indent=12).lstrip()
    DEFINITION = _wrap_line_by_word(DEFINITION, length=79, indent=12).lstrip()
    KEYWORDS   = _wrap_line_by_word(KEYWORDS  , length=79, indent=12).lstrip()
    SOURCE     = _wrap_line_by_word(SOURCE    , length=79, indent=12).lstrip()

    # If there are multiple lines in ORGANISM
    if '\n' in ORGANISM:
        lines = ORGANISM.split('\n')
        # The 1st line of ORGANISM should be the species name
        # Starting from the 2nd line, it should be the taxanomy
        # Make everything from the 2nd line a single line without wrapping
        ORGANISM = lines[0] + '\n' + ' ' * 12 + ''.join(lines[1:])

    header = f"""\
LOCUS       {ACCESSION}        {length} bp    {molecule}    {shape}   {division} {_today()}
DEFINITION  {DEFINITION}
ACCESSION   {ACCESSION}
KEYWORDS    {KEYWORDS}
SOURCE      {SOURCE}
  ORGANISM  {ORGANISM}\n"""

    return header


def _init_feature(length, organism, mol_type='genomic DNA'):
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


def _format_ref_seq(seq):
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


def write_genbank(data, file, DEFINITION='.', KEYWORDS='.', SOURCE='.',
                  ORGANISM='.', genbank_division='ENV'):
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

        genbank_division: str,
            The three-letter tag for one of the 18 divisions in the GenBank database.
            Default 'ENV' for environmental sampling sequences
    """
    if type(data) is dict:
        data = data.values()

    with GenbankWriter(file) as writer:
        for chromosome in data:
            sequence = chromosome.sequence

            # --- LOCUS text section --- #
            if chromosome.genbank_LOCUS == '':
                LOCUS = _make_header(
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
            else:
                LOCUS = chromosome.genbank_LOCUS

            # --- FEATURES text section --- #
            FEATURES = _init_feature(
                length=len(sequence),
                organism=ORGANISM.split('\n')[0]
            )

            for f in chromosome.feature_array:
                if f.type == 'source': continue

                # For CDS, re-make 'translation', 'codon_start' and 'transl_table'
                if f.type == 'CDS':
                    if f.strand == '+':
                        cds = sequence[(f.start - 1):f.end]
                    else:  # strand == '-'
                        cds = rev_comp(sequence[(f.start - 1):f.end])
                    # Always start with M regardless of GUG or other alternate start codons
                    # Remove stop codon *
                    aa = 'M' + translate(cds)[1:-1]
                    f.set_attribute('translation', aa)
                    f.set_attribute('codon_start', f.frame)
                    f.set_attribute('transl_table', 11)
                FEATURES += _generic_feature_to_genbank_text(f)

            # --- ORIGIN text sections --- #
            ORIGIN = _format_ref_seq(sequence)

            # --- Write into file --- #
            writer.write(
                GenbankItem(
                    LOCUS=LOCUS,
                    FEATURES=FEATURES,
                    ORIGIN=ORIGIN
                )
            )


def make_genbank(fasta, gtf, output, shape='linear', DEFINITION='.',
                 KEYWORDS='.', SOURCE='.', ORGANISM='.', genbank_division='ENV'):
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

