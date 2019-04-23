from .gtftools import read_gtf
from .fasta import FastaParser
from .dnatools import translate, rev_comp
from datetime import date
from .data_class import gtf_to_generic_feature


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
            L = L + ' ' + '\n' + ' '*indent + word
    return L[1:]


def _wrap_line_by_char(line, length, indent):
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


def _generic_feature_to_genbank_text(generic_feature):
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
    else:  # len(self.regions) >= 2
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
            newline = ' ' * 21 + f"/{key}={val}\n"
        else:
            newline = ' ' * 21 + f"/{key}=\"{val}\"\n"

        if key == 'translation':
            text += _wrap_line_by_char(newline, length=79, indent=21)
        else:
            text += _wrap_line_by_word(newline, length=79, indent=21)

    return text


# HEADER
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

    # The 1st line of ORGANISM should be the species name
    ORGANISM_1 = ORGANISM.split('\n')[0]
    # Starting from the 2nd line, it should be the taxanomy
    # Make everything from the 2nd line a single line
    ORGANISM_2 = ' '.join(ORGANISM.split('\n')[1:])

    text = f"""\
LOCUS       {ACCESSION}        {length} bp    {molecule}    {shape}   {division} {_today()}
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
    block = _wrap_paragraph(block, length=79, indent=12)
    # Put the first line and the block together -> genbank header section
    header = lines[0] + '\n' + block

    return header


# FEATURE
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


# ORIGIN
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


# The master public function
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

    # Group the features in the GTF file in a dictionary
    # {'seqname': [GtfFeature, ...]}
    feature_dict = read_gtf(gtf, as_dict=True)

    gbk = open(output, 'w')
    
    for genome_id, genome_seq in FastaParser(fasta):
        # --- HEADER --- #
        text = _make_header(
            molecule = 'DNA',
            length = len(genome_seq),
            shape = shape,
            ACCESSION = genome_id,
            DEFINITION = DEFINITION,
            KEYWORDS = KEYWORDS,
            SOURCE = SOURCE,
            ORGANISM = ORGANISM,
            division=genbank_division
        )
        gbk.write(text)

        # --- FEATURE --- #
        text = _init_feature(
            length=len(genome_seq),
            organism=ORGANISM.split('\n')[0],
        )
        gbk.write(text)

        # Get the features belonging to the current genome
        features = feature_dict.get(genome_id, [])
        for f in features:
            f = gtf_to_generic_feature(f)

            # If the feature is a CDS, need to add:
            #   protein sequence (translation)
            #   frame (codon start)
            if f.type == 'CDS':
                if f.strand == '+':
                    cds = genome_seq[(f.start-1):f.end]
                else:  # strand == '-'
                    cds = rev_comp(genome_seq[(f.start-1):f.end])
                f.add_attribute('translation', translate(cds))
                f.add_attribute('codon_start', f.frame)
            text = _generic_feature_to_genbank_text(f)
            gbk.write(text)

        # --- ORIGIN --- #
        text = _format_ref_seq(genome_seq)
        gbk.write(text)

        # Finally, write '//' to separate genomes
        gbk.write('//\n')

    gbk.close()
