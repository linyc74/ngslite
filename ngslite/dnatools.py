def rev_comp(seq):
    """
    Returns reverse complementary sequence of the input DNA string.
    """
    comp = {'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A',
            'N': 'N',
            'M': 'K',
            'K': 'M',
            'R': 'Y',
            'Y': 'R',

            'a': 't',
            'c': 'g',
            'g': 'c',
            't': 'a',
            'n': 'n',
            'm': 'k',
            'k': 'm',
            'r': 'y',
            'y': 'r'}
    return ''.join([comp[base] for base in seq[::-1]])


def translate(dna):
    """
    Translate a DNA sequence.
    If the DNA length is not a multiple of 3, then leave the last 1 or 2 bases untranslated.

    Args:
        dna: str,
            the DNA sequence

    Returns: str,
        the translated amino acid sequence
    """
    codon = {
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'CTT': 'L',
    'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'TTA': 'L',
    'TTG': 'L', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
    'GTG': 'V', 'TTT': 'F', 'TTC': 'F', 'ATG': 'M',
    'TGT': 'C', 'TGC': 'C', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'GGT': 'G', 'GGC': 'G',
    'GGA': 'G', 'GGG': 'G', 'CCT': 'P', 'CCC': 'P',
    'CCA': 'P', 'CCG': 'P', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TGG': 'W', 'CAA': 'Q',
    'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'CAT': 'H',
    'CAC': 'H', 'GAA': 'E', 'GAG': 'E', 'GAT': 'D',
    'GAC': 'D', 'AAA': 'K', 'AAG': 'K', 'CGT': 'R',
    'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R',
    'AGG': 'R', 'TAA': '*', 'TAG': '*', 'TGA': '*'}

    dna = dna.upper()
    peptide = []
    for i in range(int(len(dna)/3)):
        aa = codon.get(dna[(i*3):(i*3+3)], 'X')
        peptide.append(aa)
    return ''.join(peptide)


def base_content(seq, base):
    """
    Compute the fraction of <base> in the <seq>.

    Args
        seq: str

        base: str
            e.g., 'A' gives the A content; 'GC' gives the GC content; 'ACGT' should give 1.0

    Returns: float
    """
    seq = seq.upper()
    base = base.upper()

    ret = 0
    for b in set(base):
        ret += seq.count(b) / len(seq)
    return ret

