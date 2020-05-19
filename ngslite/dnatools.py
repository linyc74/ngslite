def rev_comp(seq: str) -> str:
    """
    Returns reverse complementary sequence of the input DNA string
    """
    comp = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'N': 'N',
        'M': 'K',  # M = A C
        'K': 'M',  # K = G T
        'R': 'Y',  # R = A G
        'Y': 'R',  # Y = C T
        'S': 'S',  # S = C G
        'W': 'W',  # W = A T
        'B': 'V',  # B = C G T
        'V': 'B',  # V = A C G
        'D': 'H',  # D = A G T
        'H': 'D',  # H = A C T

        'a': 't',
        'c': 'g',
        'g': 'c',
        't': 'a',
        'n': 'n',
        'm': 'k',
        'k': 'm',
        'r': 'y',
        'y': 'r',
        's': 's',
        'w': 'w',
        'b': 'v',
        'v': 'b',
        'd': 'h',
        'h': 'd'
    }
    return ''.join([comp[base] for base in seq[::-1]])


def translate(dna: str) -> str:
    """
    Translate a DNA sequence
    If the DNA length is not a multiple of 3,
        then leave the last 1 or 2 bases untranslated
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
        'AGG': 'R', 'TAA': '*', 'TAG': '*', 'TGA': '*',

        'GTK': 'V'
    }
    dna = dna.upper()
    peptide = []
    for i in range(int(len(dna)/3)):
        aa = codon.get(dna[(i*3):(i*3+3)], 'X')
        peptide.append(aa)
    return ''.join(peptide)


def base_content(seq: str, base: str) -> float:
    """
    Compute the fraction of <base> in the <seq>

    Args
        seq

        base:
            e.g., 'A' gives the A content; 'GC' gives the GC content; 'ACGT' should give 1.0
    """
    seq = seq.upper()
    base = base.upper()

    ret = 0
    for b in set(base):
        ret += seq.count(b) / len(seq)
    return ret
