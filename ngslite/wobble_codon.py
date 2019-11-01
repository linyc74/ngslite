from itertools import combinations

IUPAC = {
    ('A', 'C'): 'M',
    ('G', 'T'): 'K',
    ('A', 'G'): 'R',
    ('C', 'T'): 'Y',
    ('C', 'G'): 'S',
    ('A', 'T'): 'W',
    ('C', 'G', 'T'): 'B',
    ('A', 'C', 'G'): 'V',
    ('A', 'G', 'T'): 'D',
    ('A', 'C', 'T'): 'H',
    ('A', 'C', 'G', 'T'): 'N',
}

CODON = {
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'TTA': 'L', 'TTG': 'L',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TTT': 'F', 'TTC': 'F',
    'ATG': 'M',
    'TGT': 'C', 'TGC': 'C',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'TAT': 'Y', 'TAC': 'Y',
    'TGG': 'W',
    'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N',
    'CAT': 'H', 'CAC': 'H',
    'GAA': 'E', 'GAG': 'E',
    'GAT': 'D', 'GAC': 'D',
    'AAA': 'K', 'AAG': 'K',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'TAA': '*', 'TAG': '*', 'TGA': '*',
}


def build_aa_dict(codon_dict):
    """
    Invert the key, value relationship in <codon_dict>

    <codon_dict> has codon (str) as key, and single aminoacid (str) as value

    The returned <aa_dict> has aminoacid (str) as key, and list of codons as value

    Args:
        codon_dict: dict
            e.g. {'ATT': 'I', 'ATC': 'I', 'ATA': 'I', ...}

    Returns: aa_dict
        e.g. {'I': ['ATT', 'ATC', 'ATA'], ...}
    """
    aa_dict = {}
    for codon, aa in codon_dict.items():
        aa_dict.setdefault(aa, [])
        aa_dict[aa].append(codon)
    return aa_dict


def generate_wobble_codons(codon_list):
    """
    Given a list of codons, generate all possible wobble codons
    with the third base represented by IUPAC ambiguity codes

    For example,

        ['AAT', 'AAC'] --> ['AAY']

        ['AAT', 'AAC', 'AAA'] --> ['AAY', 'AAM', 'AAW', 'AAH']

        ['AAT', 'AAC', 'GGT', 'GGC'] --> ['AAY', 'GGY']

    Args:
        codon_list: list of str

    Returns: list of str
    """
    wobble_codons = []

    # Build a dict with (the first two bases) as key, and (a list of third bases) as value
    # For example, codon_list = ['AAT', 'AAC',  -->  d = {'AA': ['T', 'C'],
    #                            'GGT', 'GGC']            'GG': ['T', 'C']}
    d = {}
    for codon in codon_list:
        d.setdefault(codon[0:2], [])
        d[codon[0:2]].append(codon[2])

    # Each value is a list, sort the list
    # d = {'AA': ['T', 'C'], --> d = {'AA': ['C', 'T'],
    #      'GG': ['T', 'C']}          'GG': ['C', 'T']}
    d = {key: sorted(val) for key, val in d.items()}

    for first_two, third_bases in d.items():
        if len(third_bases) == 1:
            # No ambiguity at the third base, skip
            continue

        # For example, if third_bases = ['A', 'C', 'G', 'T'],
        # then need to choose k = 2, 3, 4 from the third_bases
        for k in range(2, len(third_bases) + 1):
            for comb in combinations(third_bases, k):
                degenerated_base = IUPAC[comb]
                wobble_codons.append(first_two + degenerated_base)

    return wobble_codons


def build_wobble_codon_dict(codon_dict):
    """
    Build a dictionary with wobble codon (str) as key, and aminoacid (str) as value

    Wobble codon uses IUPAC ambiguity codes for the third base

    Args:
        codon_dict: dict
            e.g. {'ATT': 'I', 'ATC': 'I', 'ATA': 'I', ...}

    Returns: dict
        e.g. {'ATY': 'I', 'ATM': 'I', 'ATW': 'I', 'ATH': 'I', ...}
    """
    ret = {}

    aa_dict = build_aa_dict(codon_dict=codon_dict)

    for aa, codons in aa_dict.items():
        for wobble_codon in generate_wobble_codons(codon_list=codons):
            ret[wobble_codon] = aa

    return ret


def pretty_print_wobble_codon_dict():
    wobble_codon_dict = build_wobble_codon_dict(codon_dict=CODON)

    text = '{\n   '
    cur_prefix, cur_aa = None, None
    for wobble_codon, aa in wobble_codon_dict.items():

        if cur_prefix is None:
            cur_prefix = wobble_codon[:2]
        if cur_aa is None:
            cur_aa = aa

        if wobble_codon[:2] != cur_prefix or aa != cur_aa:
            char0 = '\n    '
            cur_prefix = wobble_codon[:2]
            cur_aa = aa
        else:
            char0 = ' '

        text += f'{char0}{wobble_codon}: {aa},'

    print(text + '\n}')


if __name__ == '__main__':
    pretty_print_wobble_codon_dict()
