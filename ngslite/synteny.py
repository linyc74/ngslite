from .fasta_gtf import read_fasta_gtf
from .genbank_parse import read_genbank
from .genbank_write import write_genbank
from .data_class import GenericFeature
import numpy as np
import re
from random import randint
from itertools import combinations
from scipy.cluster import hierarchy


score = {
    'gap': -1,
    'common_tag': +100,
    'CDS_intergenic': -1,
    'CDS_CDS_same_strand': +2,
    'CDS_CDS_opposite_strand': -5,
    'intergenic_intergeic': 0,
    'else': 0
}

score_scalar = 0.01


def _tag_feature(feature, qualifiers, skip='orf[0-9]+'):
    """
    Use all of the words in the <qualifiers> to tag the feature

    Args:
        feature: GenericFeature object

        qualifiers: str, or list of str

        skip: str, regex pattern
            Attribute words matching this pattern will be skipped
            Usually we want to skip very generic patterns like 'orf#####'
    """
    for key, val in feature.attributes:
        if key in qualifiers:
            for word in re.split('[|]|[;]|[,]', str(val)):
                if not re.search(skip, word):
                    feature.tags.append(word.strip())
    # Remove duplicates and make the same combination of words always identical
    feature.tags = list(set(feature.tags))


def _has_common_tag(tags1, tags2, skip='orf[0-9]+'):
    """
    Args:
        tags1: list of str

        tags2: list of str

        skip: str, regex pattern
            Tags matching this pattern will be skipped
            Usually we want to skip very generic patterns like 'orf#####'

    Returns: bool
    """
    for t in tags1:
        if t in tags2:
            if not re.search(skip, t):
                return True
    return False


def _make_color_dict(chromosome_dict):
    """
    From all Feature objects of all Chromosome objects,
        use tags to build a color dictionary, which has the following data structure:
        {
            tuple(feature.tags): color, ...
        {

        in which color is str and has the format: #ffffff (white)

    Args:
        chromosome_dict: dict
        {
            chromosome1.seqname: chromosome1,
            chromosome2.seqname: chromosome2, ...
        }

    Returns: dict
    """
    color_dict = {}
    for chromosome in chromosome_dict.values():
        for feature in chromosome.feature_array:
            # no tag -> skip
            if len(feature.tags) == 0: continue
            # feature.tags is list, convert to tuple for hashing
            key = tuple(feature.tags)
            color_dict.setdefault(key, f'#{randint(0, 256**3):06x}')
    return color_dict


def _set_feature_color(chromosome_dict, color_dict):
    """
    For each feature, use its tag to get a color from color_dict
        and set the color attribute with key='Color', val=color

    For example, red color = #ff0000

    Args:
        chromosome_dict: dict

        color_dict: dict
    """
    for chromosome in chromosome_dict.values():
        for feature in chromosome.feature_array:
            if len(feature.tags) == 0:
                # no tag -> default white color
                color = '#ffffff'
            else:
                # Get color with tag
                color = color_dict[tuple(feature.tags)]
            # Genieous format: /Color="#ff0000"
            # Old color (if any) will be replaced
            feature.set_attribute(key='Color', val=color)


def _compare(feature1, feature2):
    """
    Args:
        feature1: GenericFeature object

        feature2: GenericFeature object

    Returns: float
    """
    X, Y = feature1, feature2

    if X is None:
        return score['gap']

    elif Y is None:
        return score['gap']

    elif (X.type == 'CDS' and Y.type == 'intergenic') or \
            (X.type == 'intergenic' and Y.type == 'CDS'):
        return score['CDS_intergenic']

    elif X.type == 'CDS' and Y.type == 'CDS':
        if X.strand == Y.strand:
            if _has_common_tag(X.tags, Y.tags):
                return score['common_tag']
            else:
                return score['CDS_CDS_same_strand']
        elif X.strand != Y.strand:
            return score['CDS_CDS_opposite_strand']

    return score['else']


def _align(feature_array_1, feature_array_2):
    """
    Aligns two FeatureArray objects and return the optimal score of alignment

    Args:
        feature_array_1: FeatureArray object

        feature_array_2: FeatureArray object

    Returns: float
    """
    L1 = feature_array_1
    L2 = feature_array_2

    M = np.zeros((len(L1)+1, len(L2)+1), dtype=np.int)

    for r in range(len(L1)):
        M[r+1, 0] = M[r, 0] + score['gap']

    for c in range(len(L2)):
        M[0, c+1] = M[0, c] + score['gap']

    for r in range(len(L1)):
        for c in range(len(L2)):
            down_right = M[r, c] + _compare(L1[r], L2[c])
            right = M[r+1, c] + _compare(None, L2[c])
            down = M[r, c+1] + _compare(L1[r], None)

            M[r+1, c+1] = max(down_right, right, down)

    return M[-1, -1]


def _insert_intergenic(feature_array, seqname):
    """
    Args:
        feature_array: FeatureArray object

        seqname: str
    """
    L = []
    last_pos = 0
    for f in feature_array:
        if f.start > last_pos + 1:
            L.append(
                GenericFeature(
                    seqname=seqname,
                    type_='intergenic',
                    start=last_pos + 1,
                    end=f.start - 1,
                    strand='+',
                )
            )
            L.append(f)
        else:  # feature.start <= last_pos + 1
            L.append(f)

        last_pos = f.end
    return L


def get_qualifier(genbank, keywords):
    """
    From a <genbank> file, get the qualifiers that contain <keywords>

    Args:
        genbank: str, path-like
            The input genbank file

        keywords: str, or list of str

    Returns: list of str
    """
    if type(keywords) is str:
        keywords = [keywords]

    qualifiers = {}
    for chromosome in read_genbank(genbank):
        for feature in chromosome.feature_array:
            for qualifier, content in feature.attributes:
                # No need to look for keywords if content is not str
                if type(content) is not str: continue
                for word in keywords:
                    if word in content:
                        qualifiers.setdefault(qualifier, True)
    return list(qualifiers.keys())


def synteny(genbank=None, fasta=None, gtf=None, output=None, ax=None, qualifiers=('name', 'gene', 'label')):
    """
    Args:
        genbank: str, path-like

        fasta: str, path-like
            fasta and gtf together overrides genbank

        gtf: str, path-like

        output: str, path-like
            The output genbank file
            Has to be genbank file because the order of chromosomes is
                the only thing that matters for graphic display

        ax: matplotlib Axes instance

        qualifiers: str or list of str
            Qualifiers (e.g. gene, db_xref) of the genbank feature to be used
                to tag and align features

    Returns: list of str
        The order of seqname after synteny alignment and clustering
    """
    if fasta is not None and gtf is not None:
        chromosome_dict = read_fasta_gtf(fasta, gtf, as_dict=True, circular=False)
    else:
        chromosome_dict = read_genbank(genbank, as_dict=True)

    # For each feature, tag it with its own attributes using the <qualifiers>
    if type(qualifiers) is str:
        qualifiers = [qualifiers]
    for seqname, chromosome in chromosome_dict.items():
        for feature in chromosome.feature_array:
            _tag_feature(feature, qualifiers)

    # Generate a color dictionary from all features of all chromosomes
    # Each unique tag or combination of tags gets a unique color
    # Tags (keys) are hashed by tuple, colors are represented by hexadecimal strings, e.g.
    # {
    #   ('zinc finger', 'DNA binding'): #0f0f0f
    # }
    color_dict = _make_color_dict(chromosome_dict)

    # Set the color of each feature using the color_dict
    _set_feature_color(chromosome_dict, color_dict)

    # simply get a list of seqnames
    seqname_list = list(chromosome_dict.keys())

    # scores is in the format of condensed distance matrix
    #   which is an 1-D array of seqname pairs following the order given by combinations()
    scores = []
    for seqname1, seqname2 in combinations(seqname_list, 2):

        arr1 = chromosome_dict[seqname1].feature_array
        arr2 = chromosome_dict[seqname2].feature_array

        arr1 = _insert_intergenic(arr1, seqname1)
        arr2 = _insert_intergenic(arr2, seqname2)

        scores.append(_align(arr1, arr2))  # compute alignment score

    # Condensed distance matrix, which is literally a list of distance (must be positive)
    # distance = e^(-score), the higher the score, the shorter the distance
    dist_mat = np.exp(
        -np.array(scores) * score_scalar
    )

    # Compute linkage matrix
    linkage_matrix = hierarchy.linkage(dist_mat, 'average')

    # Plot dendrogram and get id order
    id_order = hierarchy.dendrogram(
        Z=linkage_matrix,
        ax=ax,
        orientation='right' # root on the right
    )['ivl'] # 'ivl': a list of labels corresponding to the leaf nodes

    # Get actual seqname order with id order
    seqname_order = [seqname_list[int(i)] for i in id_order]

    # Set the contig labels on the y axis
    if ax is not None:
        ax.set_yticklabels(seqname_order, rotation=90, verticalalignment='top')

    # Re-order Chromosome objects in a list
    data = [chromosome_dict[seqname] for seqname in seqname_order]

    if output is None: output = 'synteny.gb'
    write_genbank(data=data, file=output)

    return seqname_order
