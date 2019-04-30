from .fasta_gtf import read_fasta_gtf
from .genbank_parse import read_genbank
from .genbank_write import write_genbank
from .data_class import GenericFeature
import numpy as np
import re
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


def _tag_feature(feature, attr_keys, skip='orf[0-9]+'):
    """
    Args:
        feature: GenericFeature object

        attr_keys: str, or list of str

        skip: str, regex pattern
            Attribute words matching this pattern will be skipped
            Usually we want to skip very generic patterns like 'orf#####'
    """
    for key, val in feature.attributes:
        if key in attr_keys:
            for word in re.split('[|]|[;]|[,]', str(val)):
                if not re.search(skip, word):
                    feature.tags.append(word.strip())


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


def synteny(genbank=None, fasta=None, gtf=None, output=None, ax=None):
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

    Returns: list of str
        The order of seqname after synteny alignment and clustering
    """
    if fasta is not None and gtf is not None:
        chromosome_dict = read_fasta_gtf(fasta, gtf, as_dict=True, circular=False)
    else:
        chromosome_dict = read_genbank(genbank, as_dict=True)

    # For each feature, tag it with its own attributes
    attr_keys = ['name', 'gene']
    for seqname, chromosome in chromosome_dict.items():
        for feature in chromosome.feature_array:
            _tag_feature(feature, attr_keys)

    # Simply a list of seqnames
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

    linkage_matrix = hierarchy.linkage(dist_mat, 'average')
    id_order = hierarchy.dendrogram(linkage_matrix, ax=ax)['ivl']
    seqname_order = [seqname_list[int(i)] for i in id_order]

    # Re-order Chromosome objects in a list
    data = [chromosome_dict[seqname] for seqname in seqname_order]

    if output is None: output = 'synteny.gb'
    write_genbank(data=data, file=output)

    return seqname_order
