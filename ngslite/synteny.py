import re
import random
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.cluster import hierarchy
from typing import Optional, Union, List, Dict
from .fasta_gtf import read_fasta_gtf
from .genbank_parse import read_genbank
from .genbank_write import write_genbank
from .dataclass import Chromosome, FeatureArray, GenericFeature


SCORE = {
    'gap': -1,
    'feature_intergenic': -1,
    'different_feature': -2,
    'same_feature_same_strand': +2,
    'same_feature_opposite_strand': -5,
    'common_tag': +100
}


SCORE_SCALAR = 0.01


class Aligner:

    def has_common_tag(
            self,
            tags1: List[str],
            tags2: List[str],
            skip: str = 'orf[0-9]+') -> bool:
        """
        skip: regex pattern
            Tags matching this pattern will be skipped
            Usually we want to skip very generic patterns like 'orf#####'
        """
        for t in tags1:
            if t in tags2:
                if not re.search(skip, t):
                    return True
        return False

    def compare(
            self,
            feature1: Optional[GenericFeature],
            feature2: Optional[GenericFeature]) -> float:

        X, Y = feature1, feature2

        if X is None:
            return SCORE['gap']

        if Y is None:
            return SCORE['gap']

        if X.type == Y.type == 'intergenic':
            return 0

        if X.type == 'intergenic' or Y.type == 'intergenic':
            return SCORE['feature_intergenic']

        if X.type != Y.type:
            return SCORE['different_feature']

        else:  # X.type == Y.type but not 'intergenic'

            if X.strand == Y.strand:

                if self.has_common_tag(X.tags, Y.tags):
                    return SCORE['common_tag']
                else:
                    return SCORE['same_feature_same_strand']

            else:  # X.strand != Y.strand
                return SCORE['same_feature_opposite_strand']

    def run(self,
            arr1: FeatureArray,
            arr2: FeatureArray) -> float:

        L1 = arr1
        L2 = arr2

        M = np.zeros((len(L1) + 1, len(L2) + 1), dtype=np.int)

        for r in range(len(L1)):
            M[r + 1, 0] = M[r, 0] + SCORE['gap']

        for c in range(len(L2)):
            M[0, c + 1] = M[0, c] + SCORE['gap']

        for r in range(len(L1)):
            for c in range(len(L2)):
                down_right = M[r, c] + self.compare(L1[r], L2[c])
                right = M[r + 1, c] + self.compare(None, L2[c])
                down = M[r, c + 1] + self.compare(L1[r], None)

                M[r + 1, c + 1] = max(down_right, right, down)

        return M[-1, -1]


class Synteny:

    genbank: Optional[str] = None
    fasta: Optional[str] = None
    gtf: Optional[str] = None
    output: str = 'synteny.gb'
    ax: Optional[plt.Axes] = None
    qualifiers: List[str] = ['name', 'gene', 'label']

    chromosome_dict: Dict[str, Chromosome]

    skip_pattern: str = 'orf[0-9]+'
    default_color = '#ffffff'

    def set_output(self, output: Optional[str]):
        if output is not None:
            self.output = output

    def set_qualifiers(
            self,
            qualifiers: Optional[Union[str, List[str]]]):

        if qualifiers is None:
            return

        if type(self.qualifiers) is str:
            self.qualifiers = [qualifiers]
        else:
            self.qualifiers = qualifiers

    def set_chromosome_dict(self):

        if self.fasta is not None and self.gtf is not None:
            self.chromosome_dict = read_fasta_gtf(
                fasta=self.fasta,
                gtf=self.gtf,
                as_dict=True,
                circular=False)
        else:
            self.chromosome_dict = read_genbank(
                file=self.genbank,
                as_dict=True)

    def __tag(self, feature: GenericFeature):

        for key, val in feature.attributes:
            if key in self.qualifiers:
                for word in re.split('[|]|[;]|[,]', str(val)):
                    if not re.search(self.skip_pattern, word):
                        feature.tags.append(word.strip())

        feature.tags = list(set(feature.tags))

    def tag_features(self):
        for seqname, chromosome in self.chromosome_dict.items():
            for feature in chromosome.feature_array:
                self.__tag(feature=feature)

    def __get_color_dict(self) -> Dict[tuple, str]:

        random.seed(1)

        color_dict = {}

        for chromosome in self.chromosome_dict.values():
            for feature in chromosome.feature_array:
                if len(feature.tags) == 0:  # no tag -> skip
                    continue
                key = tuple(set(feature.tags))
                color = f'#{random.randint(0, 256 ** 3):06x}'
                color_dict[key] = color

        return color_dict

    def set_feature_color(self):

        color_dict = self.__get_color_dict()

        for chromosome in self.chromosome_dict.values():
            for feature in chromosome.features:
                key = tuple(set(feature.tags))
                color = color_dict.get(key, self.default_color)
                feature.set_attribute(key='Color', val=color)  # Genieous format: /Color="#ff0000"

    def __insert_intergenic(
            self, feature_array: FeatureArray) -> FeatureArray:

        features = []

        last_pos = 0
        for f in feature_array:
            intergenic_gap = f.start > last_pos + 1
            if intergenic_gap:
                intergenic_feature = GenericFeature(
                    seqname=feature_array.seqname,
                    type_='intergenic',
                    start=last_pos + 1,
                    end=f.start - 1,
                    strand='+')
                features.append(intergenic_feature)
                features.append(f)
            else:
                features.append(f)

            last_pos = f.end

        return FeatureArray(
            features=features,
            seqname=feature_array.seqname,
            chromosome_size=feature_array.chromosome_size,
            circular=False)

    def get_score_matrix(self) -> np.ndarray:

        d = self.chromosome_dict

        seqnames = list(d.keys())

        # score_matrix is in the format of condensed distance matrix
        #   which is an 1-D array of seqname pairs following the order given by combinations()
        score_matrix = []

        for seqname1, seqname2 in combinations(seqnames, 2):
            chromosome1 = d[seqname1]
            chromosome2 = d[seqname2]

            arr1 = chromosome1.feature_array
            arr2 = chromosome2.feature_array

            arr1 = self.__insert_intergenic(feature_array=arr1)
            arr2 = self.__insert_intergenic(feature_array=arr2)

            score = Aligner().run(arr1=arr1, arr2=arr2)

            score_matrix.append(score)

        return np.array(score_matrix)

    def get_distance_matrix(
            self, score_matrix: np.ndarray) -> np.ndarray:
        # Distance must > 0
        # Define distance = e^(-score), the higher the score, the shorter the distance
        return np.exp(- score_matrix * SCORE_SCALAR)

    def plot_tree_and_get_idx_order(
            self, linkage_matrix: np.ndarray) -> List[int]:

        idx_order = hierarchy.dendrogram(
            Z=linkage_matrix,
            ax=self.ax,
            orientation='top'  # root on top
        )['ivl']  # 'ivl': a list of labels corresponding to the leaf nodes

        return idx_order

    def get_ordered_seqnames(
            self, idx_order: List[int]) -> List[str]:

        seqnames = list(self.chromosome_dict.keys())

        return [seqnames[int(i)] for i in idx_order]

    def get_ordered_chromosomes(
            self, ordered_seqnames: List[str]) -> List[Chromosome]:

        ordered_chromosomes = [
            self.chromosome_dict[seqname] for seqname in ordered_seqnames]

        return ordered_chromosomes

    def set_x_axis_contig_labels(self, labels: List[str]):
        if self.ax is not None:
            self.ax.set_xticklabels(
                labels,
                rotation=90,
                horizontalalignment='right')

    def main(
            self,
            genbank: Optional[str] = None,
            fasta: Optional[str] = None,
            gtf: Optional[str] = None,
            output: Optional[str] = None,
            ax: Optional[plt.Axes] = None,
            qualifiers: Optional[Union[str, List[str]]] = None) -> List[str]:

        self.genbank = genbank
        self.fasta = fasta
        self.gtf = gtf
        self.set_output(output=output)
        self.ax = ax
        self.set_qualifiers(qualifiers=qualifiers)

        self.set_chromosome_dict()

        self.tag_features()

        self.set_feature_color()

        score_matrix = self.get_score_matrix()

        distance_matrix = self.get_distance_matrix(score_matrix=score_matrix)

        linkage_matrix = hierarchy.linkage(distance_matrix, 'average')

        idx_order = self.plot_tree_and_get_idx_order(
            linkage_matrix=linkage_matrix)

        ordered_seqnames = self.get_ordered_seqnames(idx_order=idx_order)

        ordered_chromosomes = self.get_ordered_chromosomes(
            ordered_seqnames=ordered_seqnames)

        self.set_x_axis_contig_labels(labels=ordered_seqnames)

        write_genbank(data=ordered_chromosomes, file=self.output)

        return ordered_seqnames


def synteny(
        genbank: Optional[str] = None,
        fasta: Optional[str] = None,
        gtf: Optional[str] = None,
        output: Optional[str] = None,
        ax: Optional[plt.Axes] = None,
        qualifiers: Optional[Union[str, List[str]]] = None) -> List[str]:
    """
    Args:
        genbank: path-like

        fasta: path-like
            fasta and gtf together overrides genbank

        gtf: path-like

        output: path-like
            The output genbank file
            Has to be genbank file because the order of chromosomes is
                the only thing that matters for graphic display

        ax: matplotlib Axes instance

        qualifiers:
            Qualifiers (e.g. gene, db_xref) of the genbank feature to be used
                to tag and align features

    Returns:
        The order of seqname after synteny alignment and clustering
    """
    return Synteny().main(
        genbank=genbank,
        fasta=fasta,
        gtf=gtf,
        output=output,
        ax=ax,
        qualifiers=qualifiers)


def get_qualifiers(
        genbank: str,
        keywords: Union[str, List[str]]) -> List[str]:
    """
    From a <genbank> file, get the qualifiers that contain <keywords>

    Args:
        genbank: path-like
            The input genbank file

        keywords
    """
    if type(keywords) is str:
        keywords = [keywords]

    qualifiers = {}
    for chromosome in read_genbank(genbank):
        for feature in chromosome.feature_array:
            for qualifier, content in feature.attributes:
                # No need to look for keywords if content is not str
                if type(content) is not str:
                    continue
                for word in keywords:
                    if word in content:
                        qualifiers.setdefault(qualifier, True)
    return list(qualifiers.keys())
