from typing import Dict, List, Union
from .fasta import FastaParser
from .gff import read_gff
from .dataclass import Chromosome, FeatureArray, gff_to_generic_feature


def read_fasta_gff(
        fasta: str,
        gff: str,
        as_dict: bool = False,
        circular: bool = False) \
        -> Union[List[Chromosome], Dict[str, Chromosome]]:
    """
    Read fasta and gff simultaneously into Chromosome objects, i.e. annotated genomes

    Args:
        fasta: path-like

        gff: path-like

        as_dict:
            Return as dictionary or not

        circular:
            Shape of the DNA molecule

    Returns: list of Chromosome, or dict

        [Chromosome_1, Chromosome_2, ...]

        or

        {
            Chromosome_1.seqname: Chromosome_1,
            Chromosome_2.seqname: Chromosome_2, ...
        }
    """
    chromosomes = []

    feature_dict = read_gff(gff, as_dict=True)
    # feature_dict has the data structure:
    # {
    #   seqname: [GffFeature, ...], ...
    # }

    with FastaParser(fasta) as parser:
        for seqname, sequence in parser:

            gff_features = feature_dict.get(seqname, [])

            features = list(map(gff_to_generic_feature, gff_features))

            features = FeatureArray(
                seqname=seqname,
                chromosome_size=len(sequence),
                features=features,
                circular=circular)

            c = Chromosome(
                seqname=seqname,
                sequence=sequence,
                features=features,
                circular=False)

            chromosomes.append(c)

    if as_dict:
        return {chrom.seqname: chrom for chrom in chromosomes}

    return chromosomes
