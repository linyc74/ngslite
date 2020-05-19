from typing import Dict, List, Union
from .fasta import FastaParser
from .gtftools import read_gtf
from .dataclass import Chromosome, FeatureArray


def read_fasta_gtf(
        fasta: str,
        gtf: str,
        as_dict: bool = False,
        circular: bool = False) \
        -> Union[List[Chromosome], Dict[str, Chromosome]]:
    """
    Read fasta and gtf simultaneously into Chromosome objects, i.e. annotated genomes

    Args:
        fasta: path-like

        gtf: path-like

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

    feature_dict = read_gtf(gtf, as_dict=True)
    # feature_dict has the data structure:
    # {
    #   seqname: [GtfFeature, ...], ...
    # }

    with FastaParser(fasta) as parser:
        for seqname, sequence in parser:

            gtf_features = feature_dict.get(seqname, [])

            features = FeatureArray(
                seqname=seqname,
                genome_size=len(sequence),
                features=gtf_features,
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
