from .fasta import FastaParser
from .gtftools import read_gtf
from .data_class import Chromosome, FeatureArray


def read_fasta_gtf(fasta, gtf, as_dict=False, circular=False):
    """
    Read fasta and gtf simultaneously into Chromosome ojbects, i.e. annotated genomes

    Args:
        fasta: str, path-like

        gtf: str, path-like

        as_dict: bool

        circular: bool
            Shape of the DNA molecule

    Returns: list of Chromosome objects, or dict

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
            feature_array = FeatureArray(
                seqname=seqname,
                genome_size=len(sequence),
                features=gtf_features,
                circular=circular
            )
            chromosomes.append(
                Chromosome(
                    seqname=seqname,
                    sequence=sequence,
                    feature_array=feature_array,
                    circular=False
                )
            )
    if as_dict:
        return {chrom.seqname: chrom for chrom in chromosomes}
    return chromosomes
