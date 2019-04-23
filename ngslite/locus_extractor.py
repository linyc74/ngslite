from.data_class import FeatureArray
from .fasta import read_fasta, write_fasta
from .gtftools import read_gtf, write_gtf
from .dnatools import rev_comp
from .genbank_parse import genbank_to_fasta, genbank_to_gtf
from .genbank_write import make_genbank
from .lowlevel import __temp
import os


def _find_feature(feature_arr, keywords):
    """
    Args:
        feature_arr: list of GenericFeature objects

        keywords: list of str
    """
    for f in feature_arr:
        for word in keywords:
            for key, val in f.attributes:
                if word in str(val):
                    return f
    return None


def locus_extractor(fasta, gtf, keywords, flank, fasta_out, gtf_out):
    """
    Use <keywords> to find a 'seed' feature (i.e. CDS) in each contig and
        extract the genomic region from (seed.start - flank) to (seed.end + flank)

    If the seed feature is on the reverse strand,
        reverse the genomic region

    Does not support circular genome, thus if a GtfFeature has start > end position, then
        that feature will be discarded

    Args:
        fasta: str, path-like
            The input fasta file

        gtf: str, path-like
            The input GTF file

        keywords: list of str

        flank: int

        fasta_out: str, path-like
            The output fasta file

        gtf_out: str, path-like
            The output GTF file
    """
    if type(keywords) is str:
        keywords = [keywords]

    fasta_dict = read_fasta(fasta, as_dict=True)
    gtf_dict = read_gtf(gtf, as_dict=True)

    __sequence_dict = dict()
    __feature_dict = dict()

    # Iterate through each contig
    for seqname, sequence in fasta_dict.items():
        gtf_features = gtf_dict.get(seqname, [])
        if len(gtf_features) == 0: continue  # no annotation, skip this contig

        feature_array = FeatureArray(
            seqname=seqname,
            genome_size=len(sequence),
            features=gtf_features,
            circular=False
        )

        seed = _find_feature(feature_array, keywords)
        if seed is None: continue  # seed not found, skip this contig

        region_start = seed.start - flank
        region_end = seed.end + flank

        region_start = max(region_start, 1)  # out of bound

        # Extract sequence from fasta
        sequence = sequence[region_start-1: region_end]

        # Crop features from the feature array
        feature_array.crop(region_start, region_end)

        # Reverse sequence and features if seed is on the reverse strand
        if seed.strand == '-':
            sequence = rev_comp(sequence)
            feature_array.reverse()

        __sequence_dict[seqname] = sequence
        __feature_dict[seqname] = feature_array

    # Write output files
    write_fasta(data=__sequence_dict, file=fasta_out)
    write_gtf(data=__feature_dict, file=gtf_out)


def genbank_locus_extractor(genbank, keywords, flank, output):
    """
    A wrapper function of locus_extractor() that does:
        genbank to (fasta, gtf)
        locus_extractor()
        (fasta, gtf) to genbank

    Args:
        genbank: str, path-like
            The input genbank file

        keywords: list of str

        flank: int

        output: str, path-like
            The output genbank file
    """
    temp_fa = __temp('temp', '.fa')
    temp_gtf = __temp('temp', '.gtf')

    genbank_to_fasta(
        file=genbank,
        output=temp_fa
    )

    genbank_to_gtf(
        file=genbank,
        output=temp_gtf,
        skip_attributes=['translation', 'codon_start', 'transl_table']
    )

    locus_extractor(
        fasta=temp_fa,
        gtf=temp_gtf,
        keywords=keywords,
        flank=flank,
        fasta_out=temp_fa,
        gtf_out=temp_gtf
    )

    make_genbank(
        fasta=temp_fa,
        gtf=temp_gtf,
        output=output
    )

    os.remove(temp_fa)
    os.remove(temp_gtf)
