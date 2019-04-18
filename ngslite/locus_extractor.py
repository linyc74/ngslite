from .fasta import read_fasta, FastaWriter
from .gtftools import gtf_to_dict, dict_to_gtf
from .dnatools import rev_comp
from .genbank_parse import genbank_to_fasta, genbank_to_gtf
from .genbank_write import make_genbank
from .lowlevel import __temp
import os


def __find_feature(feature_arr, keywords):
    """
    In an array of GtfFeature objects, find one that contains <keywords> in GtfFeature.attribute

    Args:
        feature_arr: list [GtfFeature, ...]

        keywords: list [str, ...]

    Returns: GtfFeature
    """
    for f in feature_arr:
        for k in keywords:
            if k in f.attribute:
                return f
    return None


def __subset_features(feature_arr, start, end):
    """
    In an array of GtfFeature objects, return those located within <start> to <end>

    Args:
        feature_arr: list [GtfFeature, ...]

        start: int
            1-based, inclusive

        end: int
            1-based, inclusive

    Returns:
        list [GtfFeature, ...]
    """
    __feature_arr = []
    for f in feature_arr:
        if f.start >= start and f.end <= end:
            __feature_arr.append(f)
    return __feature_arr


def __move_features(feature_arr, offset):
    """
    For each GtfFeature objects in an array,
        GtfFeature.start += offset
        GtfFeature.end += offset

    Args:
        feature_arr: list [GtfFeature, ...]

        offset: int

    Returns:
        list [GtfFeature, ...]
    """
    __feature_arr = []
    for f in feature_arr:
        __f = f._replace(start=f.start + offset, end=f.end + offset)
        __feature_arr.append(__f)
    return __feature_arr


def __crop_features(feature_arr, start, end):
    """
    From an array of GtfFeature objects, crop out those located within <start> to <end>

    Move the locations of GtfFeature by
        GtfFeature.start += -start + 1
        GtfFeature.end += -start + 1

    Args:
        feature_arr: list [GtfFeature, ...]

        start: int
            1-based, inclusive

        end: int
            1-base, inclusive
    """
    __feature_arr = __subset_features(feature_arr, start, end)
    return __move_features(__feature_arr, offset=-start+1)


def __reverse_features(feature_arr, length):
    """
    In an array of GtfFeature objects, reverse the order of those objects and
        invert start and end of each object by
            start = length - end + 1
            end = length - start + 1

    Args:
        feature_arr: list [GtfFeature, ...]

        length: int,
            The length (bp) of the genome
    """
    __feature_arr = []
    for f in feature_arr[::-1]:
        __feature_arr.append(
            f._replace(
                start = length - f.end + 1,
                end = length - f.start + 1,
                strand = ['+', '-'][f.strand == '+']
            )
        )
    return __feature_arr


def locus_extractor(fasta, gtf, keywords, flank, fasta_out, gtf_out):
    """
    Use <keywords> to find a 'seed' feature (i.e. CDS) in each contig and
        extract the genomic region from (seed.start - flank) to (seed.end + flank)

    If the seed feature is on the reverse strand,
        reverse the genomic region

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

    fasta = read_fasta(fasta, as_dict=True)
    gtf = gtf_to_dict(gtf)

    __fasta = dict()
    __gtf = dict()

    # Iterate through each contig
    for seqname in fasta.keys():
        feature_arr = gtf.get(seqname, None)

        if feature_arr is None:  # no annotation for the contig
            continue

        seed = __find_feature(feature_arr, keywords)

        if seed is None:  # seed not found, skip this contig
            continue

        region_start = seed.start - flank
        region_end = seed.end + flank

        region_start = max(region_start, 1)  # out of bound
        region_end = min(region_end, len(fasta[seqname]))

        # Extract sequence from fasta
        __fasta[seqname] = fasta[seqname][region_start-1: region_end]

        # Crop features from the feature array
        __gtf[seqname] = __crop_features(
            feature_arr=gtf[seqname],
            start=region_start,
            end=region_end
        )

        # Reverse sequence and features if seed is on the reverse strand
        if seed.strand == '-':
            __fasta[seqname] = rev_comp(__fasta[seqname])
            __gtf[seqname] = __reverse_features(
                feature_arr=__gtf[seqname],
                length=len(__fasta[seqname])
            )

    with FastaWriter(fasta_out) as writer:
        for head, seq in __fasta.items():
            writer.write(head, seq)

    dict_to_gtf(dict_=__gtf, output=gtf_out)


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
        skip_attributes=['translation', 'codon_start']
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
