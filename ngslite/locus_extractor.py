from.data_class import FeatureArray, Chromosome
from .fasta import read_fasta, write_fasta
from .gtftools import read_gtf, write_gtf
from .dnatools import rev_comp
from .genbank_parse import read_genbank
from .genbank_write import write_genbank


def _find_feature(feature_arr, keywords):
    """
    In a FeatureArray object, if any one of the GenericFeature contains any one of the <keywords>,
        then return the GenericFeature object

    Args:
        feature_arr: FeatureArray object

        keywords: list of str
    """
    # For each GenericFeature
    for f in feature_arr:
        # For each keyword
        for word in keywords:
            # For each value in the attributes [(key_1, val_1), (key_2, val_2), ...]
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
    Similar to locus_extractor() except the genomic sequence and annotation are
        read from a genbank file

    Args:
        genbank: str, path-like
            The input genbank file

        keywords: str, or list of str

        flank: int

        output: str, path-like
            The output genbank file
    """
    if type(keywords) is str:
        keywords = [keywords]

    loci = []  # list of Chromosome objects

    # Iterate through each contig
    for chromosome in read_genbank(genbank):  # read_genbank() -> list of Chromosome objects
        seqname = chromosome.seqname
        feature_array = chromosome.feature_array
        sequence = chromosome.sequence

        seed = _find_feature(feature_array, keywords)
        if seed is None: continue  # seed not found, skip this chromosome

        region_start = seed.start - flank
        region_end = seed.end + flank

        region_start = max(region_start, 1)  # out of bound
        region_end = min(region_end, len(sequence))

        # Add information to the seqname
        seqname = seqname + f';from={region_start};to={region_end};strand={seed.strand}'

        # Extract sequence from fasta
        sequence = sequence[region_start-1: region_end]

        # Crop features from the feature array
        feature_array.crop(region_start, region_end)

        # Reverse sequence and features if seed is on the reverse strand
        if seed.strand == '-':
            sequence = rev_comp(sequence)
            feature_array.reverse()

        loci.append(
            Chromosome(
                seqname=seqname,
                sequence=sequence,
                feature_array=feature_array,
                circular=False  # Must be linear molecule after locus extraction
            )
        )

    # Write output files
    write_genbank(data=loci, file=output)