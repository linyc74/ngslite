from copy import deepcopy
from typing import List, Tuple, Union
from .fasta import write_fasta
from .dnatools import rev_comp
from .gtftools import write_gtf
from .gfftools import write_gff
from .fasta_gtf import read_fasta_gtf
from .fasta_gff import read_fasta_gff
from .genbank_parse import read_genbank
from .genbank_write import write_genbank
from .dataclass import GenericFeature, FeatureArray, Chromosome


class LocusExtractor:

    keywords: List[str]
    flank: int

    def set_keywords(self, keywords: Union[str, List[str]]):
        self.keywords = [keywords] if type(keywords) is str else keywords

    def get_seed_features(
            self,
            features: FeatureArray) -> List[GenericFeature]:

        ret = []

        for f in features:
            for word in self.keywords:
                for key, val in f.attributes:
                    if word in str(val):
                        ret.append(f)

        return ret

    def get_region_start_end(
            self,
            chromosome: Chromosome,
            seed_feature: GenericFeature) -> Tuple[int, int]:

        s = seed_feature

        region_start = max(s.start - self.flank, 1)
        region_end = min(s.end + self.flank, len(chromosome.sequence))

        return region_start, region_end

    def get_locus_seqname(
            self,
            chromosome: Chromosome,
            seed_feature: GenericFeature) -> str:

        region_start, region_end = self.get_region_start_end(
            chromosome=chromosome,
            seed_feature=seed_feature)

        return f'{chromosome.seqname}:{region_start}-{region_end}'

    def get_locus_sequence(
            self,
            chromosome: Chromosome,
            seed_feature: GenericFeature) -> str:

        region_start, region_end = self.get_region_start_end(
            chromosome=chromosome,
            seed_feature=seed_feature)

        seq = chromosome.sequence[region_start - 1: region_end]

        if seed_feature.strand == '-':
            seq = rev_comp(seq)

        return seq

    def __set_locus_seqname(
            self,
            features: FeatureArray,
            locus_seqname: str):
        features.seqname = locus_seqname
        for f in features:
            f.seqname = locus_seqname

    def get_locus_features(
            self,
            chromosome: Chromosome,
            seed_feature: GenericFeature,
            locus_seqname: str) -> FeatureArray:

        region_start, region_end = self.get_region_start_end(
            chromosome=chromosome,
            seed_feature=seed_feature)

        features = deepcopy(chromosome.features)
        features.crop(start=region_start, end=region_end)

        self.__set_locus_seqname(
            features=features,
            locus_seqname=locus_seqname)

        if seed_feature.strand == '-':
            features.reverse()

        return features

    def get_loci(
            self,
            chromosomes: List[Chromosome]) -> List[Chromosome]:

        loci: List[Chromosome] = []

        for chromosome in chromosomes:

            seed_features = self.get_seed_features(
                features=chromosome.features)

            for seed_feature in seed_features:

                locus_seqname = self.get_locus_seqname(
                    chromosome=chromosome,
                    seed_feature=seed_feature)

                locus_sequence = self.get_locus_sequence(
                    chromosome=chromosome,
                    seed_feature=seed_feature)

                locus_features = self.get_locus_features(
                    chromosome=chromosome,
                    seed_feature=seed_feature,
                    locus_seqname=locus_seqname)

                locus = Chromosome(
                    seqname=locus_seqname,
                    sequence=locus_sequence,
                    features=locus_features,
                    circular=False)

                loci.append(locus)

        return loci

    def write_fasta(
            self, loci: List[Chromosome], file: str):
        data = {l.seqname: l.sequence for l in loci}
        write_fasta(data=data, file=file)

    def write_gtf(
            self, loci: List[Chromosome], file: str):
        data = {l.seqname: l.features for l in loci}
        write_gtf(data=data, file=file)

    def write_gff(
            self, loci: List[Chromosome], file: str):

        data = {l.seqname: l.features for l in loci}

        write_gff(data=data, file=file)

    def use_fasta_gtf(
            self,
            fasta: str,
            gtf: str,
            keywords: Union[str, List[str]],
            flank: int,
            fasta_out: str,
            gtf_out: str):

        self.set_keywords(keywords=keywords)
        self.flank = flank

        chromosomes = read_fasta_gtf(
            fasta=fasta,
            gtf=gtf,
            as_dict=False,
            circular=False)

        loci = self.get_loci(chromosomes=chromosomes)

        self.write_fasta(loci=loci, file=fasta_out)
        self.write_gtf(loci=loci, file=gtf_out)

    def use_fasta_gff(
            self,
            fasta: str,
            gff: str,
            keywords: Union[str, List[str]],
            flank: int,
            fasta_out: str,
            gff_out: str):

        self.set_keywords(keywords=keywords)
        self.flank = flank

        chromosomes = read_fasta_gff(
            fasta=fasta,
            gff=gff,
            as_dict=False,
            circular=False)

        loci = self.get_loci(chromosomes=chromosomes)

        self.write_fasta(loci=loci, file=fasta_out)
        self.write_gff(loci=loci, file=gff_out)

    def use_gbk(
            self,
            gbk: str,
            keywords: Union[str, List[str]],
            flank: int,
            output: str):

        self.set_keywords(keywords=keywords)
        self.flank = flank

        chromosomes = read_genbank(file=gbk)

        loci = self.get_loci(chromosomes=chromosomes)

        write_genbank(data=loci, file=output)


def locus_extractor(
        fasta: str,
        gtf: str,
        keywords: Union[str, List[str]],
        flank: int,
        fasta_out: str,
        gtf_out: str):
    """
    Use <keywords> to find a 'seed' feature (i.e. CDS) in each contig and
        extract the genomic region from (seed.start - flank) to (seed.end + flank)

    If the seed feature is on the reverse strand,
        reverse the genomic region

    Does not support circular genome, thus if a GtfFeature has start > end position, then
        that feature will be discarded

    Args:
        fasta: path-like
            The input fasta file

        gtf: path-like
            The input GTF file

        keywords

        flank

        fasta_out: path-like
            The output fasta file

        gtf_out: path-like
            The output GTF file
    """
    LocusExtractor().use_fasta_gtf(
        fasta=fasta,
        gtf=gtf,
        keywords=keywords,
        flank=flank,
        fasta_out=fasta_out,
        gtf_out=gtf_out)


def genbank_locus_extractor(
        genbank: str,
        keywords: Union[str, List[str]],
        flank: int,
        output: str):
    """
    Similar to locus_extractor() except the genomic sequence and annotation are
        read from a genbank file

    Args:
        genbank: path-like
            The input genbank file

        keywords

        flank

        output: path-like
            The output genbank file
    """
    LocusExtractor().use_gbk(
        gbk=genbank,
        keywords=keywords,
        flank=flank,
        output=output)
