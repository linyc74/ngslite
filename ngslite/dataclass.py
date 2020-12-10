from copy import deepcopy
from typing import List, Tuple, Optional, Union
from collections import namedtuple
from .dna import rev_comp


GtfFeature = namedtuple('GtfFeature', 'seqname source feature start end score strand frame attribute')


GffFeature = namedtuple('GffFeature', 'seqid source type start end score strand phase attributes')


class GenericFeature:
    """
    Instance attributes:
        seqname:
            The unique name of chromosome, genome or contig to which this feature belongs

        type:
            The type of the feature, e.g. 'CDS', 'tRNA', 'misc_feature'

        start:
            1-based inclusive

        end:
            1-based inclusive

        strand:
            {'+' or '-'}

        frame:
            {1, 2, 3}

        attributes:
            [
                ('gene', 'dam'),
                ('locus_tag', 'T4p032')
                ('note', 'produces N6-methyladenine in GATC sites (in modified or non-modified cytosine DNA)'),
                ('codon_start', 1)
            ]

        tags

        regions:
            Default just one region

        partial_start
            Partial 5' end

        partial_end
            Partial 3' end
    """

    def __init__(
            self,
            seqname: str,
            type_: str,
            start: int,
            end: int,
            strand: str,
            attributes: Optional[List[Tuple[str, Union[str, int, float]]]] = None,
            tags: Optional[List[str]] = None,
            regions: Optional[List[Tuple[int, int, str]]] = None,
            frame: int = 1,
            partial_start: bool = False,
            partial_end: bool = False,
            chromosome_size: Optional[int] = None):

        self.seqname = seqname
        self.type = type_
        self.start = start
        self.end = end
        self.strand = strand
        self.frame = frame
        self.partial_start = partial_start
        self.partial_end = partial_end
        self.chromosome_size = chromosome_size

        if attributes is None:
            attributes = []
        self.attributes = attributes

        if tags is None:
            tags = []
        self.tags = tags

        if regions is None:
            regions = [(self.start, self.end, self.strand)]
        self.regions = regions

    def __len__(self):
        if self.end >= self.start:
            return self.end - self.start + 1
        else:
            assert self.chromosome_size is not None, \
                f'Cannot get len(feature) because end < start but chromosome_size is None'
            return self.chromosome_size - self.start + 1 + self.end

    def __str__(self):
        ps = ['', ' (partial)'][self.partial_start]
        pe = ['', ' (partial)'][self.partial_end]

        attr_str = ''
        for key, val in self.attributes:
            if type(val) is str:
                attr_str += f"\n    {key}: '{val}'"
            else:  # type(val) is int or float
                attr_str += f"\n    {key}: {val}"

        return f'''\
GenericFeature
  seqname: '{self.seqname}'
  type: '{self.type}'
  start: {self.start}{ps}
  end: {self.end}{pe}
  strand: '{self.strand}'
  frame: {self.frame}
  attributes: {attr_str}
  tags: {','.join(self.tags)}
  regions: {self.regions}'''

    def __repr__(self):
        return f"""GenericFeature(seqname='{self.seqname}', type_='{self.type}'\
, start={self.start}, end={self.end}, strand='{self.strand}', attributes=\
{self.attributes}, tags={self.tags}, regions={self.regions}, frame={self.frame}\
, partial_start={self.partial_start}, partial_end={self.partial_end})"""

    def get_attribute(self, key: str) -> \
            Optional[Union[str, int, float, List[Union[str, int, float]]]]:

        vals = [v for k, v in self.attributes if k == key]
        if len(vals) == 0:
            return None
        elif len(vals) == 1:
            return vals[0]
        else:
            return vals

    def set_attribute(self, key: str, val: Union[str, int, float]):
        """
        If more than one attributes found, then set the first one, and remove the rest
        If no attribute found, then add an attribute with key, val
        """
        new = []  # new attributes list
        done = False  # done setting the attribute or not

        for k, v in self.attributes:
            if k == key:
                if not done:
                    new.append((key, val))
                    done = True
                else:
                    # Already done, don't append the (key, val)
                    #   to the new list, i.e. remove the attribute
                    pass
            else:
                # Keep existing attributes where k != key
                new.append((k, v))

        if not done:
            # No key found, add a new attribute
            new.append((key, val))

        self.attributes = new

    def add_attribute(self, key: str, val: Union[str, int, float]):
        self.attributes.append((key, val))

    def remove_attribute(self, key: str):
        self.attributes = [(k, v) for (k, v) in self.attributes if k != key]


class FeatureArray:
    """
    A list of GenericFeature that behaves like a list, but more than a list
    """

    seqname: str
    chromosome_size: int
    circular: bool
    __pos: int

    # The main data structure, always remain sorted
    __arr: List[GenericFeature]

    def __init__(
            self,
            seqname: str,
            chromosome_size: int,
            features: Optional[Union[GenericFeature, List[GenericFeature]]] = None,
            circular: bool = False):

        self.seqname = seqname
        self.chromosome_size = chromosome_size
        self.circular = circular
        self.__pos = 0

        if features is None:
            features = []
        elif type(features) is GenericFeature:
            features = [features]

        self.__arr = features

        self.sort()

    def __getitem__(self, item: Union[int, slice]):
        if isinstance(item, int):
            return self.__arr[item]
        elif isinstance(item, slice):
            return FeatureArray(
                seqname=self.seqname,
                chromosome_size=self.chromosome_size,
                features=self.__arr[item],
                circular=self.circular)
        else:
            raise TypeError('FeatureArray indices must be integers or slices')

    def __add__(self, other):
        if type(other) is not FeatureArray:
            raise TypeError(f"unsupported operand type(s) for +: '{type(self)}' and '{type(other)}'")

        for attr in ['seqname', 'chromosome_size', 'circular']:
            assert getattr(self, attr) == getattr(other, attr), f"where attr = '{attr}'"

        return FeatureArray(
            seqname=self.seqname,
            chromosome_size=self.chromosome_size,
            features=self.__arr + other.__arr,
            circular=self.circular)

    def __iter__(self):
        self.__pos = 0
        return self

    def __next__(self):
        if self.__pos < len(self.__arr):
            r = self.__arr[self.__pos]
            self.__pos += 1
            return r
        else:
            self.__pos = 0
            raise StopIteration

    def __len__(self):
        return len(self.__arr)

    def __str__(self):
        shape = 'circular' if self.circular else 'linear'
        return f'''\
FeatureArray
  seqname: '{self.seqname}'
  genome size: {self.chromosome_size} bp ({shape})
  number of features: {len(self.__arr)}'''

    def __repr__(self):
        return f'''\
FeatureArray(seqname={self.seqname}, genome_size=\
{self.chromosome_size}, features={repr(self.__arr)}, circular={self.circular})'''

    def __wrap(self):
        """
        For circular genomes the feature range can be out of bound (< 1 or > genome_size)
        This method wraps those positions back between 1..genome_size
        """
        size = self.chromosome_size
        for f in self.__arr:
            if f.start > size:
                f.start -= size
            elif f.start < 1:
                f.start += size

            if f.end > size:
                f.end -= size
            elif f.end < 1:
                f.end += size

            for i in range(len(f.regions)):
                start, end, strand = f.regions[i]

                if start > size:
                    start -= size
                elif start < 1:
                    start += size

                if end > size:
                    end -= size
                elif end < 1:
                    end += size

                f.regions[i] = (start, end, strand)

    def sort(self):
        """
        Sort features by each feature's start position
        """
        self.__arr = sorted(self.__arr, key=lambda f: f.start)

    def append(self, feature: GenericFeature):
        """
        Append a new feature
        """
        self.__arr.append(feature)
        self.sort()

    def subset(self, start: int, end: int):
        """
        Subset (or filter) features without changing their positions

        Args:
            start:
                1-based, inclusive

            end:
                1-based, inclusive
        """

        if not self.circular and end < start:
            raise AssertionError('For linear molecule end position must be greater or equal than start position')

        # Circular range -> Set end to greater than genome_size
        if end < start:
            end += self.chromosome_size

        arr = []
        for f in self.__arr:
            fs, fe = f.start, f.end
            # Circular feature -> Set feature end to greater than genome_size
            if fe < fs:
                fe += self.chromosome_size
            # Main criteria to be satisfied
            if start <= fs and fe <= end:
                arr.append(f)
        self.__arr = arr

    def offset(self, offset: int):
        """
        Offset (shift) features by bp
        """
        size = self.chromosome_size
        for f in self.__arr:
            f.start += offset
            f.end += offset
            for i, region in enumerate(f.regions):
                start, end, strand = region
                start += offset
                end += offset
                f.regions[i] = (start, end, strand)

            if self.circular:
                self.__wrap()
            # After wrapping every feature will be in-bound

        in_bound = lambda x: 1 <= x.start <= size and 1 <= x.end <= size
        self.__arr = list(filter(in_bound, self.__arr))

    def crop(self, start: int, end: int):
        """
        Crop features from <start> to <end>,
            which resets positions of features and the genome size

        Args:
            start:
                1-based, inclusive

            end:
                1-based, inclusive
        """
        # Covers whole genome -> Do not subset
        if start == 1 and end == self.chromosome_size:
            return

        if not self.circular and end < start:
            return

        self.subset(start, end)
        self.offset(-start + 1)

        self.chromosome_size = min(end, self.chromosome_size) - max(start, 1) + 1

    def reverse(self):
        """
        Reverse features for reverse complementary of the genome sequence,
            which makes start as end, and end as start
        """
        size = self.chromosome_size
        for f in self.__arr:
            # Simultaneously update start and end
            f.start, f.end = (size - f.end + 1), (size - f.start + 1)
            f.strand = '-' if f.strand == '+' else '+'

            for i in range(len(f.regions)):
                start, end, strand = f.regions[i]
                # Simultaneously update start and end
                start, end = size - end + 1, size - start + 1
                strand = '-' if strand == '+' else '+'
                f.regions[i] = (start, end, strand)
            # Use start position to sort the regions (i.e. exons)
            f.regions = sorted(f.regions, key=lambda x: x[0])

        self.sort()

    def pop(self, index: int = -1):
        """
        Pop out a feature
        """
        return self.__arr.pop(index)

    def set_seqname(self, seqname: str):
        self.seqname = seqname
        for f in self.__arr:
            f.seqname = seqname

    def copy(self):
        return deepcopy(self)


class Chromosome:
    """
    Instance attributes:
        seqname:
            The unique name of chromosome, genome or contig to which this feature_array belons

        sequence:
            The genome sequence

        feature_array:
            The annotation

        circular:
            Is circular genome or not

        genbank_locus_text:
            The complete LOCUS section of a genbank file
    """

    seqname: str
    sequence: str
    features: FeatureArray
    feature_array: FeatureArray
    circular: bool
    genbank_locus_text: str

    def __init__(
            self,
            seqname: str,
            sequence: str,
            features: FeatureArray,
            circular: bool = False,
            genbank_locus_text: str = ''):

        self.seqname = seqname
        self.sequence = sequence
        self.features = features
        self.feature_array = self.features
        self.circular = circular
        self.genbank_locus_text = genbank_locus_text

    def __repr__(self):
        return f"""Chromosome(seqname='{self.seqname}', sequence='{self.sequence}'\
, features={repr(self.features)}, circular={self.circular}, genbank_locus_text=\
'{self.genbank_locus_text})'"""

    def crop(self, start: int, end: int):
        """
        Args:
            start: 1-based, inclusive

            end: 1-based, inclusive
        """

        if not self.circular and end < start:
            return

        if start < 1:
            start = 1

        self.sequence = self.sequence[start-1:end]
        self.features.crop(start=start, end=end)

    def reverse(self):
        """
        Reverse the features for reverse complementary of the genome sequence,
            which simply makes start as end, and end as start
        """
        self.sequence = rev_comp(self.sequence)
        self.features.reverse()

    def set_seqname(self, seqname: str):
        self.seqname = seqname
        self.features.set_seqname(seqname=seqname)
        
    def copy(self):
        return deepcopy(self)


def gtf_to_generic_feature(gtf_feature: GtfFeature) -> GenericFeature:

    assert type(gtf_feature) is GtfFeature

    f = gtf_feature

    attr_list = []
    for a in f.attribute.split(';'):
        key = a.split(' "')[0]
        val = a[len(key)+2:-1]
        attr_list.append((key, val))

    return GenericFeature(
        seqname=f.seqname,
        type_=f.feature,
        start=f.start,
        end=f.end,
        strand=f.strand,
        attributes=attr_list,
        frame=f.frame + 1
    )


def gff_to_generic_feature(gff_feature: GffFeature) -> GenericFeature:

    assert type(gff_feature) is GffFeature

    f = gff_feature
    items = [item for item in f.attributes.split(';') if item]

    attr_list = []
    for item in items:
        key, val = item.split('=')
        attr_list.append((key, val))

    return GenericFeature(
        seqname=f.seqid,
        type_=f.type,
        start=f.start,
        end=f.end,
        strand=f.strand,
        attributes=attr_list,
        frame=1 if f.phase == '.' else f.phase + 1
    )


def generic_to_gtf_feature(generic_feature: GenericFeature) -> GtfFeature:
    """
    Covert GenericFeature to GtfFeature (namedtuple)

    Args:
        generic_feature: GenericFeature object
    """
    assert type(generic_feature) is GenericFeature

    f = generic_feature

    # Pack attributes into a single line of str
    attr_str = ''
    for key, val in f.attributes:
        if type(val) is int or type(val) is float:
            attr_str += f"{key} {val};"
        else:  # type(val) is str -> Add quote ""
            val = val.replace(';', '<semicolon>')  # semicolon is not allowed in the GTF attribute field
            attr_str += f"{key} \"{val}\";"

    return GtfFeature(
        seqname=f.seqname,
        source='.',
        feature=f.type,
        start=f.start,
        end=f.end,
        score='.',
        strand=f.strand,
        frame=f.frame - 1,  # GTF frame is 0, 1, 2
        attribute=attr_str[:-1]  # Remove trailing ';'
    )


def generic_to_gff_feature(generic_feature: GenericFeature) -> GffFeature:
    """
    Covert GenericFeature to GffFeature (namedtuple)

    Args:
        generic_feature: GenericFeature object
    """
    assert type(generic_feature) is GenericFeature

    f = generic_feature

    # Pack attributes into a single line of str
    attr_str = ''
    for key, val in f.attributes:
        # Escape characters not allowed
        if type(val) is str:
            val = val.replace(';', '%3B').replace('=', '%3D')
        attr_str += f"{key}={val};"

    return GffFeature(
        seqid=f.seqname,
        source='.',
        type=f.type,
        start=f.start,
        end=f.end,
        score='.',
        strand=f.strand,
        phase=f.frame - 1,  # GFF frame is 0, 1, 2
        attributes=attr_str[:-1]  # Remove trailing ';'
    )
