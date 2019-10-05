import unittest
from ngslite.data_class import GenericFeature, FeatureArray


generic_features = list()

generic_features.append(
    GenericFeature(
        seqname='starwars',
        type_='DNA',
        start=74,
        end=1983,
        strand='+',
        attributes=None,
        tags=None,
        regions=None,
        frame=1,
        partial_start=False,
        partial_end=False
    )
)

generic_features.append(
    GenericFeature(
        seqname='starwars',
        type_='DNA',
        start=11,
        end=1977,
        strand='+',
        attributes=None,
        tags=None,
        regions=None,
        frame=1,
        partial_start=False,
        partial_end=False
    )
)

generic_features.append(
    GenericFeature(
        seqname='starwars',
        type_='DNA',
        start=10001,
        end=11977,
        strand='+',
        attributes=None,
        tags=None,
        regions=None,
        frame=1,
        partial_start=False,
        partial_end=False
    )
)

feature_array = FeatureArray(
    seqname='starwars',
    genome_size=1000000,
    features=generic_features,
    circular=False
)


class TestGenericFeature(unittest.TestCase):

    def test_(self):
        r = generic_features[1].get_attribute()


class TestFeatureArray(unittest.TestCase):

    def test_(self):
        feature_array.sort()


if __name__ == '__main__':
    unittest.main()

