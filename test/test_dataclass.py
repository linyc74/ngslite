import unittest
from ngslite.dataclass import GenericFeature, FeatureArray


class TestGenericFeature(unittest.TestCase):

    def setUp(self):
        self.feature = GenericFeature(
            seqname='seqname',
            type_='DNA',
            start=1,
            end=1000,
            strand='+',
            attributes=[('key_1', 'val_1'), ('key_2', 'val_2.0'), ('key_2', 'val_2.1')],
            tags=None,
            regions=[(1, 100, '+'), (200, 300, '+'), (900, 1000, '+')],
            frame=2,
            partial_start=True,
            partial_end=False)

    def test___len__(self):
        self.assertEqual(1000, len(self.feature))

        feature = GenericFeature(
            seqname='seqname',
            type_='DNA',
            start=999,
            end=10,
            strand='+',
            chromosome_size=1000)

        self.assertEqual(12, len(feature))

    def test_get_attribute(self):
        self.assertEqual(self.feature.get_attribute('key_1'), 'val_1')
        self.assertEqual(self.feature.get_attribute('key_2'), ['val_2.0', 'val_2.1'])

    def test_set_attribute(self):
        self.feature.set_attribute(key='key_1', val='new_val_1')
        self.feature.set_attribute(key='key_2', val='new_val_2')
        self.assertListEqual(
            self.feature.attributes,
            [('key_1', 'new_val_1'), ('key_2', 'new_val_2')]
        )

    def test_add_attribute(self):
        self.feature.add_attribute(key='key_3', val='val_3')
        self.assertListEqual(
            self.feature.attributes,
            [('key_1', 'val_1'), ('key_2', 'val_2.0'), ('key_2', 'val_2.1'), ('key_3', 'val_3')]
        )

    def test_remove_attribute(self):
        self.feature.remove_attribute(key='key_2')
        self.assertEqual(self.feature.attributes, [('key_1', 'val_1')])
        self.feature.remove_attribute(key='key_1')
        self.assertEqual(self.feature.attributes, [])


class TestFeatureArray(unittest.TestCase):

    def setUp(self):
        self.seqname = 'seqname'
        self.chromosome_size = 100

        self.feature_list = []

        for start, end in [
            (1, 10),
            (21, 30),
            (41, 50),
        ]:
            f = GenericFeature(
                seqname=self.seqname,
                type_='CDS',
                start=start,
                end=end,
                strand='+')
            self.feature_list.append(f)

        self.feature_array = FeatureArray(
            seqname=self.seqname,
            chromosome_size=self.chromosome_size,
            features=self.feature_list,
            circular=True)

    def test___getitem__(self):
        self.assertEqual(
            id(self.feature_array[0]),
            id(self.feature_list[0])
        )

        self.assertNotEqual(
            id(self.feature_array),
            id(self.feature_array[:])
        )

        for id1, id2, in zip(self.feature_array[:], self.feature_array):
            self.assertEqual(id1, id2)

    def test___add__(self):
        new_feature = GenericFeature(
            seqname=self.seqname,
            type_='CDS',
            start=11,
            end=15,
            strand='+')

        new_array = FeatureArray(
            seqname=self.seqname,
            chromosome_size=self.chromosome_size,
            features=new_feature,
            circular=True)

        new_array = self.feature_array + new_array

        self.assertNotEqual(id(new_array), id(self.feature_array))
        self.assertEqual(4, len(new_array))

        self.assertEqual(id(new_array[0]), id(self.feature_array[0]))
        self.assertEqual(id(new_array[2]), id(self.feature_array[1]))
        self.assertEqual(id(new_array[3]), id(self.feature_array[2]))

    def test_append_and_sort(self):
        feature = GenericFeature(
            seqname=self.seqname,
            type_='CDS',
            start=11,
            end=15,
            strand='+')

        self.feature_array.append(feature=feature)  # append includes sort

        self.assertEqual(4, len(self.feature_array))

        start_positions = [f.start for f in self.feature_array]
        self.assertListEqual(sorted(start_positions), start_positions)

    def test_subset(self):
        self.feature_array.subset(start=20, end=51)

        self.assertEqual(21, self.feature_array[0].start)
        self.assertEqual(30, self.feature_array[0].end)

        self.assertEqual(41, self.feature_array[1].start)
        self.assertEqual(50, self.feature_array[1].end)

    def test_offset(self):
        pass

    def test_crop(self):
        pass

    def test_reverse(self):
        pass

    def test_pop(self):
        pass

    def test_set_seqname(self):
        oldname = self.feature_array.seqname

        self.feature_array.set_seqname(seqname=f'{oldname}+')

        self.assertEqual(f'{oldname}+', self.feature_array.seqname)

        for feature in self.feature_array:
            self.assertEqual(f'{oldname}+', feature.seqname)
