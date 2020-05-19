import unittest
from ngslite.dataclass import GenericFeature
from ngslite.genbank_write import wrap_line_by_word, wrap_line_by_char, \
    generic_feature_to_genbank_text, translate_feature


class TestGenbankWrite(unittest.TestCase):

    def test_wrap_line_by_word(self):
        line = '''This is a good world with good intentions, where people love each other every day.'''
        wrapped = '''\
        This is a good world with
        good intentions, where
        people love each other
        every day.'''
        self.assertEqual(
            wrapped,
            wrap_line_by_word(
                line=line, length=33, indent=8, sep=' ', keep_sep=False)
        )

        line = '''This-is-a-good-world-with-good-intentions,-where-people-love-each-other-every-day.'''
        wrapped = '''\
        This-is-a-good-world-
        with-good-intentions,-
        where-people-love-each-
        other-every-day.'''
        self.assertEqual(
            wrapped,
            wrap_line_by_word(
                line=line, length=33, indent=8, sep='-', keep_sep=True)
        )

    def test_wrap_line_by_char(self):
        line = 'ABCDEFGHIJKLMNOPQRSEUVWXYZ'
        wrapped = '''\
        ABCDEF
        GHIJKL
        MNOPQR
        SEUVWX
        YZ'''
        self.assertEqual(
            wrapped,
            wrap_line_by_char(line=line, length=14, indent=8)
        )

    def test_generic_feature_to_genbank_text(self):
        feature = GenericFeature(
            seqname='yuchenglin',
            type_='DNA',
            start=74,
            end=1983,
            strand='-',
            attributes=[
                ('gene', 'engineering'),
                ('accession', 'R97b47402'),
                ('age', 36),
                ('distance', 42.195),
                ('description', 'Yu-Cheng is a computational biologist who tries to do something \
cool at the intersection of biology and programming. He also likes to run and play guitar.'),
                ('translation', 'GIBBERISHGIBBERISHGIBBERISHGIBBERISHGIBBERISHGIBBERISHGIBBERISH'),
            ],
            regions=[(74, 740, '-'), (1000, 1100, '-'), (1200, 1300, '-'), (1400, 1500, '-'), (1900, 1983, '-')],
            frame=1,
            partial_start=True,
            partial_end=True
        )

        text = '''\
     DNA             complement(join(<74..740,1000..1100,1200..1300,1400..1500,
                     1900..>1983))
                     /gene="engineering"
                     /accession="R97b47402"
                     /age=36
                     /distance=42.195
                     /description="Yu-Cheng is a computational biologist who
                     tries to do something cool at the intersection of biology
                     and programming. He also likes to run and play guitar."
                     /translation="GIBBERISHGIBBERISHGIBBERISHGIBBERISHGIBBERIS
                     HGIBBERISHGIBBERISH"
'''
        self.assertEqual(
            text,
            generic_feature_to_genbank_text(feature=feature)
        )

    def test_translate_feature(self):
        feature = GenericFeature(
            seqname='seqname',
            type_='CDS',
            start=1,
            end=250,
            strand='-',
            attributes=[('codon_start', 1)],
            regions=[(1, 50, '-'), (101, 150, '-'), (201, 250, '-')],
            partial_start=False,
            partial_end=False
        )

        sequence = '''\
ttcgtcaccggtatcgctttccactttgcctaccgcaccatggacaagca\
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
ggaagacgctctcaacctgctgcccactggtcaccttggcactgaggagc\
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\
aggctgccgacattgagaggcgccagagcgttactgccgagaagcactaa'''

        translation = '''\
MVLLGSNALAPLNVGSLLLSAKVTSGQQVESVFLLVHGAVGKVESDTGDE'''

        self.assertEqual(
            translation,
            translate_feature(feature=feature, sequence=sequence)
        )


if __name__ == '__main__':
    unittest.main()

