import shutil
from ngslite.fasta_gff import read_fasta_gff
from ngslite.dataclass import GenericFeature
from .setup import setup_dirs, TestCase


class TestFunctions(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_read_fasta_gff(self):
        test_fna = f'{self.indir}/test.fna'
        test_gff = f'{self.indir}/test.gff'

        chromosome = read_fasta_gff(fasta=test_fna, gff=test_gff)[0]

        f1 = GenericFeature(
            seqname='SEQNAME',
            type_='CDS',
            start=1,
            end=9,
            strand='+',
            attributes=[('key1', 'value 1'), ('key2', 'value 2')],
            tags=[],
            regions=[(1, 9, '+')],
            frame=1,
            partial_start=False,
            partial_end=False)

        f2 = GenericFeature(
            seqname='SEQNAME',
            type_='CDS',
            start=21,
            end=29,
            strand='+',
            attributes=[('key3', 'value 3'), ('key4', 'value 4')],
            tags=[],
            regions=[(21, 29, '+')],
            frame=1,
            partial_start=False,
            partial_end=False)

        self.assertEqual('SEQNAME', chromosome.seqname)
        self.assertEqual(repr(f1), repr(chromosome.features[0]))
        self.assertEqual(repr(f2), repr(chromosome.features[1]))
