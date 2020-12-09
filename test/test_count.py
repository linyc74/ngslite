import shutil
from ngslite.count import count_reads
from .setup import setup_dirs, TestCase


class TestCount(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_count_fq(self):
        file = f'{self.indir}/example.1.fq'
        n = count_reads(file=file)
        self.assertEqual(2546, n)

    def test_count_fq_gz(self):
        file = f'{self.indir}/example.1.fq.gz'
        n = count_reads(file=file)
        self.assertEqual(2546, n)

    def test_count_fa(self):
        file = f'{self.indir}/example.1.fa'
        n = count_reads(file=file)
        self.assertEqual(2546, n)

    def test_count_fa_gz(self):
        file = f'{self.indir}/example.1.fa.gz'
        n = count_reads(file=file)
        self.assertEqual(2546, n)

    def test_count_sam(self):
        file = f'{self.indir}/example.sam'
        n = count_reads(file=file)
        self.assertEqual(1430, n)

    def test_count_bam(self):
        file = f'{self.indir}/example.bam'
        n = count_reads(file=file)
        self.assertEqual(1430, n)
