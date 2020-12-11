import os
import shutil
from ngslite import fq_to_fa
from test.setup import setup_dirs, TestCase


class TestSeqtk(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_fq_to_fa(self):

        file = f'{self.indir}/example.1.fq'
        output = f'{self.outdir}/example.1.fa'
        correct_output = f'{self.indir}/example.1.fa'

        path = fq_to_fa(file=file, keep=True, output=output)

        self.assertTrue(os.path.exists(path))

        self.assertFileEqual(file1=output, file2=correct_output)
