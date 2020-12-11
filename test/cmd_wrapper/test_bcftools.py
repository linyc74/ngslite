import os
import shutil
from ngslite import vcf_to_bcf, bcf_to_vcf
from ..setup import setup_dirs, TestCase


class TestBcftools(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_vcf_to_bcf(self):

        file = f'{self.indir}/example.vcf'
        output = f'{self.outdir}/example.bcf'
        correct_output = f'{self.indir}/example.bcf'

        path = vcf_to_bcf(file=file, keep=True, output=output)

        self.assertTrue(os.path.exists(path))

        with open(output, 'rb') as fh1:
            with open(correct_output, 'rb') as fh2:
                self.assertTrue(fh1.read(), fh2.read())

    def test_bcf_to_vcf(self):
        file = f'{self.indir}/example.bcf'
        output = f'{self.outdir}/example.vcf'
        correct_output = f'{self.indir}/example.vcf'

        path = bcf_to_vcf(file=file, keep=True, output=output)

        self.assertTrue(os.path.exists(path))

        self.assertFileEqual(file1=output, file2=correct_output)
