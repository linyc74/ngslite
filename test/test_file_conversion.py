import os
import shutil
from ngslite.file_conversion import sam_to_bam, bam_to_sam, fq_to_fa, vcf_to_bcf, bcf_to_vcf
from ngslite.random import random_sample
from ngslite.count import count_reads
from .setup import setup_dirs, TestCase


class TestFileConversion(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_sam_to_bam(self):

        file = f'{self.indir}/example.sam'
        output = f'{self.outdir}/example.bam'
        correct_output = f'{self.indir}/example.bam'

        path = sam_to_bam(file=file, keep=True, output=output)

        self.assertTrue(os.path.exists(path))

        with open(output, 'rb') as fh1:
            with open(correct_output, 'rb') as fh2:
                self.assertTrue(fh1.read(), fh2.read())

    def test_bam_to_sam(self):

        file = f'{self.indir}/example.bam'
        output = f'{self.outdir}/example.sam'
        correct_output = f'{self.indir}/example.sam'

        path = bam_to_sam(file=file, keep=True, output=output)

        self.assertTrue(os.path.exists(path))

        self.assertFileEqual(file1=output, file2=correct_output)

    def test_fq_to_fa(self):

        file = f'{self.indir}/example.1.fq'
        output = f'{self.outdir}/example.1.fa'
        correct_output = f'{self.indir}/example.1.fa'

        path = fq_to_fa(file=file, keep=True, output=output)

        self.assertTrue(os.path.exists(path))

        self.assertFileEqual(file1=output, file2=correct_output)

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
