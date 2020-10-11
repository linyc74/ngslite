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

        infile = f'{self.indir}/example.sam'
        outfile = f'{self.outdir}/example.bam'
        checkfile = f'{self.indir}/example.bam'

        path = sam_to_bam(file=infile, keep=True, output=outfile)

        self.assertTrue(os.path.exists(path))

        with open(outfile, 'rb') as fh1:
            with open(checkfile, 'rb') as fh2:
                self.assertTrue(fh1.read(), fh2.read())

    def test_bam_to_sam(self):
        infile = f'{self.indir}/example.bam'
        outfile = f'{self.outdir}/example.sam'
        checkfile = f'{self.indir}/example.sam'

        path = bam_to_sam(file=infile, keep=True, output=outfile)

        self.assertTrue(os.path.exists(path))

        self.assertFileEqual(file1=outfile, file2=checkfile)
