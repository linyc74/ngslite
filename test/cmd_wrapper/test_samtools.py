import os
import shutil
from ngslite import sam_to_bam, bam_to_sam
from ..setup import setup_dirs, TestCase


class TestSamtools(TestCase):

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
