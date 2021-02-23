import os
import shutil
from ngslite import bedtools_coverage, sam_to_bam
from ..setup import setup_dirs, TestCase


class TestBedtools(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_bedtools_coverage(self):
        bed = f'{self.indir}/genes.bed'
        sam = f'{self.indir}/reads.sam'
        output = f'{self.outdir}/output.tsv'

        bam = sam_to_bam(file=sam, keep=True)

        bedtools_coverage(
            bed=bed,
            bam=bam,
            output=output)

        os.remove(bam)

        self.assertFileEqual(
            file1=output,
            file2=f'{self.indir}/output.tsv'
        )
