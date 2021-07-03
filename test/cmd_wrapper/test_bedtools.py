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

    def test_bed_input(self):
        bed = f'{self.indir}/genes.bed'
        sam = f'{self.indir}/reads.sam'
        output = f'{self.outdir}/output.tsv'

        bam = sam_to_bam(file=sam, keep=True)

        bedtools_coverage(
            bed_or_gff=bed,
            bam=bam,
            output=output,
            only_counts=False,
            strand=None)

        os.remove(bam)

        self.assertFileEqual(
            file1=output,
            file2=f'{self.indir}/output_bed.tsv'
        )

    def test_gff_input(self):
        gff = f'{self.indir}/genes.gff'
        sam = f'{self.indir}/reads.sam'
        output = f'{self.outdir}/output.tsv'

        bam = sam_to_bam(file=sam, keep=True)

        bedtools_coverage(
            bed_or_gff=gff,
            bam=bam,
            output=output,
            only_counts=False,
            strand='opposite')

        os.remove(bam)

        self.assertFileEqual(
            file1=output,
            file2=f'{self.indir}/output_gff.tsv'
        )
