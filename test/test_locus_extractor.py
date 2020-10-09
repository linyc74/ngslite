import os
import shutil
from ngslite.locus_extractor import LocusExtractor
from .setup import setup_dirs, TestCase


class TestLocusExtractor(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

        self.keywords = ['DNA polymerase']
        self.flank = 1000

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def remove_date(self, gbk: str):
        temp_gbk = f'{self.workdir}/temp.gbk'

        with open(gbk) as reader:
            with open(temp_gbk, 'w') as writer:
                for line in reader:
                    if line.startswith('LOCUS'):
                        line = line.split('linear')[0] + 'linear\n'
                    writer.write(line)

        os.remove(gbk)
        os.rename(temp_gbk, gbk)

    def test_use_gbk(self):
        gbk = f'{self.indir}/xp15.gb'
        output = f'{self.outdir}/loci.gb'

        LocusExtractor().use_gbk(
            gbk=gbk,
            keywords=self.keywords,
            flank=self.flank,
            output=output)

        self.remove_date(gbk=output)

        self.assertFileEqual(f'{self.indir}/loci.gb', output)

    def test_use_fasta_gtf(self):
        fasta = f'{self.indir}/xp15.fna'
        gtf = f'{self.indir}/xp15.gtf'
        fasta_out = f'{self.outdir}/loci.fna'
        gtf_out = f'{self.outdir}/loci.gtf'

        LocusExtractor().use_fasta_gtf(
            fasta=fasta,
            gtf=gtf,
            keywords=self.keywords,
            flank=self.flank,
            fasta_out=fasta_out,
            gtf_out=gtf_out)

        self.assertFileEqual(f'{self.indir}/loci.fna', fasta_out)
        self.assertFileEqual(f'{self.indir}/loci.gtf', gtf_out)

    def test_use_fasta_gff(self):
        fasta = f'{self.indir}/xp15.fna'
        gff = f'{self.indir}/xp15.gff'
        fasta_out = f'{self.outdir}/loci.fna'
        gff_out = f'{self.outdir}/loci.gff'

        LocusExtractor().use_fasta_gff(
            fasta=fasta,
            gff=gff,
            keywords=self.keywords,
            flank=self.flank,
            fasta_out=fasta_out,
            gff_out=gff_out)

        self.assertFileEqual(f'{self.indir}/loci.fna', fasta_out)
        self.assertFileEqual(f'{self.indir}/loci.gff', gff_out)
