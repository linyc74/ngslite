import shutil
from ngslite.vcf import VcfParser, VcfWriter, unpack_vcf_info
from .setup import setup_dirs, TestCase


class TestVcfParser(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_iteration(self):
        parser = VcfParser(f'{self.indir}/example.vcf')
        for variant in parser:
            self.assertTrue(type(variant) is tuple)
            self.assertEqual(12, len(variant))


class TestVcfWriter(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_write(self):
        header = '''\
##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003'''

        variant = (
            20,
            14370,
            'rs6054257',	
            'G',
            'A',
            29,
            'PASS',
            'NS=3;DP=14;AF=0.5;DB;H2',
            'GT:GQ:DP:HQ',
            '0|0:48:1:51,51',
            '1|0:48:8:51,51',
            '1/1:43:5:.,.',
        )

        output = f'{self.outdir}/output.vcf'

        with VcfWriter(file=output, header=header, mode='w') as writer:
            writer.write(variant=variant)

        line = '\t'.join(map(str, variant))
        expect = f'{header}\n{line}\n'
        with open(output) as fh:
            self.assertEqual(expect, fh.read())


class TestStaticFunctions(TestCase):

    def test_unpack_vcf_info(self):
        variant = (
            20,
            14370,
            'rs6054257',
            'G',
            'A',
            29,
            'PASS',
            'NS=3;DP=14;AF=0.5;DB;H2',
        )
        actual = unpack_vcf_info(var=variant)
        expect = {'NS': '3', 'DP': '14', 'AF': '0.5', 'DB': None, 'H2': None}
        self.assertDictEqual(expect, actual)
