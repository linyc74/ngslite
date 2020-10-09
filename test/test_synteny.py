import shutil
from ngslite.synteny import synteny
from .setup import setup_dirs, TestCase


class TestSynteny(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

        self.keywords = ['DNA polymerase']
        self.flank = 1000

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_synteny(self):
        seqname_order = synteny(
            genbank=f'{self.indir}/loci.gb',
            output=f'{self.outdir}/synteny.gb')

        expected = [
            'NC_007024:37356-41278',
            'NC_007024:27507-30313',
            'NC_007024:29897-32379']

        self.assertListEqual(expected, seqname_order)

        self.assertFileEqual(
            f'{self.indir}/synteny.gb',
            f'{self.outdir}/synteny.gb')
