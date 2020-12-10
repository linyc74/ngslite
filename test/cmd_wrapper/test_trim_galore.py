import shutil
from ngslite.cmd_wrapper.trim_galore import trim_galore
from test.setup import setup_dirs, TestCase


class TestTrimGalore(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_trim_galore(self):

        fq1 = f'{self.indir}/sample.1.fq'
        fq2 = f'{self.indir}/sample.2.fq'
        log = f'{self.outdir}/trim_galore.log'

        trim_galore(fq1=fq1, fq2=fq2, log=log)

        msg = f'''\
Multicore support not enabled. Proceeding with single-core trimming.
Path to Cutadapt set as: 'cutadapt' (default)
Cutadapt seems to be working fine (tested command 'cutadapt --version')'''

        with open(log) as fh:
            self.assertTrue(msg in fh.read())
