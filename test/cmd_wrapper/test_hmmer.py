import os
import shutil
from ngslite.cmd_wrapper.hmmer import hmmsearch, parse_hmmsearch_result
from ngslite import gzip
from test.setup import setup_dirs, TestCase


class TestHmmer(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_hmmsearch_and_parse_result(self):

        hmm = f'{self.indir}/Thymidylat_synt.hmm.gz'
        database = f'{self.indir}/xp15.fna.gz'
        output_txt = f'{self.outdir}/output.txt'
        output_gtf = f'{self.outdir}/output.gtf'

        hmmsearch(
            hmm=hmm,
            database=database,
            output=output_txt,
            cpu=6)

        gzip(database, keep=True)

        parse_hmmsearch_result(
            file=output_txt,
            output=output_gtf,
            database=database[:-3])

        os.remove(database[:-3])

        # Do not test output_txt because 'CPU time' is included in the result, which changes everytime

        self.assertFileEqual(f'{self.indir}/output.gtf', output_gtf)
