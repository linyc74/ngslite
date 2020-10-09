import shutil
from ngslite.filetools import get_files, get_dirs
from .setup import setup_dirs, TestCase


class TestCmdToolkit(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_get_files(self):

        files = get_files(
            source=f'{self.indir}',
            startswith='',
            endswith='',
            isfullpath=False)

        expected = ['1.txt', '2.txt']

        self.assertListEqual(expected, files)

    def test_get_dirs(self):

        files = get_dirs(
            source=f'{self.indir}',
            startswith='',
            endswith='',
            isfullpath=False)

        expected = ['1', '2']

        self.assertListEqual(expected, files)
