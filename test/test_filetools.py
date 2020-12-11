import shutil
from ngslite.filetools import get_files, get_dirs, get_temp_path, zip_broadcast
from .setup import setup_dirs, TestCase


class TestFiletools(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_get_files(self):

        files = get_files(
            source=f'{self.indir}/test_get_files',
            startswith='',
            endswith='',
            isfullpath=False)

        expected = ['1.txt', '2.txt']

        self.assertListEqual(expected, files)

    def test_get_dirs(self):

        files = get_dirs(
            source=f'{self.indir}/test_get_dirs',
            startswith='',
            endswith='',
            isfullpath=False)

        expected = ['1', '2']

        self.assertListEqual(expected, files)

    def test_get_temp_path(self):

        self.assertEqual('temp001', get_temp_path())

        path = get_temp_path(
            prefix=f'{self.indir}/existing',
            suffix='.txt'
        )
        expected = f'{self.indir}/existing002.txt'
        self.assertEqual(expected, path)

    def test_zip_broadcast(self):

        list_1 = [1, 2, 3, 4, 0]
        list_2 = [5, 6]

        zipped = zip_broadcast(list_1, list_2)

        expected = [(1, 5), (2, 6), (3, 5), (4, 6), (0, 5)]

        self.assertListEqual(expected, zipped)
