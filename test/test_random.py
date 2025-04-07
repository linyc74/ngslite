from ngslite.random import random_sample
from .setup import TestCase


class TestRandom(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_fa(self):

        file = f'{self.indir}/NC_000866.fa'
        fraction = 0.1
        output = f'{self.outdir}/NC_000866_{fraction}.fa'

        random_sample(
            file=file,
            fraction=fraction,
            output=output)

        i = 0
        with open(output) as fh:
            for line in fh:
                if line.startswith('>'):
                    i += 1

        self.assertEqual(27, i)
