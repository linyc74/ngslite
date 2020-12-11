import shutil
from ngslite.fasta import FastaParser, FastaWriter
from .setup import setup_dirs, TestCase


class TestFastaParser(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_iteration(self):
        with FastaParser(f'{self.indir}/test.fa') as parser:
            actual = [(head, seq) for head, seq in parser]
        expected = [
            ('HEADER_1', 'ACGT'),
            ('HEADER_2', 'TCGA'),
        ]
        self.assertListEqual(expected, actual)


class TestFastaWriter(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_write(self):
        header = 'HEADER'
        seq = 'M' * 81
        output = f'{self.outdir}/output.fa'

        with FastaWriter(file=output, mode='w') as writer:
            writer.write(header=header, sequence=seq)

        a = 'M' * 80
        b = 'M'
        expected = f'>{header}\n{a}\n{b}\n'
        with open(output) as fh:
            self.assertEqual(expected, fh.read())
