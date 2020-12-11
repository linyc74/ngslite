import shutil
from ngslite.fastq import FastqParser, FastqWriter
from .setup import setup_dirs, TestCase


class TestFastqParser(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_iteration(self):
        with FastqParser(f'{self.indir}/test.1.fq') as parser:
            actual = [item for item in parser]

        expected = [
            (
                '@M03991:M03991:000000000-BV6FD:1:1104:10085:22424 1:Y:0:GCCAAT',
                 'ATAATAATGTAGCTGGTGTTATATTTACGGGAACATGTTTAAAATCAAGCAAGCATTCTTTACTGTATGCTCAAGA',
                 '+',
                 'CC@CCE<@FC6,CCF<,,;B@@,<;CFCC;@;@CF,C6<,C<<EF<@FA9,<@F8,<,,C,CFGG<@6FC<FEFCF'
            ),
            (
                '@M03991:M03991:000000000-BV6FD:1:1104:10331:19090 1:N:0:GCCAAT',
                 'TACTTTATATTTGCTTCCTGCTATCTCTTGCATTTCTTTCTTGACCTCATCAAGGTTACTGCCATCCTTAAGTCTT',
                 '+',
                 'CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGFGGGGGGGGGF'
            ),
        ]

        self.assertListEqual(expected, actual)


class TestFastqWriter(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_write(self):
        output = f'{self.outdir}/output.fq'

        reads = [
            (
                '@M03991:M03991:000000000-BV6FD:1:1104:10085:22424 1:Y:0:GCCAAT',
                 'ATAATAATGTAGCTGGTGTTATATTTACGGGAACATGTTTAAAATCAAGCAAGCATTCTTTACTGTATGCTCAAGA',
                 '+',
                 'CC@CCE<@FC6,CCF<,,;B@@,<;CFCC;@;@CF,C6<,C<<EF<@FA9,<@F8,<,,C,CFGG<@6FC<FEFCF'
            ),
            (
                '@M03991:M03991:000000000-BV6FD:1:1104:10331:19090 1:N:0:GCCAAT',
                 'TACTTTATATTTGCTTCCTGCTATCTCTTGCATTTCTTTCTTGACCTCATCAAGGTTACTGCCATCCTTAAGTCTT',
                 '+',
                 'CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGFGGGGGGGGGF'
            ),
        ]

        with FastqWriter(file=output, mode='w') as writer:
            for read in reads:
                writer.write(read=read)

        expected = '''\
@M03991:M03991:000000000-BV6FD:1:1104:10085:22424 1:Y:0:GCCAAT
ATAATAATGTAGCTGGTGTTATATTTACGGGAACATGTTTAAAATCAAGCAAGCATTCTTTACTGTATGCTCAAGA
+
CC@CCE<@FC6,CCF<,,;B@@,<;CFCC;@;@CF,C6<,C<<EF<@FA9,<@F8,<,,C,CFGG<@6FC<FEFCF
@M03991:M03991:000000000-BV6FD:1:1104:10331:19090 1:N:0:GCCAAT
TACTTTATATTTGCTTCCTGCTATCTCTTGCATTTCTTTCTTGACCTCATCAAGGTTACTGCCATCCTTAAGTCTT
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGFGGGGGGGGGF
'''
        with open(output) as fh:
            self.assertEqual(expected, fh.read())
