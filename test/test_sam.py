import shutil
from ngslite.sam import SamParser, SamWriter, decode_flag, encode_flag
from .setup import setup_dirs, TestCase


class TestSamParser(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_iteration(self):
        with SamParser(f'{self.indir}/test.sam') as parser:

            for line in parser.header.splitlines():
                self.assertTrue(line.startswith('@'))

            for read in parser:
                self.assertTrue(type(read) is tuple)
                self.assertTrue(len(read) >= 11)


class TestSamWriter(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_write(self):
        output = f'{self.outdir}/output.sam'

        header = '@HD	VN:1.0	SO:coordinate'

        read = (
            'M05473:M05473:000000000-BYCJ5:1:1109:26192:8171',
            163,
            'dir=../18_1025_trinity_assembly/contigs/;file=Control_32mer.fa;id=TRINITY_DN87_c0_g6_i1',
            1232,
            42,
            '75M',
            '=',
            1328,
            172,
            'CCGAAGATGACAGCGACCTGATCGCCCGCGCCATTGCGGATGAGGGCTTCGTTTATCGCGGTCTTGGCGGCCTCG',
            'CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
        )

        with SamWriter(file=output, header=header, mode='w') as writer:
            writer.write(read=read)

        read_str = '\t'.join(map(str, read))
        expect = f'{header}\n{read_str}\n'

        with open(output) as fh:
            self.assertEqual(expect, fh.read())

        with SamWriter(file=output, header=header, mode='a') as writer:
            writer.write(read=read)

        expect += f'{read_str}\n'

        with open(output) as fh:
            self.assertEqual(expect, fh.read())


class TestFunctions(TestCase):

    def test_decode_flag(self):
        flag = 1365
        decoded = decode_flag(flag=flag)
        expected = {
            'read paired': True,
            'read mapped in proper pair': False,
            'read unmapped': True,
            'mate unmapped': False,
            'read reverse strand': True,
            'mate reverse strand': False,
            'first in pair': True,
            'second in pair': False,
            'not primary alignment': True,
            'read fails platform/vendor quality checks': False,
            'read is PCR or optical duplicate': True,
            'supplementary alignment': False,
        }
        self.assertDictEqual(expected, decoded)

    def test_encode_flag(self):
        flag_dict = {
            'read paired': True,
            'read mapped in proper pair': False,
            'read unmapped': True,
            'mate unmapped': False,
            'read reverse strand': True,
            'mate reverse strand': False,
            'first in pair': True,
            'second in pair': False,
            'not primary alignment': True,
            'read fails platform/vendor quality checks': False,
            'read is PCR or optical duplicate': True,
            'supplementary alignment': False,
        }
        encoded = encode_flag(flag_dict=flag_dict)
        self.assertEqual(1365, encoded)

