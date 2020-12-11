import shutil
from ngslite.gff import GffParser, GffWriter, GffFeature
from .setup import setup_dirs, TestCase


class TestGffParser(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_iteration(self):

        first = GffFeature(
            seqid='SEQF1003|ACFE01000008.1',
            source='Prodigal:2.6',
            type='CDS',
            start=969,
            end=1190,
            score='.',
            strand='+',
            phase=0,
            attributes='ID=SEQF1003_00001;inference=ab initio prediction:Prodigal:2.6;locus_tag=SEQF1003_00001;product=hypothetical protein'
        )

        last = GffFeature(
            seqid='SEQF1003|ACFE01000008.1',
            source='Prodigal:2.6',
            type='CDS',
            start=2117,
            end=2320,
            score='.',
            strand='-',
            phase=0,
            attributes='ID=SEQF1003_00004;Name=cspC;gene=cspC;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P39158;locus_tag=SEQF1003_00004;product=Cold shock protein CspC'
        )

        with GffParser(f'{self.indir}/test.gff') as parser:
            header = parser.header
            features = [f for f in parser]

        self.assertEqual('##sequence-region SEQF1003_ACFE01000008.1 1 3040\n##---', header)
        self.assertEqual(repr(first), repr(features[0]))
        self.assertEqual(repr(last), repr(features[-1]))


class TestGffWriter(TestCase):

    def setUp(self):
        self.indir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_write(self):
        output = f'{self.outdir}/output.gff'

        header = '##sequence-region SEQF1003_ACFE01000008.1 1 3040'

        first = GffFeature(
            seqid='SEQF1003|ACFE01000008.1',
            source='Prodigal:2.6',
            type='CDS',
            start=969,
            end=1190,
            score='.',
            strand='+',
            phase=0,
            attributes='ID=SEQF1003_00001;inference=ab initio prediction:Prodigal:2.6;locus_tag=SEQF1003_00001;product=hypothetical protein'
        )

        second = GffFeature(
            seqid='SEQF1003|ACFE01000008.1',
            source='Prodigal:2.6',
            type='CDS',
            start=2117,
            end=2320,
            score='.',
            strand='-',
            phase=0,
            attributes='ID=SEQF1003_00004;Name=cspC;gene=cspC;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P39158;locus_tag=SEQF1003_00004;product=Cold shock protein CspC'
        )

        with GffWriter(file=output, header=header, mode='w') as writer:
            writer.write(feature=first)

        expect = f'''\
##gff-version 3
{header}
SEQF1003|ACFE01000008.1	Prodigal:2.6	CDS	969	1190	.	+	0	ID=SEQF1003_00001;inference=ab initio prediction:Prodigal:2.6;locus_tag=SEQF1003_00001;product=hypothetical protein
'''
        with open(output) as fh:
            self.assertEqual(expect, fh.read())

        with GffWriter(file=output, header=header, mode='a') as writer:
            writer.write(feature=second)

        expect += 'SEQF1003|ACFE01000008.1	Prodigal:2.6	CDS	2117	2320	.	-	0	ID=SEQF1003_00004;Name=cspC;gene=cspC;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P39158;locus_tag=SEQF1003_00004;product=Cold shock protein CspC\n'

        with open(output) as fh:
            self.assertEqual(expect, fh.read())

