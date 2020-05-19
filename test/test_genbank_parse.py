import shutil
import unittest
from ngslite.fasta import FastaParser
from ngslite.dataclass import GenericFeature
from ngslite.genbank_parse import get_seqname, get_circular, \
    split_features_text, get_feature_type, get_feature_location, \
    get_feature_attributes, get_sequence, is_valid_first_line_of_feature, \
    construct_chromosome, genbank_to_fasta, genbank_to_gff, genbank_to_gtf, \
    read_genbank, GenbankTextParser
from .setup_dirs import setup_dirs


locus_text = '''\
LOCUS       NC_000866      168903 bp    DNA     linear   PHG 15-FEB-2017
DEFINITION  Enterobacteria phage T4, complete genome.
ACCESSION   NC_000866
KEYWORDS    .
SOURCE      Enterobacteria phage T4 (T4)
  ORGANISM  Enterobacteria phage T4
            Viruses; dsDNA viruses, no RNA stage; Caudovirales; Myoviridae;
            T4-like viruses.
REFERENCE   1  (bases 1 to 168903)
  AUTHORS   Miller,E.S., Kutter,E., Mosig,G., Arisaka,F., Kunisawa,T. and
            Ruger,W.
  TITLE     Bacteriophage T4 Genome
  JOURNAL   Microbiol. Mol. Biol. Rev. 67 (1), 86-156 (2003)
  PUBMED    12626685
REFERENCE   2  (bases 80196 to 80618)
  AUTHORS   Shcherbakov,V., Granovsky,I., Plugina,L., Shcherbakova,T.,
            Sizova,S., Pyatkov,K., Shlyapnikov,M. and Shubina,O.
  TITLE     Focused genetic recombination of bacteriophage t4 initiated by
            double-strand breaks
  JOURNAL   Genetics 162 (2), 543-556 (2002) MEDLINE 22286362
  PUBMED    12399370
'''


features_text = '''\
FEATURES             Location/Qualifiers
     source          1..168903
                     /organism="Enterobacteria phage T4"
                     /mol_type="genomic DNA"
                     /db_xref="taxon:10665"
     mRNA            join(484000..485000,486000..487000,488000..489000)
                     /gene=60
                     /locus_tag="T4p003"
                     /note="topoisomerase II in T4 and related phages an
                     insertion splits the large Topo II subunit in two parts
                     (gp60 and gp39); a 50-bp segment is inserted in T4 gene 60
                     that forms an mRNA secondary structure that is
                     translationally bypassed"
                     /product="gp60 topoisomerase II, large subunit,
                     C-terminalregion"
                     /protein_id="NP_049618.1"
                     /db_xref="GI:9632675"
                     /db_xref="GeneID:1258779"
     CDS             complement(join(<484000..485000,486000..487000,
                     488000..>489000))
                     /gene=60
                     /locus_tag="T4p003"
                     /note="topoisomerase II in T4 and related phages an
                     insertion splits the large Topo II subunit in two parts
                     (gp60 and gp39); a 50-bp segment is inserted in T4 gene 60
                     that forms an mRNA secondary structure that is
                     translationally bypassed"
                     /codon_start=1
                     /transl_table=11
                     /product="gp60 topoisomerase II, large subunit,
                     C-terminalregion"
                     /protein_id="NP_049618.1"
                     /db_xref="GI:9632675"
                     /db_xref="GeneID:1258779"
                     /translation="MKFVKIDSSSVDMKKYKLQNNVRRSIKSSSMNYANVAIMTDADH
                     DGLGSIYPSLLGFFSNWPELFEQGRIRFVKTPVIIAQVGKKQEWFYTVAEYESAKDAL
                     PKHSIRYIKGLGSLEKSEYREMIQNPVYDVVKLPENWKELFEMLMGDNADLRKEWMSQ
                     "
'''


feature_text_1 = '''\
     mRNA            join(484000..485000,486000..487000,488000..489000)
                     /gene=60
                     /locus_tag="T4p003"
                     /note="topoisomerase II in T4 and related phages an
                     insertion splits the large Topo II subunit in two parts
                     (gp60 and gp39); a 50-bp segment is inserted in T4 gene 60
                     that forms an mRNA secondary structure that is
                     translationally bypassed"
                     /product="gp60 topoisomerase II, large subunit,
                     C-terminalregion"
                     /protein_id="NP_049618.1"
                     /db_xref="GI:9632675"
                     /db_xref="GeneID:1258779"
'''


feature_text_2 = '''\
     CDS             complement(join(<484000..485000,486000..487000,
                     488000..>489000))
                     /gene=60
                     /locus_tag="T4p003"
                     /note="topoisomerase II in T4 and related phages an
                     insertion splits the large Topo II subunit in two parts
                     (gp60 and gp39); a 50-bp segment is inserted in T4 gene 60
                     that forms an mRNA secondary structure that is
                     translationally bypassed"
                     /codon_start=1
                     /transl_table=11
                     /product="gp60 topoisomerase II, large subunit,
                     C-terminalregion"
                     /protein_id="NP_049618.1"
                     /db_xref="GI:9632675"
                     /db_xref="GeneID:1258779"
                     /translation="MKFVKIDSSSVDMKKYKLQNNVRRSIKSSSMNYANVAIMTDADH
                     DGLGSIYPSLLGFFSNWPELFEQGRIRFVKTPVIIAQVGKKQEWFYTVAEYESAKDAL
                     PKHSIRYIKGLGSLEKSEYREMIQNPVYDVVKLPENWKELFEMLMGDNADLRKEWMSQ
                     "
'''


origin_text = '''\
ORIGIN      
        1 tgaatggttc attatctaaa aacttgttag caaccttaga tcaggttgaa gctgaggata
       61 tcgcacataa ccaactagaa ttagaacggg ctggtctgga tgatcgagcg tttggtggtg
      121 tta
'''


class TestGenbankParse(unittest.TestCase):

    def setUp(self):
        self.testdir, self.datadir, self.workdir, self.outdir = setup_dirs(__file__)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        shutil.rmtree(self.outdir)

    def test_get_seqname(self):
        self.assertEqual('NC_000866', get_seqname(locus_text))

    def test_get_circular(self):
        self.assertFalse(get_circular(locus_text))

    def test_split_features_text(self):
        list_ = split_features_text(features_text)
        self.assertEqual(feature_text_1, list_[1])
        self.assertEqual(feature_text_2, list_[2])

    def test_get_feature_type(self):
        self.assertEqual('mRNA', get_feature_type(feature_text_1))
        self.assertEqual('CDS', get_feature_type(feature_text_2))

    def test_get_feature_location(self):
        start, end, strand, regions, partial_start, partial_end = get_feature_location(feature_text_1)
        self.assertEqual(484000, start)
        self.assertEqual(489000, end)
        self.assertEqual('+', strand)
        self.assertEqual((484000, 485000, '+'), regions[0])
        self.assertEqual((486000, 487000, '+'), regions[1])
        self.assertEqual((488000, 489000, '+'), regions[2])
        self.assertFalse(partial_start)
        self.assertFalse(partial_end)

        start, end, strand, regions, partial_start, partial_end = get_feature_location(feature_text_2)
        self.assertEqual(484000, start)
        self.assertEqual(489000, end)
        self.assertEqual('-', strand)
        self.assertEqual((484000, 485000, '-'), regions[0])
        self.assertEqual((486000, 487000, '-'), regions[1])
        self.assertEqual((488000, 489000, '-'), regions[2])
        self.assertTrue(partial_start)
        self.assertTrue(partial_end)

    def test_get_feature_attributes(self):
        attr = get_feature_attributes(feature_text_2)
        self.assertEqual(attr[0], ('gene', 60))
        self.assertEqual(attr[1], ('locus_tag', 'T4p003'))
        self.assertEqual(attr[2], ('note', 'topoisomerase II in T4 and related phages an \
insertion splits the large Topo II subunit in two parts \
(gp60 and gp39); a 50-bp segment is inserted in T4 gene 60 \
that forms an mRNA secondary structure that is \
translationally bypassed'))
        self.assertEqual(attr[3], ('codon_start', 1))
        self.assertEqual(attr[4], ('transl_table', 11))
        self.assertEqual(attr[5], ('product', 'gp60 topoisomerase II, large subunit, C-terminalregion'))
        self.assertEqual(attr[6], ('protein_id', 'NP_049618.1'))
        self.assertEqual(attr[7], ('db_xref', 'GI:9632675'))
        self.assertEqual(attr[8], ('db_xref', 'GeneID:1258779'))
        self.assertEqual(attr[9], ('translation', 'MKFVKIDSSSVDMKKYKLQNNVRRSIKSSSMNYANVAIM\
TDADHDGLGSIYPSLLGFFSNWPELFEQGRIRFVKTPVIIAQVGKKQEWFYTVAEYESAKDALPKHSIRYIKGLGSLEKSEYREMIQNPV\
YDVVKLPENWKELFEMLMGDNADLRKEWMSQ'))

    def test_get_sequence(self):
        self.assertEqual(
            'tgaatggttcattatctaaaaacttgttagcaaccttagatcaggttgaagctgagga\
tatcgcacataaccaactagaattagaacgggctggtctggatgatcgagcgtttggtggtgtta',
            get_sequence(origin_text)
        )

    def test_is_valid_first_line_of_feature(self):
        valid_lines = [
            '     mRNA            join(1..10,100..110)'
            '     CDS             complement(1..10)'
            '     CDS             <1..10'
        ]

        invalid_lines = [
            '                     ',
            '     mRNA CDS        1..10'
        ]

        for line in valid_lines:
            self.assertTrue(is_valid_first_line_of_feature(line))

        for line in invalid_lines:
            self.assertFalse(is_valid_first_line_of_feature(line))

    def test_construct_chromosome(self):
        expected_feature = GenericFeature(
            seqname='NC_000866',
            type_='CDS',
            start=2458,
            end=2990,
            strand='-',
            attributes=[
                ('gene', '60'),
                ('locus_tag', 'T4p003'),
                ('note', 'topoisomerase II in T4 and related phages an insertion \
splits the large Topo II subunit in two parts (gp60 and gp39); a 50-bp segment is \
inserted in T4 gene 60 that forms an mRNA secondary structure that is translationally bypassed'),
                ('codon_start', 1),
                ('transl_table', 11),
                ('product', 'gp60 topoisomerase II, large subunit, C-terminal region'),
                ('protein_id', 'NP_049618.1'),
                ('db_xref', 'GeneID:1258779'),
                ('translation', 'MKFVKIDSSSVDMKKYKLQNNVRRSIKSSSMNYANVAIMTDADHD\
GLGSIYPSLLGFFSNWPELFEQGRIRFVKTPVIIAQVGKKQEWFYTVAEYESAKDALPKHSIRYIKGLGSLEKSEYRE\
MIQNPVYDVVKLPENWKELFEMLMGDNADLRKEWMSQ')
            ],
            regions=[(2458, 2802, '-'), (2853, 2990, '-')],
            frame=1,
            partial_start=False,
            partial_end=False
        )

        with GenbankTextParser(file=f'{self.datadir}/NC_000866.gbk') as parser:
            genbank_text = parser.next()
            chromosome = construct_chromosome(genbank_text=genbank_text)

        for feature in chromosome.features:
            locus_tag = feature.get_attribute('locus_tag')
            if feature.type == 'CDS' and locus_tag == 'T4p003':
                self.assertEqual(repr(expected_feature), repr(feature))
                break

    def test_genbank_to_fasta(self):
        fname = 'NC_000866'
        genbank_to_fasta(
            file=f'{self.datadir}/{fname}.gbk',
            output=f'{self.outdir}/{fname}.fna'
        )

        parser1 = FastaParser(f'{self.outdir}/{fname}.fna')
        parser2 = FastaParser(f'{self.datadir}/{fname}.fna')

        for item1, item2 in zip(parser1, parser2):
            _, seq1 = item1
            _, seq2 = item2
            self.assertEqual(seq1.lower(), seq2.lower())

        parser1.close()
        parser2.close()

    def test_genbank_to_gtf(self):
        fname = 'NC_000866'
        genbank_to_gtf(
            file=f'{self.datadir}/{fname}.gbk',
            output=f'{self.outdir}/{fname}.gtf',
            skip_types=''
        )

    def test_genbank_to_gff(self):
        fname = 'NC_000866'
        genbank_to_gff(
            file=f'{self.datadir}/{fname}.gbk',
            output=f'{self.outdir}/{fname}.gff3',
            skip_types=''
        )

    def test_read_genbank(self):
        target_feature = GenericFeature(
            seqname='NC_000866',
            type_='CDS',
            start=2458,
            end=2990,
            strand='-',
            attributes=[
                ('gene', '60'),
                ('locus_tag', 'T4p003'),
                ('note', 'topoisomerase II in T4 and related phages an insertion \
splits the large Topo II subunit in two parts (gp60 and gp39); a 50-bp segment is \
inserted in T4 gene 60 that forms an mRNA secondary structure that is translationally bypassed'),
                ('codon_start', 1),
                ('transl_table', 11),
                ('product', 'gp60 topoisomerase II, large subunit, C-terminal region'),
                ('protein_id', 'NP_049618.1'),
                ('db_xref', 'GeneID:1258779'),
                ('translation', 'MKFVKIDSSSVDMKKYKLQNNVRRSIKSSSMNYANVAIMTDADHD\
GLGSIYPSLLGFFSNWPELFEQGRIRFVKTPVIIAQVGKKQEWFYTVAEYESAKDALPKHSIRYIKGLGSLEKSEYRE\
MIQNPVYDVVKLPENWKELFEMLMGDNADLRKEWMSQ')
            ],
            regions=[(2458, 2802, '-'), (2853, 2990, '-')],
            frame=1,
            partial_start=False,
            partial_end=False
        )

        fname = 'NC_000866'
        chromosomes = read_genbank(file=f'{self.datadir}/{fname}.gbk')
        for chromosome in chromosomes:
            for feature in chromosome.features:
                locus_tag = feature.get_attribute('locus_tag')
                if feature.type == 'CDS' and locus_tag == 'T4p003':
                    self.assertEqual(repr(feature), repr(target_feature))
                    break


if __name__ == '__main__':
    unittest.main()
