import unittest
import ngslite.genbank_parse as genbank_parse
from ngslite.data_class import GenericFeature


LOCUS = '''\
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

FEATURES = '''\
FEATURES             Location/Qualifiers
     source          1..168903
                     /organism="Enterobacteria phage T4"
                     /mol_type="genomic DNA"
                     /db_xref="taxon:10665"
     gene            complement(2200..2403)
                     /gene="rIIA.1"
                     /locus_tag="T4p002"
                     /db_xref="GeneID:1258774"
     CDS             join(complement(2853..2990),complement(2458..2802))
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
     gene            join(1..100,1001..2000)
                     /gene="good"
                     /note="it's a good gene"
'''

features = [''] * 4

features[0] = '''\
     source          1..168903
                     /organism="Enterobacteria phage T4"
                     /mol_type="genomic DNA"
                     /db_xref="taxon:10665"
'''

features[1] = '''\
     gene            complement(2200..2403)
                     /gene="rIIA.1"
                     /locus_tag="T4p002"
                     /db_xref="GeneID:1258774"
'''

features[2] = '''\
     CDS             join(complement(2853..2990),complement(2458..2802))
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

features[3] = '''\
     gene            join(501..1000,1..100)
                     /gene="good"
                     /note="it's a good gene"
'''

ORIGIN = '''\
ORIGIN      
        1 tgaatggttc attatctaaa aacttgttag caaccttaga tcaggttgaa gctgaggata
       61 tcgcacataa ccaactagaa ttagaacggg ctggtctgga tgatcgagcg tttggtggtg
      121 tta
'''

generic_feature = GenericFeature(
    seqname='yuchenglin',
    type_='DNA',
    start=74,
    end=1983,
    strand='+',
    attributes=[
        ('gene', 'AAA  [M=132]'),
        ('accession', 'PF00004.29'),
        ('translation', 'GIBBERISH'),
        ('age', 36),
        ('distance', 42.195)
    ]
)


class TestGenbankParse(unittest.TestCase):

    def test__get_seqname(self):
        r = genbank_parse._get_seqname(LOCUS)
        self.assertEqual(r, 'NC_000866')

    def test__is_circular(self):
        r = genbank_parse._is_circular(LOCUS)
        self.assertFalse(r)

    def test__split_feature_text(self):
        r = genbank_parse._split_feature_text(FEATURES)
        for i in range(3):
            self.assertEqual(r[i], features[i])

    def test__get_feature_type(self):
        r = genbank_parse._get_feature_type(features[1])
        self.assertEqual(r, 'gene')
        r = genbank_parse._get_feature_type(features[2])
        self.assertEqual(r, 'CDS')

    def test__get_feature_location(self):
        start, end, strand, regions = genbank_parse._get_feature_location(features[2])
        self.assertEqual(start, 2458)
        self.assertEqual(end, 2990)
        self.assertEqual(strand, '-')
        self.assertEqual(regions[0], (2458, 2802, '-'))
        self.assertEqual(regions[1], (2853, 2990, '-'))

        start, end, strand, regions = genbank_parse._get_feature_location(features[3])
        self.assertEqual(start, 1)
        self.assertEqual(end, 1000)
        self.assertEqual(strand, '+')
        self.assertEqual(regions[0], (1, 100, '+'))
        self.assertEqual(regions[1], (501, 1000, '+'))

    def test__get_feature_attributes(self):
        attr = genbank_parse._get_feature_attributes(features[2])
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

    def test__get_sequence(self):
        seq = genbank_parse._get_sequence(ORIGIN)
        self.assertEqual(seq, 'tgaatggttcattatctaaaaacttgttagcaaccttagatcaggttgaagctgaggat\
atcgcacataaccaactagaattagaacgggctggtctggatgatcgagcgtttggtggtgtta')

    def test__pack_attributes(self):
        line = genbank_parse._pack_attributes(
            feature=generic_feature,
            skip_attributes='translation'
        )
        self.assertEqual(line, 'gene "AAA  [M=132]";accession "PF00004.29";age 36;distance 42.195')


if __name__ == '__main__':
    unittest.main()

