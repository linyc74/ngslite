import unittest
from ngslite.dnatools import rev_comp, translate, base_content


class TestDnaTools(unittest.TestCase):

    def test_rev_comp(self):
        self.assertEqual('AGTCagct', rev_comp('agctGACT'))
        self.assertEqual('anAN', rev_comp('NTnt'))

    def test_translate(self):
        dna = 'TAATCGTTTGCAATTCAGGGcTTGATcTACaCTggATTGCCAttcTCTCAAAGTATTaTGCAGGACGGCGTGCgCGTTCCATGtaaaCCTgTCATAACTT'
        prot = '*SFAIQGLIYTGLPFSQSIMQDGVRVPCKPVIT'
        self.assertEqual(prot, translate(dna))

    def test_base_content(self):
        dna = 'AaaCccGgTt'
        for base, content in [
            ('A', 0.3), ('a', 0.3),
            ('C', 0.3), ('c', 0.3),
            ('G', 0.2), ('g', 0.2),
            ('T', 0.2), ('t', 0.2),
        ]:
            self.assertEqual(content, base_content(seq=dna, base=base))


if __name__ == '__main__':
    unittest.main()
