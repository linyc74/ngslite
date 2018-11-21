import unittest


from ngslite.dnatools import *


class TestDnaTools(unittest.TestCase):
    def test_rev_comp(self):
        self.assertEqual('AGTCagct', rev_comp('agctGACT'))

    def test_translate(self):
        dna = 'TAATCGTTTGCAATTCAGGGcTTGATcTACaCTggATTGCCAttcTCTCAAAGTATTaTGCAGGACGGCGTGCgCGTTCCATGtaaaCCTgTCATAACTT'
        prot = '*SFAIQGLIYTGLPFSQSIMQDGVRVPCKPVIT'
        self.assertEqual(prot, translate(dna))

    def test_base_content(self):
        dna = 'AaaCccGgTt'
        self.assertEqual(0.3, base_content(seq=dna, base='A'))
        self.assertEqual(0.3, base_content(seq=dna, base='C'))
        self.assertEqual(0.2, base_content(seq=dna, base='G'))
        self.assertEqual(0.2, base_content(seq=dna, base='T'))
        self.assertEqual(0.3, base_content(seq=dna, base='a'))
        self.assertEqual(0.3, base_content(seq=dna, base='c'))
        self.assertEqual(0.2, base_content(seq=dna, base='g'))
        self.assertEqual(0.2, base_content(seq=dna, base='t'))

