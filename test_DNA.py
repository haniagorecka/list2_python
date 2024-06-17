import unittest
from unittest import TestCase

from DNA import ProteinSequence, RNASequence, DNASequence


class TestProteinSequence(TestCase):


    def setUp(self):
        self.tab = ['Lys', 'Ile', 'Asp', 'Ala']
        self.prot = ProteinSequence(self.tab)
    def test_create(self):
        self.assertEqual(self.prot.data, self.tab, "Blad tworzenia obiektu")
        self.assertEqual(self.prot.length, len(self.tab), "Blad tworzenia obiektu")
        self.assertEqual(self.prot.identifier, "ProteinSeq4 LIAA", "Blad tworzenia obiektu")
        with self.assertRaises(TypeError):
            ProteinSequence(3)
        with self.assertRaises(ValueError):
            ProteinSequence(["X", "Y"])
        with self.assertRaises(ValueError):
            ProteinSequence("XYZ")

    def test_mutate(self):
        self.prot.mutate(1, 'Ala')
        self.assertEqual(self.prot.data, ['Lys', 'Ala', 'Asp', 'Ala'], "Blad funkcji mutującej")
        with self.assertRaises(IndexError):
            self.prot.mutate(7, 'Ala')
        with self.assertRaises(IndexError):
            self.prot.mutate(4, 'Ala')
        with self.assertRaises(ValueError):
            self.prot.mutate(1, 'X')
        with self.assertRaises(ValueError):
            self.prot.mutate(1, '')
        with self.assertRaises(ValueError):
            prot2 = ProteinSequence("")
            prot2.mutate(1, 'Ala')

    def test_find_motif(self):
        self.assertEqual(self.prot.find_motif(['Asp']), 2, "Blad funkcji find_motif")
        self.assertEqual(self.prot.find_motif(['Asp', 'Ala']), 2, "Blad funkcji find_motif")
        self.assertEqual(self.prot.find_motif(['Lys', 'Ile', 'Asp', 'Ala']), 0, "Blad funkcji find_motif")
        with self.assertRaises(IndexError):
            self.prot.find_motif(['Lys', 'Ile', 'Asp', 'Ala', 'Ala'])
        with self.assertRaises(TypeError):
            self.prot.find_motif(2)
        with self.assertRaises(ValueError):
            self.prot.find_motif(None)
        with self.assertRaises(ValueError):
            self.prot.find_motif(['Leu'])

    def test_str(self):
        string = f'>{self.prot.identifier}\n'
        for el in self.prot.data:
            string += el
        print("S P R A W DZ A M" + string)
        string1 = str(self.prot)
        self.assertEqual(string1, string, "Błąd funkcji str()")
        with self.assertRaises(ValueError):
            self.prot.data = None
            print(self.prot)


class TestRNASequence(TestCase):

    def setUp(self):
        self.tab = "AUAgaUagga"
        self.tab1 = "aagaggauuuaaa"
        self.rna = RNASequence(self.tab)
        self.rna1 = RNASequence(self.tab1)
    def test_create(self):

        self.assertEqual(self.rna.data, self.tab.upper(), "Blad tworzenia obiektu")
        self.assertEqual(self.rna.length, len(self.tab), "Blad tworzenia obiektu")
        self.assertEqual(self.rna.identifier, f'RNA{len(self.tab)} 5A0C2U3G', "Blad tworzenia obiektu")
        self.assertEqual(self.rna1.data, self.tab1.upper(), "Blad tworzenia obiektu")
        self.assertEqual(self.rna1.length, len(self.tab1), "Blad tworzenia obiektu")
        self.assertEqual(self.rna1.identifier, f'RNA{len(self.tab1)} 7A0C3U3G', "Blad tworzenia obiektu")
        with self.assertRaises(TypeError):
            RNASequence(3)
        with self.assertRaises(ValueError):
            RNASequence("XXX")

    def test_translate(self):
        prot2 = self.rna.translate()
        self.assertEqual(prot2.data, ['Ile', 'Asp', 'Arg'], "Blad funkcji translate()1")
        self.assertEqual(self.rna1.translate().data, ['Lys', 'Arg', 'Ile'], "Blad funkcji translate()2")
        self.assertEqual(RNASequence("").translate().data, [], "Blad funkcji translate()3")
        self.assertEqual(RNASequence("AG").translate().data, [], "Blad funkcji translate()4")
        with self.assertRaises(ValueError):
            RNASequence("XYZ").translate()

    def test_mutate(self):
        self.rna.mutate(1, 'G')
        self.assertEqual(self.rna.data, "AGAgaUagga".upper(), "Blad funkcji mutującej")
        with self.assertRaises(IndexError):
            self.rna.mutate(18, 'G')
        with self.assertRaises(IndexError):
            self.rna.mutate(10, 'G')
        with self.assertRaises(ValueError):
            self.rna.mutate(1, 'X')
        with self.assertRaises(ValueError):
            self.rna.mutate(1, '')
        with self.assertRaises(ValueError):
            rna1 = RNASequence("")
            rna1.mutate(1, 'G')

    def test_find_motif(self):
        self.assertEqual(self.rna.find_motif("A"), 0, "Blad funkcji find_motif")
        self.assertEqual(self.rna.find_motif("AUA"), 0, "Blad funkcji find_motif")
        self.assertEqual(self.rna.find_motif("GG"), 7, "Blad funkcji find_motif")
        self.assertEqual(self.rna.find_motif("gg"), 7, "Blad funkcji find_motif")
        with self.assertRaises(IndexError):
            self.rna.find_motif("aaaaaaaaaaaaaaaaa")
        with self.assertRaises(TypeError):
            self.rna.find_motif(2)
        with self.assertRaises(ValueError):
            self.rna.find_motif(None)
        with self.assertRaises(ValueError):
            self.rna.find_motif("GGG")

    def test_str(self):
        string = f'>{self.rna.identifier}\n{self.rna.data}'
        self.assertEqual(str(self.rna), string, "Błąd funkcji str()")
        string1 = f'>{self.rna1.identifier}\n{self.rna1.data}'
        self.assertEqual(str(self.rna1), string1, "Błąd funkcji str()")
        with self.assertRaises(ValueError):
            self.rna.data = ""
            print(self.rna)


class TestDNASequence(TestCase):

    def setUp(self):
        self.tab = "ATAgaTagga"
        self.tab1 = "aagaggatttaaa"
        self.dna = DNASequence(self.tab)
        self.dna1 = DNASequence(self.tab1)


    def test_create(self):

        self.assertEqual(self.dna.data, self.tab.upper(), "Blad tworzenia obiektu")
        self.assertEqual(self.dna.length, len(self.tab), "Blad tworzenia obiektu")
        self.assertEqual(self.dna.identifier, f'DNA{len(self.tab)} 5A0C2T3G', "Blad tworzenia obiektu")
        self.assertEqual(self.dna1.data, self.tab1.upper(), "Blad tworzenia obiektu")
        self.assertEqual(self.dna1.length, len(self.tab1), "Blad tworzenia obiektu")
        self.assertEqual(self.dna1.identifier, f'DNA{len(self.tab1)} 7A0C3T3G', "Blad tworzenia obiektu")
        with self.assertRaises(TypeError):
            RNASequence(3)
        with self.assertRaises(ValueError):
            RNASequence("XXX")

    def test_complement(self):
        self.assertEqual(self.dna.complement().data, "TCCTATCTAT", "Blad funkcji complement()")
        self.assertEqual(self.dna1.complement().data, "TTTAAATCCTCTT", "Blad funkcji complement()")
        self.assertEqual(DNASequence("").complement().data, "", "Blad funkcji complement()")
        with self.assertRaises(ValueError):
            DNASequence("XYZ").complement()

    def test_transcribe(self):
        # self.assertEqual(self.dna.transcribe().data, "UCCUAUCUAU", "Blad funkcji transcribe()1")
        self.assertEqual(self.dna1.transcribe().data, "UUUAAAUCCUCUU", "Blad funkcji transcribe()2")
        self.assertEqual(DNASequence("").transcribe().data, "", "Blad funkcji transcribe()3")
        with self.assertRaises(ValueError):
            DNASequence("XYZ").transcribe()

    def test_mutate(self):
        self.dna.mutate(1, 'G')
        self.assertEqual(self.dna.data, "AGAgaTagga".upper(), "Blad funkcji mutującej")
        with self.assertRaises(IndexError):
            self.dna.mutate(18, 'G')
        with self.assertRaises(IndexError):
            self.dna.mutate(10, 'G')
        with self.assertRaises(ValueError):
            self.dna.mutate(1, 'X')
        with self.assertRaises(ValueError):
            self.dna.mutate(1, '')
        with self.assertRaises(ValueError):
            self.dna.data = ""
            self.dna.mutate(1, 'G')

    def test_find_motif(self):
        self.assertEqual(self.dna.find_motif("A"), 0, "Blad funkcji find_motif")
        self.assertEqual(self.dna.find_motif("AtA"), 0, "Blad funkcji find_motif")
        self.assertEqual(self.dna.find_motif("GG"), 7, "Blad funkcji find_motif")
        self.assertEqual(self.dna.find_motif("gg"), 7, "Blad funkcji find_motif")
        with self.assertRaises(IndexError):
            self.dna.find_motif("aaaaaaaaaaaaaaaaa")
        with self.assertRaises(TypeError):
            self.dna.find_motif(2)
        with self.assertRaises(ValueError):
            self.dna.find_motif(None)
        with self.assertRaises(ValueError):
            self.dna.find_motif("GGG")

    def test_str(self):
        string = f'>{self.dna.identifier}\n{self.dna.data}'
        self.assertEqual(str(self.dna), string, "Błąd funkcji str()")
        string1 = f'>{self.dna1.identifier}\n{self.dna1.data}'
        self.assertEqual(str(self.dna1), string1, "Błąd funkcji str()")
        with self.assertRaises(ValueError):
            self.dna.data = ""
            print(self.dna)
