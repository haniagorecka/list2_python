from unittest import TestCase

from Wielomian import Wielomian


class TestWielomian(TestCase):
    wsp = [1, 2, 3, 4, 5]
    wsp2 = [1.0, 2.0, 3.0]
    wsp3 = [-1, 0, 1, -5]
    wsp4 = [0, 0, 0]
    wsp5 = [0, 0, 1]
    wsp6 = [1, 0, 0]
    wiel = Wielomian(wsp)
    wiel2 = Wielomian(wsp2)
    wiel3 = Wielomian(wsp3)
    wiel4 = Wielomian(wsp4)
    wiel5 = Wielomian(wsp5)
    wiel6 = Wielomian(wsp6)

    def test_stopien(self):
        self.assertEqual(self.wiel.Stopien(), 4, "Błąd funkcji stopien()1")
        self.assertEqual(self.wiel2.Stopien(), 2, "Błąd funkcji stopien()2")
        self.assertEqual(self.wiel3.Stopien(), 3, "Błąd funkcji stopien()3")
        self.assertEqual(self.wiel4.Stopien(), 0, "Błąd funkcji stopien()4")
        self.assertEqual(self.wiel5.Stopien(), 2, "Błąd funkcji stopien()5")
        self.assertEqual(self.wiel6.Stopien(), 0, "Błąd funkcji stopien()6")

    def test_wartosc(self):
        self.assertEqual(self.wiel.wartosc(1), 15, "Błąd funkcji wartosc()1")
        self.assertEqual(self.wiel.wartosc(2), 129, "Błąd funkcji wartosc()2")
        self.assertEqual(self.wiel.wartosc(0), 1, "Błąd funkcji wartosc()3")
        self.assertEqual(self.wiel.wartosc(0.0), 1, "Błąd funkcji wartosc()4")
        self.assertEqual(self.wiel4.wartosc(11), 0, "Błąd funkcji wartosc()5")
        self.assertEqual(self.wiel.wartosc(-1), 3, "Błąd funkcji wartosc()6")

    def test_str(self):
        self.assertEqual(str(self.wiel), "5*x^4+4*x^3+3*x^2+2*x^1+1", "Błąd funkcji str()1")
        self.assertEqual(str(self.wiel6), "1", "Błąd funkcji str()2")
        self.assertEqual(str(self.wiel2), "3.0*x^2+2.0*x^1+1.0", "Błąd funkcji str()3")
        self.assertEqual(str(self.wiel3), "-5*x^3+1*x^2-1", "Błąd funkcji str()4")
        self.assertEqual(str(self.wiel4), "0", "Błąd funkcji str()5")
        self.assertEqual(str(self.wiel5), "1*x^2", "Błąd funkcji str()1")

    def test_plus(self):
        wiel7 = self.wiel+self.wiel3
        self.assertEqual(wiel7.wsp, [0, 2, 4, -1, 5], "Błąd funkcji __plus__")
        self.assertEqual((self.wiel+self.wiel4).wsp, self.wiel.wsp, "Błąd funkcji __plus__")
        self.wiel6 += self.wiel5
        self.assertEqual(self.wiel6.wsp, [1, 0, 1], "Błąd funkcji __plus__")

    def test_sub(self):
        wiel7 = self.wiel - self.wiel3
        self.assertEqual(wiel7.wsp, [2, 2, 2, 9, 5], "Błąd funkcji __sub__")
        self.assertEqual((self.wiel - self.wiel4).wsp, self.wiel.wsp, "Błąd funkcji __sub__")
        self.wiel6 -= self.wiel5
        self.assertEqual(self.wiel6.wsp, [1, 0, -1], "Błąd funkcji __sub__")

    def test_mul(self):
        wiel7 = self.wiel * self.wiel3
        self.assertEqual(wiel7.wsp, [-1, -2, -2, -7, -12, -11, -15, -25], "Błąd funkcji __mul__")
        self.assertEqual((self.wiel*self.wiel4).wsp, [], "Błąd funkcji __mul__")
        self.wiel5 *= self.wiel6
        self.assertEqual(self.wiel5.wsp, self.wiel5.wsp, "Błąd funkcji __mul__")
