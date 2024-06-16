class Wielomian:
    """
    Klasa opisująca wielomian, na podstawie listy jego współczynników
    oraz zawierająca funkcje, pozwalające na operacje na wielomianach
    param wsp modyfikowalna lista współczynników, gdzie
    indeks elementu odpowiada potędze (element o indeksie 0 to wyraz wolny,
    element o indeksie 1 to wyraz z x^1)
    :param: wsp lista współczynników, gdzie
    indeks elementu odpowiada potędze (element o indeksie 0 to wyraz wolny,
    element o indeksie 1 to wyraz z x^1)
    @author Hanna Górecka
    """

    def __init__(self, wsp):
        self.wsp = wsp

    def Stopien(self):
        """
        Funkcja zwraca stopień wielomianu
        :return: Stopien wielomianu (liczba)
        """
        i = len(self.wsp) - 1
        if self.wsp != 0.0:
            return i
        while i > 0:
            i -= 1
            if self.wsp[i] != 0.0:
                return i
        return 0

    def toString(self):
        """
        Funkcja zwraca String, z wielomianem w formacie ax^n+bx^(n-1)+...+cx+d
        :return: String z rozpisanym wielomianem
        """
        i = len(self.wsp) - 1
        text = ""

        while i >= 0:
            if len(text) == 0 and i > 1 and self.wsp[i] != 0.0:
                text += (str(self.wsp[i]) + "*x^" + str(i))
            elif len(text) == 0 and i == 1 and self.wsp[i] != 0.0:
                text += str(self.wsp[i]) + "*x"
                if self.wsp[0] > 0.0:
                    text += "+" + str(self.wsp[0])
                elif self.wsp[0] < 0.0:
                    text += str(self.wsp[0])

            elif len(text) == 0 and i == 0:
                if self.wsp[0] == 0.0:
                    text += "0"
                else:
                    text += str(self.wsp[0])

            elif i == 1 and self.wsp[i] > 0.0 and self.wsp[0] == 0.0:
                text += ("+" + str(self.wsp[1]) + "*x")

            elif i == 1 and self.wsp[i] < 0.0 and self.wsp[0] == 0.0:
                text += (str(self.wsp[1]) + "*x")

            elif i == 1 and self.wsp[i] > 0.0 and self.wsp[0] == 0.0:
                text += ("+" + str(self.wsp[1]) + "*x")

            elif i == 1 and self.wsp[i] < 0.0 and self.wsp[0] == 0.0:
                text += (str(self.wsp[1]) + "*x")

            elif i == 0 and self.wsp[0] != 0.0:
                if self.wsp[0] > 0:
                    text += "+" + str(self.wsp[0])
                else:
                    text += str(self.wsp[0])

            elif self.wsp[i] > 0.0:
                text += ("+" + str(self.wsp[i]) + "*x^" + str(i))

            elif (self.wsp[i] < 0.0):
                text += (str(self.wsp[i]) + "*x^" + str(i))

            i -= 1
        return text

    def wartosc(self, x):
        """
        Funkcja oblicza wartość wielomianu dla podanego argumentu
        :param x: argument, dla którego ma być policzony wielomian
        :return: wartość wielomianu
        """
        sum = 0.0
        i = 0

        for it in self.wsp:
            sum += it*x**i
            i += 1
        return sum


    def __add__(self, wiel2):
        """
         Przeciążenie operatora + dla wielomianów
        :param wiel2: drugi wielomian, który ma być dodany do tego wielomianu
        :return: obiekt klasy Wielomian, będący sumą dwóch wielomianów
        """
        wspList = []
        n = 0
        if len(self.wsp) >= len(wiel2.wsp):
            n = len(wiel2.wsp)
        else:
            n = len(self.wsp)

        for i in range(0, n):
            wspList.insert(i, self.wsp[i] + wiel2.wsp[i])
        if n <= len(self.wsp) - 1:
            i = n
            while True:
                wspList.insert(i, self.wsp[i])
                i += 1
                if i <= len(self.wsp.size) - 1:
                    continue
                else:
                    break

        elif n <= len(wiel2.wsp) - 1:
            i = n
            while True:
                wspList.insert(i, wiel2.wsp[i])
                i += 1
                if i <= len(wiel2.wsp) - 1:
                    continue
                else:
                    break
        wiel3 = Wielomian(wspList)
        return wiel3

    def __sub__(self, wiel2):
        """
         Przeciążenie operatora - dla wielomianów
        :param wiel2: drugi wielomian, który ma być odjęty od tego wielomianu
        :return: obiekt klasy Wielomian, będący różnicą dwóch wielomianów
        """
        wspList = []
        n = 0
        if len(self.wsp) >= len(wiel2.wsp):
            n = len(wiel2.wsp)
        else:
            n = len(self.wsp)
        for i in range(0, n):
            wspList.insert(i, self.wsp[i] - wiel2.wsp[i])
        if n < len(self.wsp):
            i = n
            while True:
                wspList.insert(i, self.wsp[i])
                i+=1
                if i < len(self.wsp):
                    continue
                else:
                    break

        elif n < len(wiel2.wsp):
            i = n
            while True:
                wspList.insert(i, (-1) * wiel2.wsp[i])
                i += 1
                if i < len(wiel2.wsp):
                    continue
                else:
                    break
        wiel3 = Wielomian(wspList)
        return wiel3

    def __mul__(self, wiel2):
        """
        Przeciążenie operatora * dla wielomianów
        :param wiel2: drugi wielomian, przez który ma być pomnożony ten wielomian
        :return: obiekt klasy Wielomian, będący iloczynem dwóch wielomianów
        """
        wspList = []
        n = len(self.wsp)
        m = len(wiel2.wsp)
        p = 0
        w = 0.0
        for i in range(0,n):
            for j in range(0,m):
                p=i+j
                w=self.wsp[i] * wiel2.wsp[j]
                if len(wspList) > p:
                    wspList[p] += w
                else:
                    wspList.insert(p, w)
        wiel3 = Wielomian(wspList)
        return wiel3

    def __iadd__(self, wiel2):
        """
        Przeciążenie operatora += dla wielomianów
        :param wiel2: drugi wielomian, który ma być dodany do tego wielomianu
        :return: obiekt klasy Wielomian, będący sumą dwóch wielomianów
        """
        return self + wiel2

    def __isub__(self, wiel2):
        """
        Przeciążenie operatora -= dla wielomianów
        :param wiel2: drugi wielomian, który ma być odjęty od tego wielomianu
        :return: obiekt klasy Wielomian, będący różnicą dwóch wielomianów
        """
        return self - wiel2
    def __imul__(self, wiel2):
        """
        Przeciążenie operatora *= dla wielomianów
        :param wiel2: drugi wielomian, przez który ma być pomnożony ten wielomian
        :return: obiekt klasy Wielomian, będący iloczynem dwóch wielomianów
        """
        return self*wiel2



tab = [1,1,1,1]
tab2 = [2,2,2,2]
wiel = Wielomian(tab)
wiel2 = Wielomian(tab2)
v = wiel.wartosc(2)
v2 = wiel2.wartosc(2)

print("\n" + str(v))
print("\n" + str(v2))
wiel3 = wiel*wiel2
print("\n" + wiel3.toString())
wiel*=wiel2
print("\n" + wiel.toString())