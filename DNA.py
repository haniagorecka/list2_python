import re


class BaseClass:
    """
    Klasa bazowa zawierająca sekwencję zasad/aminokwasów oraz dotyczące jej funkcje
    :param data: sekwencja zasad/aminokwasów
    Attributes:
        VALID_CHARS: Lista doswolonych znaków (zasad/aminokwasów)
        identifier: String identyfikujący sekwencję, będzie tworzony w dziedziczących klasach
        length: Długość sekwencji
    :author: Hanna Górecka
    """
    VALID_CHARS = []
    identifier = ""
    length = 0

    def __init__(self, data):
        """
        Inicjalizacja obiektu, definiuje właściwość length jako długość data,
        ustawia data jako właściwość klasy
        :param data: sekwencja zasad/aminokwasów
        """
        if type(data) is str or all(isinstance(x, str) for x in data):
            pass
        else:
            raise TypeError("Błędny typ danych")
        if type(data) is str:
            self.data = data.upper()
        else:
            self.data = data
        self.length = len(self.data)
        for x in self.data:
            if not self.VALID_CHARS.__contains__(x):
                raise ValueError("Błędne dane w sekwencji")

    def mutate(self, position, val):
        """
        Funkcja zmienia element o podanym indeksie na podany nowy element
        Sprawdza czy nić nie jest pusta, czy istnieje pozycja o danym indeksie i czy podana wartość jest dozwolona
        Zmienia wartość data na string po zmianie wartości
        :param position: indeks, na którym ma być zmieniona zasada
        :param val: element, na który ma być zmieniony element na indeksie o pozycji position
        """
        if self.data == "" or self.data is None:
            raise ValueError("Brak danych")
        if self.length <= position:
            raise IndexError("Brak takiej pozycji")
        if not self.VALID_CHARS.__contains__(val):
            raise ValueError("Błędny argument")

        if type(self.data) is str:
            s = list(self.data)
            s[position] = val
            string = ''.join(s)
            self.data = string
        else:
            self.data[position] = val
        self.__init__(self.data)

    def __str__(self):
        """
        Funkcja zwraca String zawierający opis sekwencji DNA w formacie FASTA jako:
        >identyfier
        data
        :return: String, opisujący sekwencję DNA/RNA/Aminokwasów
        """
        if self.data == "" or self.data is None:
            raise ValueError("Brak danych")
        if type(self.data) is str:
            return ">" + self.identifier + "\n" + self.data
        else:
            s = ">" + self.identifier + "\n"
            for el in self.data:
                s = s + el + ""
            return s

    def find_motif(self, motif):
        """
         Funkcja zwraca indeks, na którym znajduje się pierwszy element
        podanego jako argument motywu (sekwencji)
        Jeśli takiego motywu nie ma, wyrzucony jest błąd
        Wspomogłam się odpowiedzią użytkownika Rik Poggi z linku:
        https://stackoverflow.com/questions/10459493/find-indexes-of-sequence-in-list-in-python
        :param motif: motyw, szukany w data
        :return: indeks, na którym znajduje się pierwszy element motywu
        """
        if motif == "" or motif is None:
            raise ValueError("Nie podano motywu")
        if type(motif) is str or all(isinstance(x, str) for x in motif):
            pass
        else:
            raise TypeError("Błędny typ danych")
        if len(motif) > self.length:
            raise IndexError("Zbyt dlugi motyw")
        if type(motif) is str:
            for i in re.finditer(motif.upper(), self.data):
                return i.start()
        else:
            for i in range(len(self.data)):
                if self.data[i:i + len(motif)] == motif:
                    return i
        raise ValueError("Nie znaleziono szukanego motywu")


class ProteinSequence(BaseClass):
    """
        Klasa dziedzicząca po BaseClass
        Klasa opisuje sekwencję aminokwasów i zawiera metody jej dotyczące
        Attributes:
            identifier: String identyfikujący sekwencję, składa się z tekstu:
            "ProtSeq"liczba_aminokwasów pierwsze_litery_aminokwasów
    """

    def __init__(self, data):
        self.VALID_CHARS = ["Phe", "Leu", "Ile", "Met",
                            "Val", "Ser", "Pro", "Thr", "Ala", "Tyr", "His",
                            "Gln", "Asn", "Lys", "Asp", "Glu", "Cys", "Trp",
                            "Arg", "Gly"]
        super().__init__(data)
        self.identifier = "ProteinSeq" + str(self.length) + " "
        for el in self.data:
            self.identifier += list(el)[0]


class RNASequence(BaseClass):
    """
        Klasa dziedzicząca po BaseClass
        Klasa opisuje sekwencję RNA i zawiera metody jej dotyczące
        Attributes:
            identifier: String identyfikujący sekwencję, składa się z tekstu:
            "RNA"liczba_zasad liczba_adeninaA_liczba_cytozynaC_liczba_uracylU_liczba_guaninaG
     """

    def __init__(self, data):
        self.VALID_CHARS = ['A', 'C', 'U', 'G']
        super().__init__(data)
        self.identifier = "RNA" + str(self.length) + " "
        a = 0
        c = 0
        u = 0
        g = 0
        for el in self.data:
            match el:
                case 'A':
                    a += 1
                case 'C':
                    c+=1
                case 'U':
                    u+=1
                case 'G':
                    g+=1
                case _ :
                    raise ValueError("Błędny znak w nici")
        self.identifier+=f'{a}A{c}C{u}U{g}G'

    def translate(self):
        slownik_aminokwasow = {
            "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
            "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
            "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
            "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
            "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
            "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
            "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
            "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
            "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
            "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
            "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
            "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
            "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
            "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
            "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
            "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
        }
        i = 0
        aminokwasy = []
        while i + 3 <= len(self.data):
            temp = self.data[i:i + 3]
            if temp in slownik_aminokwasow.keys():
                if slownik_aminokwasow.get(temp) == "Stop":
                    break
                else:
                    aminokwasy.append(slownik_aminokwasow.get(temp))
                    i += 3
            else:
                raise ValueError("Niewłaściwy kodon w nici RNA")
        return ProteinSequence(aminokwasy)


class DNASequence(BaseClass):
    """
    Klasa dziedzicząca po BaseClass
    Klasa opisuje sekwencję DNA i zawiera metody jej dotyczące
    Attributes:
        identifier: String identyfikujący sekwencję, składa się z tekstu:
            "DNA"liczba_zasad liczba_adeninaA_liczba_cytozynaC_liczba_tyminaT_liczba_guaninaG
    """

    def __init__(self, data):
        self.VALID_CHARS = ['A', 'C', 'T', 'G']
        super().__init__(data)
        self.identifier = "DNA" + str(self.length) + " "
        a = 0
        c = 0
        t = 0
        g = 0
        for el in self.data:
            match el:
                case 'A':
                    a += 1
                case 'C':
                    c += 1
                case 'T':
                    t += 1
                case 'G':
                    g += 1
                case _:
                    raise ValueError("Błędny znak w nici")
        self.identifier += f'{a}A{c}C{t}T{g}G'

    def complement(self):
        """
            Funkcja zwraca nić komplementarną do nici kodującej DNA
            :author Hanna Górecka
            :return:komplementarna nić matrycowa (odwrócona kolejność)
            :raise ValueError jeśli znajdzie niewłaściwy znak w nici kodującej
        """
        komplementarna = ""
        for el in self.data:
            match el:
                case 'A':
                    komplementarna += 'T'
                case 'C':
                    komplementarna += 'G'
                case 'T':
                    komplementarna += 'A'
                case 'G':
                    komplementarna += 'C'
                case _:
                    raise ValueError(" Niewłaściwy symbol w nici kodujacej")
        return DNASequence(komplementarna[::-1])

    def transcribe(self):
        """
            Funkcja zwraca nić transkrybowaną RNA do nici matrycowej
            :author: Hanna Górecka
            :return: transkrybowana nić RNA
            :raise ValueError jeśli znajdzie niewłaściwy znak w nici matrycowej
            """
        transkrybowana = ""
        for el in self.data:
            match el:
                case 'A':
                    transkrybowana += 'U'
                case 'C':
                    transkrybowana += 'G'
                case 'T':
                    transkrybowana += 'A'
                case 'G':
                    transkrybowana += 'C'
                case _:
                    raise ValueError(" Niewłaściwy symbol w nici matrycowej")

        return RNASequence(transkrybowana[::-1])


kod = "ACTGA"
dna = DNASequence(kod)
print(f'Obiekt DNA dla nici ACTGA: \n{dna}')
print(f'Indeks na ktorym znajduje sie motyw GA: {dna.find_motif("GA")}')
dna.mutate(1, 'A')
print(f'Obiekt DNA zmutowany zasadą A na 1 indeksie \n{dna}')
print(f'Obiekt DNA z nicia komplementarna: \n{dna.complement()}')
rna = dna.transcribe()
print(f'Obiekt RNA powstaly jako transkrybcja DNA dla nici ACTGA:\n {rna}')
print(f'Indeks na ktorym znajduje sie motyw UU: {rna.find_motif("UU")}')
rna.mutate(1, 'U')
print(f'Obiekt RNA zmutowany zasadą U na 1 indeksie \n{rna}')
prot = rna.translate()
print(f'Obiekt ProteinSequence powstaly jako translacja RNA: \n{prot}')
prot1 = ProteinSequence(['Leu', 'Phe', 'Ala', 'Ala', 'Lys'])
print(f'Obiekt ProteinSequence: \n{prot1}')
prot1.mutate(2, 'Glu')
print(f'Obiekt ProteinSequence zmutowany aminokwasem Glu na 2 indeksie \n{prot1}')
print(f'Indeks na ktorym znajduje sie motyw  Phe Glu Ala: {prot1.find_motif(['Phe', 'Glu', 'Ala'])}')

