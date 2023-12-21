from Bio.PDB import *
import math

radii = {"O": 1.52, "P": 1.8, "N": 1.55, "H": 1.2, "C": 1.7}

def dystans(k1, k2):
    return math.sqrt(sum((k1 - k2) ** 2 for k1, k2 in zip(k1, k2)))

def liczenie_atomow(lancuch):
    return sum(1 for _ in lancuch.get_atoms())

def clash_miedzy_resztami(reszta1, reszta2, prog):
    ilosc_clash = 0
    for atom1 in reszta1.get_atoms():
        for atom2 in reszta2.get_atoms():
            dystans = dystans(atom1.coord, atom2.coord)
            if dystans <= radii[atom1.element] + radii[atom2.element] - prog:
                ilosc_clash += 1
    return ilosc_clash

def clashScore(plik, prog):
    ilosc_clash = 0
    parser = PDBParser()
    
    struktura = parser.get_structure("structure", plik)
    model = next(struktura.get_models())
    lancuch = next(model.get_chains())

    for i, reszta in enumerate(lancuch):
        for j in range(i + 2, len(lancuch)):
            ilosc_clash += clash_miedzy_resztami(lancuch[i + 1], lancuch[j + 1], prog)

    return 1000 * (ilosc_clash / liczenie_atomow(lancuch))

plik = "model_4.pdb"
prog = 0.4

wynik = clashScore(plik, prog)
print(f"Clashscore: {round(wynik, 5)}")
