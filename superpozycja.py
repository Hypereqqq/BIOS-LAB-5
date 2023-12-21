from Bio.PDB import *
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np
import math

def odleglosc(koordynaty_1, koordynaty_2):
    return sum((k1 - k2)**2 for k1, k2 in zip(koordynaty_1, koordynaty_2))

def RMSD(koordynaty_referencyjny, koordynaty_modele):
    rmsd = math.sqrt(sum(odleglosc(ref, model) for ref, model in zip(koordynaty_referencyjny, koordynaty_modele)) / len(koordynaty_referencyjny))
    return(rmsd)

def superpozycja(koordynaty_referencyjny, koordynaty_modele):
    superpozycja = SVDSuperimposer()
    superpozycja.set(koordynaty_referencyjny,koordynaty_modele)
    superpozycja.run()
    return superpozycja.get_transformed()

def pobierz_atomy(reszty_modele, reszty_referencyjny):
    pierwsza_iteracja = True
    for reszta in reszty_referencyjny.get_residues():
        for atom in reszta.get_atoms():
            koordynaty = reszty_modele[reszta.id[1]][atom.get_id()].get_coord()
            if pierwsza_iteracja:
                koordynaty_modele = np.array(koordynaty)
                pierwsza_iteracja = False
            else:
                koordynaty_modele = np.vstack([koordynaty_modele,koordynaty])
    return(koordynaty_modele)

parser = PDBParser()
pliki_z_modelami = parser.get_structure("Modele", "R1107TS081.pdb")
plik_z_referencyjnym = parser.get_structure("Referencyjny", "R1107_reference.pdb")
reszty_referencyjny = plik_z_referencyjnym[0]["0"]
koordynaty_referencyjny = pobierz_atomy(reszty_referencyjny, reszty_referencyjny)

f = open("wynik.txt", "w")
for model in pliki_z_modelami.get_models():
    reszty_modele = model["0"]
    koordynaty_modele = pobierz_atomy(reszty_modele, reszty_referencyjny)
    superpozycja_wynik = superpozycja(koordynaty_referencyjny, koordynaty_modele)
    rmsd = RMSD(koordynaty_referencyjny, superpozycja_wynik)
    f.write("Model numer " + str(model.id + 1) + "  RMSD: \n" + str(rmsd) + "\n")
