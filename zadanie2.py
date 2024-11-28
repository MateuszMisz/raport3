from math import sqrt
from Bio.PDB.PDBList import PDBList
from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB.vectors import calc_dihedral
from Bio.PDB.vectors import Vector
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--name", type=str, help="Structure name", required=False, default="4YWO")
parser.add_argument("--function", type=str, help="Function used to create the plot", required=False, default="calc_dihedral", choices=["get_phi_psi_list", "calc_dihedral"])
args = parser.parse_args()
structure_name = args.name

def calculate_phi_psi(structure):
    phi_psi_angles = []
    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            for poly in polypeptides:
                phi_psi_angles.extend(poly.get_phi_psi_list())
    return phi_psi_angles

def calculate_phi_psi_dihedral(structure):
    phi_psi_angles = []
    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            for poly in polypeptides:
                for i in range(1, len(poly) - 1):
                    res_prev = poly[i - 1]
                    res_current = poly[i]
                    res_next = poly[i + 1]
                    try:
                        c_prev = Vector(res_prev['C'].get_coord())
                        n = Vector(res_current['N'].get_coord())
                        ca = Vector(res_current['CA'].get_coord())
                        c = Vector(res_current['C'].get_coord())
                        phi = calc_dihedral(c_prev, n, ca, c)
                        n_next = Vector(res_next['N'].get_coord())
                        psi = calc_dihedral(n, ca, c, n_next)
                        phi_psi_angles.append((phi, psi))
                    except KeyError:
                        continue
    return phi_psi_angles

def Ramachandran_plot(angles):
    phi_angles = []
    psi_angles = []
    for phi, psi in angles:
        if phi is not None and psi is not None:
            phi_angles.append(np.degrees(phi))
            psi_angles.append(np.degrees(psi))
    plt.figure(figsize=(8, 6))
    plt.scatter(phi_angles, psi_angles, color='magenta', alpha=0.5)
    plt.title(f'Wykres Ramachandrana dla struktury {structure_name}')
    plt.xlabel('Kąt Phi (°)')
    plt.ylabel('Kąt Psi (°)')
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.grid(True)
    plt.show()

try:
    pdbl = PDBList()
    fetch_pdb = pdbl.retrieve_pdb_file(structure_name, file_format='pdb')
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure(structure_name, fetch_pdb)
    if args.function == "calc_dihedral":
        angles = calculate_phi_psi_dihedral(structure)
    if args.function == "get_phi_psi_list":
        angles = calculate_phi_psi(structure)
    Ramachandran_plot(angles)
except FileNotFoundError:
    print("Podano złą nazwę pliku lub istnieje problem z pobraniem pliku")
except:
    print("Wystąpił inny błąd")
