from math import sqrt
from Bio.PDB.PDBList import PDBList
from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--name", type=str, help="Structure name", required=False, default="4YWO")
args = parser.parse_args()
structure_name = args.name

def calculate_distance(structure):
    ca_atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                for atom in res.get_atoms():
                    if atom.get_name() == 'CA':
                        ca_atoms.append(atom.get_coord())
    num_atoms = len(ca_atoms)
    contact_map = np.zeros((num_atoms, num_atoms), dtype=int)
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            distance = sqrt(sum((ca_atoms[i] - ca_atoms[j]) ** 2))
            if distance < 8.0:
                contact_map[i][j] = 1
                contact_map[j][i] = 1
    return contact_map

def map_plot(contact_map):
    plt.imshow(contact_map, cmap='binary', interpolation='nearest')
    plt.title(f'Contact Map for {structure_name}')
    plt.xlabel('Atom Index')
    plt.ylabel('Atom Index')
    plt.show()

try:
    pdbl = PDBList()
    fetch_pdb = pdbl.retrieve_pdb_file(structure_name, file_format='pdb')
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure(structure_name, fetch_pdb)
    contact_map = calculate_distance(structure)
    map_plot(contact_map)
except FileNotFoundError:
    print("Podano złą nazwę pliku lub istnieje problem z pobraniem pliku")
except:
    print("Wystąpił inny błąd")
