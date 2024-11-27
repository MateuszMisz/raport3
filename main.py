import Bio.PDB
import Bio.PDB.PDBIO
from Bio.PDB.PDBList import PDBList
from Bio.PDB import PDBParser
from Bio.PDB import Residue,Structure,Atom,Model,Chain
import Bio
import Bio
from typing import List,Tuple,Dict
import argparse
representing_atoms_table={
        'purine':['N9','C2','C6'],
        'pirmidine':['N1','C2','C4'],
        'backbone':['P',"C4'"]
    }
def handle_arguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('structure_id',type=str,help='id of structure. when no --load_from_disc specified, loads .\\structure_id.pdb or downloads from rcsb')
    parser.add_argument('--load_from_disc','-l',type=str,help='file to load structure from if it already exists, but relative path isn\'t: .\\structure_id.pdb')

def load_from_RCSB(structure_id='430d'):
    pdb_list=PDBList()
    fetched_structure=pdb_list.retrieve_pdb_file(structure_id,file_format='pdb')
    parser=PDBParser()
    original_structure=parser.get_structure(structure_id,fetched_structure)
    return original_structure
def pirmidine_or_purine(residue:Residue.Residue):
    if residue.get_resname() in ['A','G']:
        return 'purine'
    elif residue.get_resname() in['C','T','U']:
        return 'pirmidine'
    else :return None
def residue_coarse_grain(residue:Residue.Residue,representing_atoms_table:dict):
    new_residue=Residue.Residue(residue.id,residue.resname,residue.segid)
    pir_or_pur=pirmidine_or_purine(residue)
    
    for atom in residue:
        print(f'{atom.name} {residue.resname} {pir_or_pur} {representing_atoms_table[pir_or_pur] if pir_or_pur is not None else ""}')
        if atom.name in representing_atoms_table['backbone']:
            new_residue.add(atom)   
        elif pir_or_pur is not None:
            if atom.name in representing_atoms_table.get(pir_or_pur):
                new_residue.add(atom)
    
    return new_residue
def get_coarse_grain_structure(structure:Structure.Structure,representing_atoms_table)->Tuple[Atom.Atom]:   
    new_structure=Structure.Structure(structure.id)
    for model in structure:
        new_structure.add(new_model:=Model.Model(model.id,model.serial_num))
        for chain in model:
            new_model.add(new_chain:=Chain.Chain(chain.id))
            for residue in chain:
                if residue.get_id()[0]!=' ':  ##if is not ' ' then its not a part of molecule probably
                    continue
                new_residue=residue_coarse_grain(residue,representing_atoms_table)
                new_chain.add(new_residue)
    return new_structure
def save_structure(structure:Structure.Structure,file_name:str|None=None):
    """default fileneame: f'{structure.id}_coarse_grained.pdb'"""
    if file_name is None:
        file_name=str(structure.id)+'_coarse_grained.pdb'
    pdb_io=Bio.PDB.PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(file_name)

original_structure=load_from_RCSB()
coarse_grained_structure=get_coarse_grain_structure(original_structure,representing_atoms_table)
save_structure(coarse_grained_structure)