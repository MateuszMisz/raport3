from Bio.PDB import Residue,Model,Chain,Atom
from Bio.PDB.Structure import Structure
from Bio.PDB import PDBParser,Superimposer,PDBIO
from typing import Dict,List,Type,Any
import argparse
from copy import deepcopy
def copy_res(residue:Residue.Residue,id,segid):
    new_residue=Residue.Residue(id,residue.resname,segid)
    for atom in residue:
        new_residue.add(deepcopy(atom))
    return new_residue
def pirmidine_or_purine(residue:Residue.Residue):
        if residue.get_resname() in ['A','G']:
            return 'purine'
        elif residue.get_resname() in['C','T','U']:
            return 'pirmidine'
        else :return None
def handle_arguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('input',type=str,help='file with coarse_grained structure')
    parser.add_argument('--output_file',type=str,help='output file_path. if no is given, pdb.file is printed on stdout')
    parser.add_argument('--oid',type=str,help='id of output structure. default is the prefix of input file +"_rebuild"')
    parser.add_argument('--adenine_templates','-A',nargs='+',type=str,help='filepath to pdb file(s) with adenine template(s)')
    parser.add_argument('--cytosine_templates','-C',nargs='+',type=str,help='filepath to pdb file(s) with cytosine template(s)')
    parser.add_argument('--uracyl_templates','-U',nargs='+',type=str,help='filepath to pdb file(s) with uracyl template(s)')
    parser.add_argument('--guanine_templates','-G',nargs='+',type=str,help='filepath to pdb file(s) with guanine template(s)')
def get_first_residue_from_structure(structure:Structure)->Residue.Residue:
    for model in structure: 
        for chain in model:
            for residue in chain:
                return residue
def load_structure(file_path:str,id:str='coarse_grained'):
    parser=PDBParser()
    coarse_grain_structure=parser.get_structure(id,file_path)
    return coarse_grain_structure
def load_templates(adenines):
    pass
class CoarseGrainRebuilder:
    def __init__(self,templates:Dict[str,List[Structure]]|Dict[str,List[str]]) -> None:
        """templates should contain all residues and backbone existing in file.
        for each residue there can be more then one template(all corresponding
        templates will be checked for each residue looking for the best RMSD)
        .it must be either, {name_of_residue : list_of_structure_objects} 
        or {name_of_residue:list_of_filepaths}(if the structures are not already loaded).
        """
        self.parser=PDBParser()
        #type_of_templates=self.check_templates_values_type(templates)
        type_of_templates='str'
        if type_of_templates=='structure':
            self.templates=templates
        elif type_of_templates=='str':
            self.templates=self.load_templates(templates)
        self.representing_atoms={
            
        }
    
    def mutual_atom_names(self, res1:Residue.Residue,res2:Residue.Residue):  
        res1_atom_names={atom.name for atom in res1}
        res2_atom_names={atom.name for atom in res2}
        intersecting_atom_names=res1_atom_names.intersection(res2_atom_names)
        return intersecting_atom_names
    def rebuild_element(self,residue:Residue.Residue,backbone:bool=False):
        pir_or_pur=pirmidine_or_purine(residue)
        if backbone:
            template=get_first_residue_from_structure(self.templates['backbone'][0])
        else:
            template=get_first_residue_from_structure(self.templates[residue.resname][0])
        new_res=copy_res(template,residue.id,residue.segid)
        intersecting_atom_names=self.mutual_atom_names(residue,new_res)
        fixed=[atom for atom in residue if atom.name in intersecting_atom_names]
        moving=[atom for atom in new_res if atom.name in intersecting_atom_names]
        super_imposer=Superimposer()
        super_imposer.set_atoms(fixed=fixed,moving=moving)
        super_imposer.apply(new_res.get_atoms())
        return new_res
    #def check_templates_values_type(self,templates:Dict[Any,Any])->Type:
    #    templates_value_type='structure'
    #    for value in templates.values():
    #        if not isinstance(value,Structure):
    #            templates_value_type=type(None)
    #            break
    #    if(templates_value_type=='structure'):
    #        return templates_value_type
    #    templates_value_type='str'
    #    for value in templates.values():
    #        if not isinstance(value,str):
    #            templates_value_type=type(None)
    #            break
    #    if templates_value_type=='str':
    #        return templates_value_type
    #    raise ValueError('all values in templates must be a list of either a str or a Bio.PDB.Structure. all values and values elements must be of same type. at least one element must be in each list')
    def load_templates(self,templates:Dict[str,List[str]]):
        tmp_templates=dict()
        i=0
        for key in templates.keys():
            tmp_templates.update({key:[]})
            tmp_templates[key].append(self.parser.get_structure(str(i),templates.get(key)[0]))
            i=i+1
        return tmp_templates
    def rebuild_residue(self,residue:Residue.Residue):
        nitrogen_base=self.rebuild_element(residue)
        backbone=self.rebuild_element(residue,backbone=True)
        new_residue=nitrogen_base
        for atom in backbone.get_atoms():
            new_residue.add(atom)
        return new_residue
    def rebuild_structure(self,coarse_grain_structure:Structure,output_id:str='test')->Structure:
        i=0
        rebuilded_structure=Structure(output_id)
        for model in coarse_grain_structure:
            new_model=Model.Model(model.id,model.serial_num)
            rebuilded_structure.add(new_model)
            for chain in model:
                new_chain=Chain.Chain(chain.id)
                new_model.add(new_chain)
                for residue in chain:
                    i+=1
                    print(i)
                    new_residue=self.rebuild_residue(residue)
                    new_chain.add(new_residue)
        return rebuilded_structure
if __name__=='__main__':

    templates={'A':['adenine1.pdb'],'C':['cytozine1.pdb'],'G':['guanine1.pdb'],'U':['uracyl1.pdb'],'backbone':['backbone1.pdb']}

    rebuilder=CoarseGrainRebuilder(templates)
    #for key,value in rebuilder.templates.items():
    #    for model in value[0]:
    #        for chain in model:
    #            for res in chain:
    #                for atom in res:
    #                    print(f'{atom.name},{res.resname},{res.id}')
    #        print('\n\n')       
    structure=load_structure('430d_coarse_grained.pdb')            
    new_structure=rebuilder.rebuild_structure(structure)
    saver=PDBIO()
    saver.set_structure(new_structure)
    saver.save('testowa_odbudowa.pdb')
    print('dupa')