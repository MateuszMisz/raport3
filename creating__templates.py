from Bio.PDB import PDBParser,PDBIO
from Bio.PDB import Structure,Chain,Model,Residue
import Bio
import argparse

arg_parser=argparse.ArgumentParser()
arg_parser.add_argument('input_file',type=str)
arg_parser.add_argument('--symbol',type=str)
arg_parser.add_argument('--res_id',type=str)
arg_parser.add_argument('--part','-p',type=str,choices=['backbone','base','all'],help='which part of residue gets added to the tamplate')
arg_parser.add_argument('--output','-o',type=str,required=True)
arg_parser.add_argument('--oid',type=str)
args=arg_parser.parse_args()

parser=PDBParser()
og_structure=parser.get_structure('ogs',args.input_file)
structure_extracted=False
for model in og_structure:
    for chain in model:
        for residue in chain:
            if (residue.id==args.res_id or residue.resname==args.symbol) and not residue.get_id()[0]!=' ' and not structure_extracted:
                new_res=Residue.Residue(residue.id,residue.resname,residue.segid)
                for atom in residue:
                    if args.part=='backbone':
                        if atom.name[-1]=="'"or atom.name in ['P','OP1','OP2']:
                            new_res.add(atom)
                    elif args.part=='base':
                        if atom.name[-1]!="'" and atom.name not in ['P','OP1','OP2']:
                            new_res.add(atom)
                    elif args.part=='all':
                        new_res.add(atom)
                structure_extracted=True
new_structure=Structure.Structure('structure')
new_model=Model.Model('model')
new_chain=Chain.Chain('A')
new_chain.add(new_res)
new_model.add(new_chain)
new_structure.add(new_model)
saver=PDBIO()
saver.set_structure(new_structure)
saver.save(args.output)
