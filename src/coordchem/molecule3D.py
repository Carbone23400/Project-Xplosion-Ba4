from coordchem.parser import ParsedComplex
from rdkit.Chem import AllChem 
from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd
from pathlib import Path
import os


#pour la molecule entrée, il faut qu'il donne le smile et apres la conformation 3d: 

#il faut d'abbord ajouter les H :
#molecule = Chem.AddHs(molecule)
#AllChem.EmbedMolecule(molecule)
#molecule


import py3Dmol

# Helpful function from a blog (https://rdkit.blogspot.com/2016/07/a-recent-post-on-in-pipeline-talked.html)
def drawit(m,p=None,confId=-1):
        mb = Chem.MolToMolBlock(m,confId=confId) # Returns SDF file as a string
        if p is None:
            p = py3Dmol.view(width=400,height=400)
        p.removeAllModels()
        p.addModel(mb,'sdf') 
        p.setStyle({'stick':{}})
        p.setBackgroundColor('white')
        p.zoomTo()
        return p.show()

#créer une molecule rdkit pour un ligand en ayant deja le smile:
def ligand_3D(ligand_smile):
      ligand_mol=Chem.MolFromSmiles(ligand_smile)
      if ligand_mol is None:
        raise ValueError("Invalid SMILES string.")
      ligand_mol=Chem.AddHs(ligand_mol)
      ligand_3D=AllChem.EmbedMolecule(ligand_mol)
      if ligand_3D != 0:
        raise ValueError("3D embedding failed.")
      return ligand_mol

if __name__ == "__main__":
    thalidomide_smiles = "O=C1c2ccccc2C(=O)N1[C@H]3CCC(=O)NC3=O"

    thalidomide_3d = ligand_3D(thalidomide_smiles)

    drawit(thalidomide_3d)
print(thalidomide_3d.GetNumConformers())