from src.coordchem.parser import ParsedComplex
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

drawit(aspirin_3d)
