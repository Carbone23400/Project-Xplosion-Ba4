# give the information about mol from its name

import re

from .parser import KNOWN_LIGANDS
from .parser import ParsedComplex
from .parser import _enrich


def ligand_data(name: str) -> dict[str, int]:
    ligand = {}
    name = _normalize_name(name)
    #full_name=name
    for ligand_symbol, data in sorted(KNOWN_LIGANDS.items(), key=lambda item: len(item[1][0]), reverse=True):
        ligand_name = _normalize_name(data[0])
       # if ligand_name in name:
        #    for prefixe, number in PREFIXE.items():
         #       ligand[ligand_symbol]=number
              #  break
      #  else :
          #  continue
        found=False
        for prefixe, number in sorted(PREFIXE.items(), key=lambda x: len(x[0]), reverse=True):
            pattern_1 = prefixe + ligand_name
            if pattern_1 in name:
                ligand[ligand_symbol] = ligand.get(ligand_symbol, 0) + number
                name = name.replace(pattern_1, "", 1)
                found = True
                break
        if found:
            continue

        pattern_2 = ligand_name
        if pattern_2 in name:
            ligand[ligand_symbol] = ligand.get(ligand_symbol, 0) + 1
            name = name.replace(pattern_2, "", 1)
            continue
             
    return ligand


METALS_NAME: dict[str, tuple[str, str]] = {
    # symbol        : (cation_name, anion_name)
    "Li": ("lithium", "lithiate"),
    "Be": ("beryllium", "beryllate"),
    "Na": ("sodium", "sodiate"),
    "Mg": ("magnesium", "magnesate"),
    "Al": ("aluminum", "aluminate"),
    "K": ("potassium", "potassiate"),
    "Ca": ("calcium", "calciate"),
    "Sc": ("scandium", "scandate"),
    "Ti": ("titanium", "titanate"),
    "V": ("vanadium", "vanadate"),
    "Cr": ("chromium", "chromate"),
    "Mn": ("manganese", "manganate"),
    "Fe": ("iron", "ferrate"),
    "Co": ("cobalt", "cobaltate"),
    "Ni": ("nickel", "nickelate"),
    "Cu": ("copper", "cuprate"),
    "Zn": ("zinc", "zincate"),
    "Ga": ("gallium", "gallate"),
    "Rb": ("rubidium", "rubidiate"),
    "Sr": ("strontium", "strontiate"),
    "Y": ("yttrium", "yttrate"),
    "Zr": ("zirconium", "zirconate"),
    "Nb": ("niobium", "niobate"),
    "Mo": ("molybdenum", "molybdate"),
    "Tc": ("technetium", "technetate"),
    "Ru": ("ruthenium", "ruthenate"),
    "Rh": ("rhodium", "rhodate"),
    "Pd": ("palladium", "palladate"),
    "Ag": ("silver", "argentate"),
    "Cd": ("cadmium", "cadmate"),
    "In": ("indium", "indate"),
    "Sn": ("tin", "stannate"),
    "Cs": ("cesium", "cesiate"),
    "Ba": ("barium", "bariate"),
    "La": ("lanthanum", "lanthanate"),
    "Ce": ("cerium", "cerate"),
    "Pr": ("praseodymium", "praseodymate"),
    "Nd": ("neodymium", "neodymate"),
    "Sm": ("samarium", "samariate"),
    "Eu": ("europium", "europate"),
    "Gd": ("gadolinium", "gadolinate"),
    "Tb": ("terbium", "terbate"),
    "Dy": ("dysprosium", "dysprosate"),
    "Ho": ("holmium", "holmate"),
    "Er": ("erbium", "erbate"),
    "Tm": ("thulium", "thulate"),
    "Yb": ("ytterbium", "ytterbate"),
    "Lu": ("lutetium", "lutetate"),
    "Hf": ("hafnium", "hafnate"),
    "Ta": ("tantalum", "tantalate"),
    "W": ("tungsten", "tungstate"),
    "Re": ("rhenium", "rhenate"),
    "Os": ("osmium", "osmate"),
    "Ir": ("iridium", "iridate"),
    "Pt": ("platinum", "platinate"),
    "Au": ("gold", "aurate"),
    "Hg": ("mercury", "mercurate"),
    "Tl": ("thallium", "thallate"),
    "Pb": ("lead", "plumbate"),
    "Bi": ("bismuth", "bismuthate"),
    "Ac": ("actinium", "actinate"),
    "Th": ("thorium", "thorate"),
    "U": ("uranium", "uranate")
}      
                
ROMAN_NUMBER={"I":1,"II":2,"III":3,"IV":4,"V":5,"VI":6,"VII":7,"VIII":8}

PREFIXE={"di":2,"bi":2,"tri":3,"tetra":4,"penta":5,"hexa":6,"hepta":7,"octa":8,"bis":2,"tris":3,"tetrakis":4,"pentakis":5,"hexakis":6,"heptakis":7,"octakis":8}



def extract_complex_charge_from_name(name: str) -> int | None:
   match = re.search(r'\((.*?)\)', name)
   if not match:
        return None
   roman = match.group(1).upper()
   return ROMAN_NUMBER.get(roman)

def metal_data(name: str) -> str:
        name = _normalize_name(name)
        metals_by_longest_name = sorted(
            METALS_NAME.items(),
            key=lambda item: max(len(item[1][0]), len(item[1][1])),
            reverse=True,
        )
        for symbols, metal_names in metals_by_longest_name:
            cation_name = _normalize_name(metal_names[0])
            anion_name = _normalize_name(metal_names[1])
            if cation_name in name:
                 return symbols
            elif anion_name in name:
                return symbols
        raise ValueError(f"Could not identify a metal name in '{name}'.")
#le but est de transformer notre nom en formule car le programme sait faire avec la formule

def parse_name(name: str) -> ParsedComplex:
     raw_name= name
     metal= metal_data(name)
     oxydation_state= extract_complex_charge_from_name(name)
     ligands= ligand_data(name)

     result= ParsedComplex(metal=metal, ligands=ligands, complex_charge=0, counter_ions={},raw_formula=raw_name)
     _enrich(result)

     if oxydation_state is not None:
          result.oxidation_state=oxydation_state
          result.complex_charge=oxydation_state + result.total_ligand_charge
     return result


def _normalize_name(name: str) -> str:
     """Normalize a coordination compound name for simple substring matching."""
     return re.sub(r"[\s_-]+", "", name).lower()

#test
if __name__=="__main__":
     name = input ("complex name is:")
     complex_info=parse_name(name)
     print(complex_info)

