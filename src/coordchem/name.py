# give the information about mol from its name

#run with python -m coordchem.name

import re

from .parser import COUNTER_IONS
from .parser import COUNTER_ION_NAMES
from .parser import KNOWN_LIGANDS
from .parser import ParsedComplex
from .parser import _apply_ambidentate_donor_assignments
from .parser import _enrich
from .parser import get_iupac_name


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
   for content in re.findall(r"\(([^)]*)\)", name):
        token = content.strip().upper()
        if token == "0":
            return 0
        if token in ROMAN_NUMBER:
            return ROMAN_NUMBER[token]
   return None

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

def build_formula(
     metal: str,
     ligands: dict[str, int],
     oxidation_state: int,
     total_ligand_charge: int,
     counter_ions: dict[str, int] | None = None,
) -> str:
     """Build a coordination complex formula from parsed name components."""
     ligand_parts = []
     for ligand, count in ligands.items():
          count_text = "" if count == 1 else str(count)
          if re.fullmatch(r"[A-Z][a-z]?", ligand):
               ligand_parts.append(f"{ligand}{count_text}")
          else:
               ligand_parts.append(f"({ligand}){count_text}")

     complex_charge = oxidation_state + total_ligand_charge
     if complex_charge == 0:
          charge_text = ""
     elif complex_charge == 1:
          charge_text = "+"
     elif complex_charge == -1:
          charge_text = "-"
     elif complex_charge > 1:
          charge_text = f"{complex_charge}+"
     else:
          charge_text = f"{abs(complex_charge)}-"

     inner = f"[{metal}{''.join(ligand_parts)}]"
     if counter_ions:
          leading = []
          trailing = []
          for ion, count in counter_ions.items():
               ion_text = _format_counter_ion_formula(ion, count)
               if COUNTER_IONS.get(ion, 0) > 0:
                    leading.append(ion_text)
               else:
                    trailing.append(ion_text)
          return f"{''.join(leading)}{inner}{''.join(trailing)}"

     return f"{inner}{charge_text}"


def _format_counter_ion_formula(ion: str, count: int) -> str:
     count_text = "" if count == 1 else str(count)
     if count == 1:
          return ion
     if re.fullmatch(r"[A-Z][a-z]?", ion):
          return f"{ion}{count_text}"
     return f"({ion}){count_text}"


def _extract_counter_ions_from_name(name: str) -> tuple[str, dict[str, int | None]]:
     """Remove known counter-ion names and return possible explicit counts."""
     normalized_name = _normalize_name(name)
     counter_ions: dict[str, int | None] = {}
     counter_names = sorted(
          COUNTER_ION_NAMES.items(),
          key=lambda item: len(_normalize_name(item[1])),
          reverse=True,
     )
     prefixes = sorted(PREFIXE.items(), key=lambda item: len(item[0]), reverse=True)

     for ion, ion_name in counter_names:
          normalized_ion_name = _normalize_name(ion_name)
          found = False
          for prefix, count in prefixes:
               pattern = prefix + normalized_ion_name
               if pattern in normalized_name:
                    counter_ions[ion] = count
                    normalized_name = normalized_name.replace(pattern, "", 1)
                    found = True
                    break
          if found:
               continue
          if normalized_ion_name in normalized_name:
               counter_ions[ion] = None
               normalized_name = normalized_name.replace(normalized_ion_name, "", 1)

     return normalized_name, counter_ions


def _resolve_counter_ion_counts(
     counter_ions: dict[str, int | None],
     complex_charge: int,
) -> dict[str, int]:
     if not counter_ions:
          return {}

     if len(counter_ions) == 1:
          ion, explicit_count = next(iter(counter_ions.items()))
          if explicit_count is not None:
               return {ion: explicit_count}
          ion_charge = COUNTER_IONS.get(ion, 0)
          if ion_charge and complex_charge and ion_charge * complex_charge < 0:
               count = abs(complex_charge) // abs(ion_charge)
               if count > 0 and count * abs(ion_charge) == abs(complex_charge):
                    return {ion: count}
          return {ion: 1}

     return {ion: count if count is not None else 1 for ion, count in counter_ions.items()}


def parse_name(name: str) -> ParsedComplex:
     raw_name= name
     name_without_counter_ions, counter_ion_names = _extract_counter_ions_from_name(name)
     metal= metal_data(name_without_counter_ions)
     oxydation_state= extract_complex_charge_from_name(name)
     ligands= ligand_data(name_without_counter_ions)

     result= ParsedComplex(metal=metal, ligands=ligands, complex_charge=0, counter_ions={},raw_formula=raw_name)
     _enrich(result)

     if oxydation_state is not None:
          result.oxidation_state=oxydation_state
          result.complex_charge=oxydation_state + result.total_ligand_charge
          result.counter_ions = _resolve_counter_ion_counts(
               counter_ion_names,
               result.complex_charge,
          )
          result.raw_formula = build_formula(
               metal,
               ligands,
               oxydation_state,
               result.total_ligand_charge,
               result.counter_ions,
          )
          _apply_ambidentate_donor_assignments(result)
          result.iupac_name = get_iupac_name(result)
     return result


def _normalize_name(name: str) -> str:
     """Normalize a coordination compound name for simple substring matching."""
     return re.sub(r"[\s_\-()]+", "", name).lower()

#test
if __name__=="__main__":
     name = input ("complex name is:")
     complex_info=parse_name(name)
     print(complex_info)
