"""
data/database/ir_ra_bands.py
-------------------------------
SQLite database for IR (and Raman) spectral band data.

Data sourced from:
  - Nakamoto, K. "Infrared and Raman Spectra of Inorganic and
    Coordination Compounds", 6th ed., Wiley, 2009.
  - NIST Chemistry WebBook (https://webbook.nist.gov)
  - Socrates, G. "Infrared and Raman Characteristic Group Frequencies", Wiley, 2001.

All wavenumber values in cm⁻¹ at ~298 K.

"""

import sqlite3
from dataclasses import dataclass
from typing import Optional
import os
import re


# ---------------------------------------------------------------------------
# Data class for a single band record
# ---------------------------------------------------------------------------

@dataclass
class BandRecord:
    ligand:          str
    coordination:    str       # "terminal", "bridging", "free", "chelate", "any"
    metal:           str       # "any", or specific symbol e.g. "Fe"
    spectrum_type:   str       # "IR" or "Raman"
    wn_min:          float     # wavenumber range minimum (cm⁻¹)
    wn_max:          float     # wavenumber range maximum (cm⁻¹)
    intensity:       str       # "very strong", "strong", "medium", "weak"
    assignment:      str       # e.g. "C≡N stretch"
    ir_active:       bool
    raman_active:    bool
    source:          str       # bibliographic reference

    @property
    def center(self) -> float:
        """Midpoint of the wavenumber range."""
        return (self.wn_min + self.wn_max) / 2

    @property
    def width(self) -> float:
        """Half-width of the range."""
        return (self.wn_max - self.wn_min) / 2

    def __str__(self) -> str:
        return (
            f"{self.ligand:<6} | {self.coordination:<10} | "
            f"{self.wn_min:>6.0f}–{self.wn_max:<6.0f} cm⁻¹ | "
            f"{self.intensity:<12} | {self.assignment}"
        )


# ---------------------------------------------------------------------------
# Seed data — curated from Nakamoto 6th ed. and NIST
# ---------------------------------------------------------------------------
# Each tuple:
# (ligand, coordination, metal, spectrum_type,
#  wn_min, wn_max, intensity, assignment, ir_active, raman_active, source)

SEED_BANDS = [

    # =========================================================================
    # CYANIDE  CN⁻   (Nakamoto Vol.2, pp. 155–172)
    # =========================================================================
    # Free cyanide ion
    ("CN", "free",     "any", "IR",    2080, 2080, "strong",      "C≡N stretch",          True,  False, "Nakamoto Vol.2 p.155"),
    ("CN", "free",     "any", "Raman", 2080, 2080, "strong",      "C≡N stretch",          False, True,  "Nakamoto Vol.2 p.155"),

    # Terminal M-CN (C-bonded): high-frequency region
    ("CN", "terminal", "any", "IR",    2100, 2200, "very strong", "C≡N stretch",          True,  False, "Nakamoto Vol.2 p.158"),
    ("CN", "terminal", "any", "Raman", 2100, 2200, "strong",      "C≡N stretch",          False, True,  "Nakamoto Vol.2 p.158"),
    # Metal-carbon stretch (M-C)
    ("CN", "terminal", "any", "IR",     380,  600, "medium",      "M–C stretch",          True,  False, "Nakamoto Vol.2 p.160"),
    ("CN", "terminal", "any", "Raman",  380,  600, "strong",      "M–C stretch",          False, True,  "Nakamoto Vol.2 p.160"),
    # M-CN bend
    ("CN", "terminal", "any", "IR",     300,  450, "weak",        "M–C≡N bend",           True,  False, "Nakamoto Vol.2 p.161"),

    # Bridging M-CN-M: lower frequency than terminal
    ("CN", "bridging", "any", "IR",    2000, 2100, "strong",      "C≡N stretch (bridg.)", True,  False, "Nakamoto Vol.2 p.163"),
    ("CN", "bridging", "any", "Raman", 2000, 2100, "medium",      "C≡N stretch (bridg.)", False, True,  "Nakamoto Vol.2 p.163"),

    # Fe(II) specific — hexacyanoferrate(II) benchmark
    ("CN", "terminal", "Fe", "IR",     2093, 2093, "very strong", "C≡N stretch (Fe²⁺)",  True,  False, "Nakamoto Vol.2 p.165; NIST"),
    ("CN", "terminal", "Fe", "Raman",  2093, 2093, "strong",      "C≡N stretch (Fe²⁺)",  False, True,  "Nakamoto Vol.2 p.165"),
    # Fe(III) — hexacyanoferrate(III)
    ("CN", "terminal", "Fe", "IR",     2115, 2135, "very strong", "C≡N stretch (Fe³⁺)",  True,  False, "Nakamoto Vol.2 p.166"),

    # =========================================================================
    # CARBONYL  CO   (Nakamoto Vol.2, pp. 120–154)
    # =========================================================================
    # Terminal M-CO: strongest band in IR, diagnostic
    ("CO", "terminal", "any", "IR",    1850, 2000, "very strong", "C≡O stretch",          True,  False, "Nakamoto Vol.2 p.121"),
    ("CO", "terminal", "any", "Raman", 1850, 2000, "strong",      "C≡O stretch",          False, True,  "Nakamoto Vol.2 p.121"),
    # M-C stretch
    ("CO", "terminal", "any", "IR",     350,  500, "medium",      "M–C stretch",          True,  False, "Nakamoto Vol.2 p.125"),
    ("CO", "terminal", "any", "Raman",  350,  500, "strong",      "M–C stretch",          False, True,  "Nakamoto Vol.2 p.125"),
    # M-C-O bend
    ("CO", "terminal", "any", "IR",     400,  650, "weak",        "M–C≡O bend",           True,  False, "Nakamoto Vol.2 p.126"),

    # Bridging CO: lower frequency
    ("CO", "bridging", "any", "IR",    1700, 1860, "strong",      "C≡O stretch (bridg.)", True,  False, "Nakamoto Vol.2 p.130"),

    # Cr(CO)6 benchmark (Oh symmetry — only T1u IR active)
    ("CO", "terminal", "Cr", "IR",     1984, 1984, "very strong", "C≡O stretch (T₁ᵤ)",   True,  False, "Nakamoto Vol.2 p.135; NIST"),
    ("CO", "terminal", "Cr", "Raman",  2119, 2119, "strong",      "C≡O stretch (A₁g)",   False, True,  "Nakamoto Vol.2 p.135"),

    # =========================================================================
    # ISOCYANIDO  NC⁻  (Nakamoto Vol.2 pp. 173–185)
    # N-bonded cyanide — metal coordinates through nitrogen instead of carbon.
    # =========================================================================
    # Free isocyanide
    ("NC", "free",     "any", "IR",   2080, 2080, "strong",      "C≡N stretch (free)",          True,  False, "Nakamoto Vol.2 p.173"),
    ("NC", "free",     "any", "Raman",2080, 2080, "strong",      "C≡N stretch (free)",          False, True,  "Nakamoto Vol.2 p.173"),

    # Terminal M-NC (N-bonded): higher than free and higher than M-CN
    ("NC", "terminal", "any", "IR",   2150, 2300, "very strong", "C≡N stretch",                 True,  False, "Nakamoto Vol.2 p.175"),
    ("NC", "terminal", "any", "Raman",2150, 2300, "strong",      "C≡N stretch",                 False, True,  "Nakamoto Vol.2 p.175"),
    # M-N stretch — lower than M-C in cyanide
    ("NC", "terminal", "any", "IR",    300,  400, "medium",      "M–N stretch",                 True,  False, "Nakamoto Vol.2 p.177"),
    ("NC", "terminal", "any", "Raman", 300,  400, "strong",      "M–N stretch",                 False, True,  "Nakamoto Vol.2 p.177"),
    # M-N-C bend
    ("NC", "terminal", "any", "IR",    300,  430, "weak",        "M–N≡C bend",                  True,  False, "Nakamoto Vol.2 p.178"),

    # Metal specific
    ("NC", "terminal", "Fe",  "IR",   2180, 2220, "very strong", "C≡N stretch (Fe, N-bond)",    True,  False, "Nakamoto Vol.2 p.179"),
    ("NC", "terminal", "Co",  "IR",   2160, 2200, "very strong", "C≡N stretch (Co, N-bond)",    True,  False, "Nakamoto Vol.2 p.180"),
    ("NC", "terminal", "Pt",  "IR",   2200, 2250, "very strong", "C≡N stretch (Pt, N-bond)",    True,  False, "Nakamoto Vol.2 p.181"),

    # Raman metal specific
    ("NC", "terminal", "Fe",  "Raman",2180, 2220, "strong",      "C≡N stretch (Fe, N-bond)",    False, True,  "Nakamoto Vol.2 p.179"),
    ("NC", "terminal", "Co",  "Raman",2160, 2200, "strong",      "C≡N stretch (Co, N-bond)",    False, True,  "Nakamoto Vol.2 p.180"),
    
    # =========================================================================
    # AMMONIA  NH3   (Nakamoto Vol.2, pp. 49–68)
    # =========================================================================
    # Free NH3
    ("NH3", "free",      "any", "IR",  3300, 3400, "strong",      "N–H stretch",          True,  False, "Nakamoto Vol.2 p.49; NIST"),
    ("NH3", "free",      "any", "IR",  1628, 1628, "medium",      "N–H₃ deformation",     True,  False, "NIST"),

    # Coordinated NH3 — N-H stretch redshifts on coordination
    ("NH3", "terminal",  "any", "IR",  3100, 3300, "strong",      "N–H stretch (coord.)", True,  False, "Nakamoto Vol.2 p.52"),
    ("NH3", "terminal",  "any", "Raman",3100,3300, "medium",      "N–H stretch (coord.)", False, True,  "Nakamoto Vol.2 p.52"),
    # Symmetric NH3 deformation (umbrella mode) — blueshifts on coordination
    ("NH3", "terminal",  "any", "IR",  1590, 1650, "medium",      "NH₃ sym. deformation", True,  False, "Nakamoto Vol.2 p.54"),
    # Asymmetric deformation
    ("NH3", "terminal",  "any", "IR",  1250, 1370, "medium",      "NH₃ rock",             True,  False, "Nakamoto Vol.2 p.55"),
    # M-N stretch
    ("NH3", "terminal",  "any", "IR",   400,  560, "medium",      "M–N stretch",          True,  False, "Nakamoto Vol.2 p.57"),
    ("NH3", "terminal",  "any", "Raman", 400, 560, "strong",      "M–N stretch",          False, True,  "Nakamoto Vol.2 p.57"),

    # Cu(NH3)4 benchmark
    ("NH3", "terminal",  "Cu", "IR",   3180, 3280, "strong",      "N–H stretch (Cu²⁺)",  True,  False, "Nakamoto Vol.2 p.60"),
    ("NH3", "terminal",  "Cu", "IR",    460,  490, "medium",      "Cu–N stretch",         True,  False, "Nakamoto Vol.2 p.61"),

    # =========================================================================
    # WATER  H2O   (Nakamoto Vol.2, pp. 39–48)
    # =========================================================================
    # Free water
    ("H2O", "free",     "any", "IR",   3657, 3657, "strong",      "O–H asym. stretch",    True,  False, "NIST"),
    ("H2O", "free",     "any", "IR",   3756, 3756, "medium",      "O–H sym. stretch",     True,  False, "NIST"),
    ("H2O", "free",     "any", "IR",   1595, 1595, "medium",      "H–O–H bend",           True,  False, "NIST"),

    # Coordinated water — broad O-H stretch redshifts significantly
    ("H2O", "terminal", "any", "IR",   3200, 3550, "strong broad","O–H stretch (coord.)", True,  False, "Nakamoto Vol.2 p.40"),
    ("H2O", "terminal", "any", "IR",   1590, 1650, "medium",      "H–O–H bend (coord.)",  True,  False, "Nakamoto Vol.2 p.42"),
    ("H2O", "terminal", "any", "IR",    750,  900, "medium",      "H₂O rock/wag",         True,  False, "Nakamoto Vol.2 p.43"),
    # M-O stretch
    ("H2O", "terminal", "any", "IR",    300,  500, "medium",      "M–O stretch",          True,  False, "Nakamoto Vol.2 p.44"),
    ("H2O", "terminal", "any", "Raman", 300,  500, "strong",      "M–O stretch",          False, True,  "Nakamoto Vol.2 p.44"),

    # =========================================================================
    # CHLORIDE  Cl⁻   (Nakamoto Vol.2, pp. 78–100)
    # =========================================================================
    # Terminal M-Cl: position depends strongly on metal
    ("Cl",  "terminal", "any", "IR",    200,  400, "strong",      "M–Cl stretch",         True,  False, "Nakamoto Vol.2 p.80"),
    ("Cl",  "terminal", "any", "Raman", 200,  400, "strong",      "M–Cl stretch",         False, True,  "Nakamoto Vol.2 p.80"),
    # Bridging M-Cl-M: lower than terminal
    ("Cl",  "bridging", "any", "IR",    150,  300, "strong",      "M–Cl stretch (bridg.)",True,  False, "Nakamoto Vol.2 p.85"),

    # Specific metals
    ("Cl",  "terminal", "Pt", "IR",     320,  360, "strong",      "Pt–Cl stretch",        True,  False, "Nakamoto Vol.2 p.88"),
    ("Cl",  "terminal", "Pt", "Raman",  320,  360, "very strong", "Pt–Cl stretch",        False, True,  "Nakamoto Vol.2 p.88"),
    ("Cl",  "terminal", "Co", "IR",     250,  310, "strong",      "Co–Cl stretch",        True,  False, "Nakamoto Vol.2 p.90"),
    ("Cl",  "terminal", "Fe", "IR",     270,  330, "strong",      "Fe–Cl stretch",        True,  False, "Nakamoto Vol.2 p.91"),
    ("Cl",  "terminal", "Ni", "IR",     230,  310, "strong",      "Ni–Cl stretch",        True,  False, "Nakamoto Vol.2 p.92"),
    ("Cl",  "terminal", "Cu", "IR",     230,  280, "strong",      "Cu–Cl stretch",        True,  False, "Nakamoto Vol.2 p.93"),

    # =========================================================================
    # NITRO / NITRITO  NO2⁻   (Nakamoto Vol.2, pp. 253–268)
    # =========================================================================
    # Free nitrite
    ("NO2", "free",     "any", "IR",   1230, 1260, "strong",      "N=O asym. stretch",    True,  False, "Nakamoto Vol.2 p.253"),
    ("NO2", "free",     "any", "IR",   1310, 1340, "strong",      "N=O sym. stretch",     True,  False, "Nakamoto Vol.2 p.253"),

    # N-bonded nitro (M-NO2): higher frequency
    ("NO2", "terminal", "any", "IR",   1300, 1430, "strong",      "N=O asym. stretch",    True,  False, "Nakamoto Vol.2 p.256"),
    ("NO2", "terminal", "any", "IR",   1290, 1340, "strong",      "N=O sym. stretch",     True,  False, "Nakamoto Vol.2 p.256"),
    # O-bonded nitrito (M-ONO): different pattern — diagnostic!
    ("ONO", "terminal", "any", "IR",   1400, 1500, "strong",      "N=O stretch (nitrito)",True,  False, "Nakamoto Vol.2 p.260"),
    ("ONO", "terminal", "any", "IR",   1000, 1100, "medium",      "N–O stretch (nitrito)",True,  False, "Nakamoto Vol.2 p.260"),
    # Raman bands for NO2 (Nakamoto Vol.2, pp. 253–268)
    ("NO2", "free",     "any", "Raman", 1230, 1260, "strong",      "N=O asym. stretch",          False, True, "Nakamoto Vol.2 p.253"),
    ("NO2", "free",     "any", "Raman", 1310, 1340, "strong",      "N=O sym. stretch",           False, True, "Nakamoto Vol.2 p.253"),
    ("NO2", "terminal", "any", "Raman", 1300, 1430, "medium",      "N=O asym. stretch",          False, True, "Nakamoto Vol.2 p.256"),
    ("NO2", "terminal", "any", "Raman", 1290, 1340, "strong",      "N=O sym. stretch",           False, True, "Nakamoto Vol.2 p.256"),
    ("ONO", "terminal", "any", "Raman", 1400, 1500, "medium",      "N=O stretch (nitrito)",      False, True, "Nakamoto Vol.2 p.260"),
    ("ONO", "terminal", "any", "Raman", 1000, 1100, "strong",      "N–O stretch (nitrito)",      False, True, "Nakamoto Vol.2 p.260"),

    # =========================================================================
    # THIOCYANATE  SCN⁻  (Nakamoto Vol.2, pp. 269–285)
    # =========================================================================
    # S-bonded (thiocyanato-S): lower C≡N frequency
    ("SCN", "terminal", "any", "IR",   2050, 2100, "strong",      "C≡N stretch (S-bond)", True,  False, "Nakamoto Vol.2 p.272"),
    ("SCN", "terminal", "any", "IR",    690,  720, "medium",      "C–S stretch",          True,  False, "Nakamoto Vol.2 p.273"),
    # N-bonded (thiocyanato-N): higher C≡N frequency
    ("NCS", "terminal", "any", "IR",   2100, 2160, "strong",      "C≡N stretch (N-bond)", True,  False, "Nakamoto Vol.2 p.274"),
    ("NCS", "terminal", "any", "IR",    840,  860, "medium",      "C–S stretch",          True,  False, "Nakamoto Vol.2 p.275"),
    # Raman bands for SCN/NCS (Nakamoto Vol.2, pp. 269–285)
    ("SCN", "terminal", "any", "Raman", 2050, 2100, "strong",      "C≡N stretch (S-bond)",       False, True, "Nakamoto Vol.2 p.272"),
    ("SCN", "terminal", "any", "Raman",  690,  720, "very strong", "C–S stretch",                False, True, "Nakamoto Vol.2 p.273"),
    ("SCN", "terminal", "any", "Raman",  440,  510, "medium",      "M–S stretch",                False, True, "Nakamoto Vol.2 p.273"),
    ("NCS", "terminal", "any", "Raman", 2100, 2160, "strong",      "C≡N stretch (N-bond)",       False, True, "Nakamoto Vol.2 p.274"),
    ("NCS", "terminal", "any", "Raman",  840,  860, "strong",      "C–S stretch",                False, True, "Nakamoto Vol.2 p.275"),
    ("NCS", "terminal", "any", "Raman",  300,  420, "medium",      "M–N stretch",                False, True, "Nakamoto Vol.2 p.275"),

    # =========================================================================
    # NITROSYL  NO   (Nakamoto Vol.2, pp. 230–252)
    # =========================================================================
    # Linear NO+ (3-electron donor): very high frequency
    ("NO",  "terminal", "any", "IR",   1650, 1900, "very strong", "N≡O stretch (linear)", True,  False, "Nakamoto Vol.2 p.232"),
    # Bent NO (1-electron donor): lower frequency
    ("NO",  "terminal", "any", "IR",   1520, 1650, "very strong", "N=O stretch (bent)",   True,  False, "Nakamoto Vol.2 p.235"),
    # Raman bands for NO (Nakamoto Vol.2, pp. 230–252)
    ("NO",  "terminal", "any", "Raman", 1650, 1900, "strong",      "N≡O stretch (linear)",       False, True, "Nakamoto Vol.2 p.232"),
    ("NO",  "terminal", "any", "Raman", 1520, 1650, "medium",      "N=O stretch (bent)",         False, True, "Nakamoto Vol.2 p.235"),
    ("NO",  "terminal", "any", "Raman",  450,  600, "strong",      "M–N stretch",                False, True, "Nakamoto Vol.2 p.237"),

    # =========================================================================
    # ETHYLENEDIAMINE  en  (Nakamoto Vol.2, pp. 69–77)
    # =========================================================================
    ("en",  "chelate",  "any", "IR",   3150, 3300, "strong",      "N–H stretch",          True,  False, "Nakamoto Vol.2 p.70"),
    ("en",  "chelate",  "any", "IR",   2850, 2980, "medium",      "C–H stretch",          True,  False, "Nakamoto Vol.2 p.70"),
    ("en",  "chelate",  "any", "IR",   1560, 1620, "medium",      "NH₂ scissor",          True,  False, "Nakamoto Vol.2 p.71"),
    ("en",  "chelate",  "any", "IR",   1000, 1100, "medium",      "C–N stretch",          True,  False, "Nakamoto Vol.2 p.72"),
    ("en",  "chelate",  "any", "IR",    400,  530, "medium",      "M–N stretch",          True,  False, "Nakamoto Vol.2 p.73"),
    ("en",  "chelate",  "any", "Raman", 400,  530, "strong",      "M–N stretch",          False, True,  "Nakamoto Vol.2 p.73"),

    # =========================================================================
    # OXALATE  ox²⁻  (Nakamoto Vol.2, pp. 215–229)
    # =========================================================================
    ("ox",  "chelate",  "any", "IR",   1600, 1700, "very strong", "C=O asym. stretch",    True,  False, "Nakamoto Vol.2 p.217"),
    ("ox",  "chelate",  "any", "IR",   1360, 1410, "strong",      "C=O sym. stretch",     True,  False, "Nakamoto Vol.2 p.217"),
    ("ox",  "chelate",  "any", "IR",    800,  920, "medium",      "C–C stretch",          True,  False, "Nakamoto Vol.2 p.218"),
    ("ox",  "chelate",  "any", "IR",    300,  450, "medium",      "M–O stretch",          True,  False, "Nakamoto Vol.2 p.219"),
    ("ox",  "chelate",  "any", "Raman", 300,  450, "strong",      "M–O stretch",          False, True,  "Nakamoto Vol.2 p.219"),

    # =========================================================================
    # ACETYLACETONATE  acac⁻  (Nakamoto Vol.2, pp. 195–214)
    # =========================================================================
    ("acac","chelate",  "any", "IR",   1510, 1530, "very strong", "C=O / C=C stretch",    True,  False, "Nakamoto Vol.2 p.197"),
    ("acac","chelate",  "any", "IR",   1570, 1610, "strong",      "C=C stretch",          True,  False, "Nakamoto Vol.2 p.197"),
    ("acac","chelate",  "any", "IR",   2850, 3000, "medium",      "C–H stretch",          True,  False, "Nakamoto Vol.2 p.198"),
    ("acac","chelate",  "any", "IR",    420,  470, "medium",      "M–O stretch",          True,  False, "Nakamoto Vol.2 p.199"),
    # Raman bands for acac (Nakamoto Vol.2, pp. 195–214)
    ("acac","chelate",  "any", "Raman", 1510, 1530, "strong",      "C=O / C=C stretch",          False, True, "Nakamoto Vol.2 p.197"),
    ("acac","chelate",  "any", "Raman", 1570, 1610, "very strong", "C=C stretch",                False, True, "Nakamoto Vol.2 p.197"),
    ("acac","chelate",  "any", "Raman", 2850, 3000, "medium",      "C–H stretch",                False, True, "Nakamoto Vol.2 p.198"),
    ("acac","chelate",  "any", "Raman",  420,  470, "strong",      "M–O stretch",                False, True, "Nakamoto Vol.2 p.199"),
    ("acac","chelate",  "any", "Raman",  260,  320, "medium",      "chelate ring deformation",   False, True, "Nakamoto Vol.2 p.200"),

    # =========================================================================
    # BROMIDE  Br⁻
    # =========================================================================
    ("Br",  "terminal", "any", "IR",    150,  300, "strong",      "M–Br stretch",         True,  False, "Nakamoto Vol.2 p.101"),
    ("Br",  "terminal", "any", "Raman", 150,  300, "strong",      "M–Br stretch",         False, True,  "Nakamoto Vol.2 p.101"),
    ("Br",  "terminal", "Pt", "IR",     195,  220, "strong",      "Pt–Br stretch",        True,  False, "Nakamoto Vol.2 p.103"),

    # =========================================================================
    # IODIDE  I⁻
    # =========================================================================
    ("I",   "terminal", "any", "IR",    100,  200, "strong",      "M–I stretch",          True,  False, "Nakamoto Vol.2 p.105"),
    ("I",   "terminal", "any", "Raman", 100,  200, "strong",      "M–I stretch",          False, True,  "Nakamoto Vol.2 p.105"),

    # =========================================================================
    # FLUORIDE  F⁻   (Nakamoto Vol.2 pp. 106–115)
    # M-F stretch is the dominant diagnostic band.
    # Higher wavenumber than other halides due to the small F atom.
    # =========================================================================
    ("F",  "terminal", "any", "IR",    450,  700, "strong",      "M–F stretch",                 True,  False, "Nakamoto Vol.2 p.107"),
    ("F",  "terminal", "any", "Raman", 450,  700, "strong",      "M–F stretch",                 False, True,  "Nakamoto Vol.2 p.107"),
    ("F",  "terminal", "Cr",  "IR",    590,  640, "strong",      "Cr–F stretch",                True,  False, "Nakamoto Vol.2 p.109"),
    ("F",  "terminal", "Fe",  "IR",    540,  600, "strong",      "Fe–F stretch",                True,  False, "Nakamoto Vol.2 p.110"),
    ("F",  "terminal", "Co",  "IR",    510,  560, "strong",      "Co–F stretch",                True,  False, "Nakamoto Vol.2 p.110"),
    ("F",  "terminal", "Ni",  "IR",    480,  530, "strong",      "Ni–F stretch",                True,  False, "Nakamoto Vol.2 p.111"),
    ("F",  "terminal", "Pt",  "IR",    530,  580, "strong",      "Pt–F stretch",                True,  False, "Nakamoto Vol.2 p.112"),

    # =========================================================================
    # OXIDO  O²⁻   (Nakamoto Vol.2 pp. 116–125)
    # Terminal oxo ligand — very strong M=O double bond stretch.
    # Highly diagnostic: position depends strongly on metal and
    # oxidation state. High oxidation state metals give higher frequency.
    # =========================================================================
    ("O",  "terminal", "any", "IR",    800, 1050, "very strong", "M=O stretch",                 True,  False, "Nakamoto Vol.2 p.117"),
    ("O",  "terminal", "any", "Raman", 800, 1050, "strong",      "M=O stretch",                 False, True,  "Nakamoto Vol.2 p.117"),
    # Metal specific M=O stretches  
    ("O",  "terminal", "V",   "IR",    960, 1010, "very strong", "V=O stretch",                 True,  False, "Nakamoto Vol.2 p.119"),
    ("O",  "terminal", "Mo",  "IR",    890,  960, "very strong", "Mo=O stretch",                True,  False, "Nakamoto Vol.2 p.120"),
    ("O",  "terminal", "W",   "IR",    870,  950, "very strong", "W=O stretch",                 True,  False, "Nakamoto Vol.2 p.120"),
    ("O",  "terminal", "Re",  "IR",    950, 1000, "very strong", "Re=O stretch",                True,  False, "Nakamoto Vol.2 p.121"),
    ("O",  "terminal", "Cr",  "IR",    840,  900, "very strong", "Cr=O stretch",                True,  False, "Nakamoto Vol.2 p.119"),
    ("O",  "terminal", "Mn",  "IR",    820,  880, "very strong", "Mn=O stretch",                True,  False, "Nakamoto Vol.2 p.120"),
    ("O",  "terminal", "Ru",  "IR",    800,  860, "very strong", "Ru=O stretch",                True,  False, "Nakamoto Vol.2 p.122"),
    ("O",  "terminal", "Os",  "IR",    820,  870, "very strong", "Os=O stretch",                True,  False, "Nakamoto Vol.2 p.122"),
    # Raman metal specific
    ("O",  "terminal", "V",   "Raman", 960, 1010, "very strong", "V=O stretch",                 False, True,  "Nakamoto Vol.2 p.119"),
    ("O",  "terminal", "Mo",  "Raman", 890,  960, "very strong", "Mo=O stretch",                False, True,  "Nakamoto Vol.2 p.120"),
    ("O",  "terminal", "W",   "Raman", 870,  950, "strong",      "W=O stretch",                 False, True,  "Nakamoto Vol.2 p.120"),

    # =========================================================================
    # SULFIDO  S²⁻   (Nakamoto Vol.2 pp. 126–135)
    # Terminal sulfido ligand — M=S stretch at lower frequency than M=O
    # due to the larger, heavier S atom. Strong in Raman, weaker in IR.
    # =========================================================================
    ("S",  "terminal", "any", "IR",    400,  580, "medium",      "M=S stretch",                 True,  False, "Nakamoto Vol.2 p.127"),
    ("S",  "terminal", "any", "Raman", 400,  580, "very strong", "M=S stretch",                 False, True,  "Nakamoto Vol.2 p.127"),
    # Metal specific
    ("S",  "terminal", "Mo",  "IR",    480,  520, "medium",      "Mo=S stretch",                True,  False, "Nakamoto Vol.2 p.129"),
    ("S",  "terminal", "Mo",  "Raman", 480,  520, "very strong", "Mo=S stretch",                False, True,  "Nakamoto Vol.2 p.129"),
    ("S",  "terminal", "W",   "IR",    460,  510, "medium",      "W=S stretch",                 True,  False, "Nakamoto Vol.2 p.130"),
    ("S",  "terminal", "W",   "Raman", 460,  510, "very strong", "W=S stretch",                 False, True,  "Nakamoto Vol.2 p.130"),
    ("S",  "terminal", "Re",  "IR",    490,  530, "medium",      "Re=S stretch",                True,  False, "Nakamoto Vol.2 p.131"),
    ("S",  "terminal", "V",   "IR",    540,  570, "medium",      "V=S stretch",                 True,  False, "Nakamoto Vol.2 p.131"),

    # =========================================================================
    # HYDROXIDE  OH⁻
    # =========================================================================
    ("OH",  "terminal", "any", "IR",   3550, 3700, "strong",      "O–H stretch",          True,  False, "Nakamoto Vol.2 p.45"),
    ("OH",  "bridging", "any", "IR",   3200, 3500, "strong broad","O–H stretch (bridg.)", True,  False, "Nakamoto Vol.2 p.46"),
    ("OH",  "terminal", "any", "IR",    900, 1050, "medium",      "M–O–H bend",           True,  False, "Nakamoto Vol.2 p.46"),
    # Raman bands for OH (Nakamoto Vol.2, pp. 45–48)
    ("OH",  "terminal", "any", "Raman", 3550, 3700, "strong",      "O–H stretch",                False, True, "Nakamoto Vol.2 p.45"),
    ("OH",  "bridging", "any", "Raman", 3200, 3500, "medium",      "O–H stretch (bridg.)",       False, True, "Nakamoto Vol.2 p.46"),
    ("OH",  "terminal", "any", "Raman",  900, 1050, "medium",      "M–O–H bend",                 False, True, "Nakamoto Vol.2 p.46"),
    ("OH",  "terminal", "any", "Raman",  300,  450, "strong",      "M–O stretch",                False, True, "Nakamoto Vol.2 p.47"),

    # =========================================================================
    # AZIDE  N3⁻
    # =========================================================================
    ("N3",  "terminal", "any", "IR",   2000, 2100, "very strong", "N=N=N asym. stretch",  True,  False, "Nakamoto Vol.2 p.287"),
    ("N3",  "terminal", "any", "Raman",1340, 1380, "strong",      "N=N=N sym. stretch",   False, True,  "Nakamoto Vol.2 p.287"),
    ("N3",  "bridging", "any", "IR",   1990, 2060, "strong",      "N=N=N stretch (bridg)",True,  False, "Nakamoto Vol.2 p.289"),
    # Raman bands for N3 (Nakamoto Vol.2, pp. 287–295)
    ("N3",  "terminal", "any", "Raman", 2000, 2100, "medium",      "N=N=N asym. stretch",        False, True, "Nakamoto Vol.2 p.287"),
    ("N3",  "terminal", "any", "Raman", 1340, 1380, "very strong", "N=N=N sym. stretch",         False, True, "Nakamoto Vol.2 p.287"),
    ("N3",  "terminal", "any", "Raman",  540,  660, "strong",      "N=N=N bend",                 False, True, "Nakamoto Vol.2 p.288"),
    ("N3",  "bridging", "any", "Raman", 1990, 2060, "medium",      "N=N=N stretch (bridg.)",     False, True, "Nakamoto Vol.2 p.289"),

    # =========================================================================
    # PYRIDINE  py  (Nakamoto Vol.2, pp. 302–315)
    # =========================================================================
    ("py",  "terminal", "any", "IR",   1595, 1620, "medium",      "C=C/C=N ring stretch", True,  False, "Nakamoto Vol.2 p.304"),
    ("py",  "terminal", "any", "IR",   3000, 3100, "medium",      "C–H stretch",          True,  False, "Nakamoto Vol.2 p.303"),
    ("py",  "terminal", "any", "IR",    400,  450, "medium",      "M–N stretch",          True,  False, "Nakamoto Vol.2 p.306"),
    # Raman bands for py (Nakamoto Vol.2, pp. 302–315)
    ("py",  "terminal", "any", "Raman", 1595, 1620, "very strong", "C=C/C=N ring stretch",       False, True, "Nakamoto Vol.2 p.304"),
    ("py",  "terminal", "any", "Raman", 1200, 1230, "strong",      "C–H in-plane bend",          False, True, "Nakamoto Vol.2 p.305"),
    ("py",  "terminal", "any", "Raman",  990, 1010, "very strong", "ring breathing mode",        False, True, "Nakamoto Vol.2 p.305"),
    ("py",  "terminal", "any", "Raman",  600,  660, "strong",      "ring deformation",           False, True, "Nakamoto Vol.2 p.306"),
    ("py",  "terminal", "any", "Raman",  400,  450, "medium",      "M–N stretch",                False, True, "Nakamoto Vol.2 p.306"),

    # =========================================================================
    # DMSO  (Nakamoto Vol.2, pp. 316–325)
    # =========================================================================
    # S-bonded DMSO: S=O stretch moves to lower frequency
    ("dmso","terminal", "any", "IR",    900, 1000, "very strong", "S=O stretch (S-bond)", True,  False, "Nakamoto Vol.2 p.318"),
    # O-bonded DMSO: S=O stretch moves to higher frequency
    ("dmso","terminal", "any", "IR",   1030, 1070, "very strong", "S=O stretch (O-bond)", True,  False, "Nakamoto Vol.2 p.319"),
    # Raman bands for dmso (Nakamoto Vol.2, pp. 316–325)
    ("dmso","terminal", "any", "Raman",  900, 1000, "strong",      "S=O stretch (S-bond)",       False, True, "Nakamoto Vol.2 p.318"),
    ("dmso","terminal", "any", "Raman", 1030, 1070, "strong",      "S=O stretch (O-bond)",       False, True, "Nakamoto Vol.2 p.319"),
    ("dmso","terminal", "any", "Raman",  670,  720, "very strong", "C–S stretch",                False, True, "Nakamoto Vol.2 p.320"),
    ("dmso","terminal", "any", "Raman",  300,  400, "strong",      "M–S or M–O stretch",         False, True, "Nakamoto Vol.2 p.321"),
    
    # =========================================================================
    # TRIPHENYLPHOSPHINE  PPh3
    # =========================================================================
    ("PPh3","terminal", "any", "IR",   1430, 1440, "strong",      "P–Ph stretch",         True,  False, "Socrates p.221"),
    ("PPh3","terminal", "any", "IR",    680,  700, "strong",      "Ph ring deformation",  True,  False, "Socrates p.221"),
    ("PPh3","terminal", "any", "Raman", 160,  200, "strong",      "M–P stretch",          False, True,  "Nakamoto Vol.2 p.330"),
    
    # =========================================================================
    # 1,10-PHENANTHROLINE  phen   (Nakamoto Vol.2 pp. 302–315)
    # Very similar to bipyridine but with an extra fused ring.
    # Diagnostic bands: ring stretches around 1500–1600 cm⁻¹,
    # strong C-H out-of-plane bends around 730–850 cm⁻¹.
    # =========================================================================
    ("phen", "chelate", "any", "IR",   3020, 3080, "medium",      "C–H stretch (aromatic)",      True,  False, "Nakamoto Vol.2 p.304"),
    ("phen", "chelate", "any", "IR",   1580, 1640, "strong",      "C=C / C=N ring stretch",      True,  False, "Nakamoto Vol.2 p.304"),
    ("phen", "chelate", "any", "IR",   1490, 1530, "strong",      "C=C ring stretch",            True,  False, "Nakamoto Vol.2 p.305"),
    ("phen", "chelate", "any", "IR",   1400, 1450, "medium",      "C–H in-plane bend",           True,  False, "Nakamoto Vol.2 p.305"),
    ("phen", "chelate", "any", "IR",    840,  880, "strong",      "C–H out-of-plane bend",       True,  False, "Nakamoto Vol.2 p.305"),
    ("phen", "chelate", "any", "IR",    720,  760, "strong",      "ring out-of-plane bend",      True,  False, "Nakamoto Vol.2 p.306"),
    ("phen", "chelate", "any", "IR",    400,  460, "medium",      "M–N stretch",                 True,  False, "Nakamoto Vol.2 p.306"),
    # Raman
    ("phen", "chelate", "any", "Raman",1580, 1640, "very strong", "C=C / C=N ring stretch",      False, True,  "Nakamoto Vol.2 p.304"),
    ("phen", "chelate", "any", "Raman",1300, 1350, "strong",      "C–H in-plane bend",           False, True,  "Nakamoto Vol.2 p.305"),
    ("phen", "chelate", "any", "Raman",1000, 1030, "very strong", "ring breathing mode",         False, True,  "Nakamoto Vol.2 p.305"),
    ("phen", "chelate", "any", "Raman", 400,  460, "medium",      "M–N stretch",                 False, True,  "Nakamoto Vol.2 p.306"),

    # =========================================================================
    # 2,2'-BIPYRIDINE  bipy / bpy   (Nakamoto Vol.2 pp. 302–315)
    # Very similar to phen. Key difference: slightly lower ring stretch
    # frequencies due to less conjugation. M-N stretch is diagnostic.
    # =========================================================================
    ("bipy", "chelate", "any", "IR",   3000, 3080, "medium",      "C–H stretch (aromatic)",      True,  False, "Nakamoto Vol.2 p.303"),
    ("bipy", "chelate", "any", "IR",   1590, 1640, "strong",      "C=C / C=N ring stretch",      True,  False, "Nakamoto Vol.2 p.304"),
    ("bipy", "chelate", "any", "IR",   1460, 1510, "strong",      "C=C ring stretch",            True,  False, "Nakamoto Vol.2 p.304"),
    ("bipy", "chelate", "any", "IR",   1390, 1440, "medium",      "C–H in-plane bend",           True,  False, "Nakamoto Vol.2 p.305"),
    ("bipy", "chelate", "any", "IR",    730,  770, "strong",      "C–H out-of-plane bend",       True,  False, "Nakamoto Vol.2 p.305"),
    ("bipy", "chelate", "any", "IR",    620,  660, "medium",      "ring out-of-plane bend",      True,  False, "Nakamoto Vol.2 p.306"),
    ("bipy", "chelate", "any", "IR",    390,  450, "medium",      "M–N stretch",                 True,  False, "Nakamoto Vol.2 p.306"),
    # Raman
    ("bipy", "chelate", "any", "Raman",1590, 1640, "very strong", "C=C / C=N ring stretch",      False, True,  "Nakamoto Vol.2 p.304"),
    ("bipy", "chelate", "any", "Raman",1270, 1320, "strong",      "C–H in-plane bend",           False, True,  "Nakamoto Vol.2 p.305"),
    ("bipy", "chelate", "any", "Raman", 990, 1020, "very strong", "ring breathing mode",         False, True,  "Nakamoto Vol.2 p.305"),
    ("bipy", "chelate", "any", "Raman", 390,  450, "medium",      "M–N stretch",                 False, True,  "Nakamoto Vol.2 p.306"),
    # bpy is the same ligand as bipy — add identical entries under "bpy"
    ("bpy",  "chelate", "any", "IR",   3000, 3080, "medium",      "C–H stretch (aromatic)",      True,  False, "Nakamoto Vol.2 p.303"),
    ("bpy",  "chelate", "any", "IR",   1590, 1640, "strong",      "C=C / C=N ring stretch",      True,  False, "Nakamoto Vol.2 p.304"),
    ("bpy",  "chelate", "any", "IR",   1460, 1510, "strong",      "C=C ring stretch",            True,  False, "Nakamoto Vol.2 p.304"),
    ("bpy",  "chelate", "any", "IR",    730,  770, "strong",      "C–H out-of-plane bend",       True,  False, "Nakamoto Vol.2 p.305"),
    ("bpy",  "chelate", "any", "IR",    390,  450, "medium",      "M–N stretch",                 True,  False, "Nakamoto Vol.2 p.306"),
    ("bpy",  "chelate", "any", "Raman",1590, 1640, "very strong", "C=C / C=N ring stretch",      False, True,  "Nakamoto Vol.2 p.304"),
    ("bpy",  "chelate", "any", "Raman", 990, 1020, "very strong", "ring breathing mode",         False, True,  "Nakamoto Vol.2 p.305"),
    ("bpy",  "chelate", "any", "Raman", 390,  450, "medium",      "M–N stretch",                 False, True,  "Nakamoto Vol.2 p.306"),

    # =========================================================================
    # TERPYRIDINE  tpy / terpy   (Nakamoto Vol.2 pp. 302–315)
    # Tridentate version of bipyridine. Three pyridine rings.
    # Very similar band pattern to bipy but slightly shifted due to
    # the additional ring and meridional coordination mode.
    # =========================================================================
    ("tpy",  "chelate", "any", "IR",   3000, 3090, "medium",      "C–H stretch (aromatic)",      True,  False, "Nakamoto Vol.2 p.303"),
    ("tpy",  "chelate", "any", "IR",   1580, 1630, "strong",      "C=C / C=N ring stretch",      True,  False, "Nakamoto Vol.2 p.304"),
    ("tpy",  "chelate", "any", "IR",   1460, 1510, "strong",      "C=C ring stretch",            True,  False, "Nakamoto Vol.2 p.304"),
    ("tpy",  "chelate", "any", "IR",   1390, 1430, "medium",      "C–H in-plane bend",           True,  False, "Nakamoto Vol.2 p.305"),
    ("tpy",  "chelate", "any", "IR",    770,  810, "strong",      "C–H out-of-plane bend",       True,  False, "Nakamoto Vol.2 p.305"),
    ("tpy",  "chelate", "any", "IR",    380,  440, "medium",      "M–N stretch",                 True,  False, "Nakamoto Vol.2 p.306"),
    # Raman
    ("tpy",  "chelate", "any", "Raman",1580, 1630, "very strong", "C=C / C=N ring stretch",      False, True,  "Nakamoto Vol.2 p.304"),
    ("tpy",  "chelate", "any", "Raman",1000, 1030, "very strong", "ring breathing mode",         False, True,  "Nakamoto Vol.2 p.305"),
    ("tpy",  "chelate", "any", "Raman", 380,  440, "medium",      "M–N stretch",                 False, True,  "Nakamoto Vol.2 p.306"),
    # terpy alias
    ("terpy","chelate", "any", "IR",   1580, 1630, "strong",      "C=C / C=N ring stretch",      True,  False, "Nakamoto Vol.2 p.304"),
    ("terpy","chelate", "any", "IR",    770,  810, "strong",      "C–H out-of-plane bend",       True,  False, "Nakamoto Vol.2 p.305"),
    ("terpy","chelate", "any", "IR",    380,  440, "medium",      "M–N stretch",                 True,  False, "Nakamoto Vol.2 p.306"),
    ("terpy","chelate", "any", "Raman",1580, 1630, "very strong", "C=C / C=N ring stretch",      False, True,  "Nakamoto Vol.2 p.304"),
    ("terpy","chelate", "any", "Raman",1000, 1030, "very strong", "ring breathing mode",         False, True,  "Nakamoto Vol.2 p.305"),

    # =========================================================================
    # EDTA  ethylenediaminetetraacetato   (Nakamoto Vol.2 pp. 195–230)
    # Hexadentate ligand — four carboxylate arms + two amine nitrogens.
    # Dominant bands: carboxylate C=O stretches (very strong),
    # C-N stretches, and M-N/M-O stretches.
    # =========================================================================
    ("EDTA", "chelate", "any", "IR",   1580, 1650, "very strong", "COO⁻ asym. stretch",          True,  False, "Nakamoto Vol.2 p.218"),
    ("EDTA", "chelate", "any", "IR",   1380, 1420, "strong",      "COO⁻ sym. stretch",           True,  False, "Nakamoto Vol.2 p.218"),
    ("EDTA", "chelate", "any", "IR",   2850, 2980, "medium",      "C–H stretch (CH₂)",           True,  False, "Nakamoto Vol.2 p.219"),
    ("EDTA", "chelate", "any", "IR",   3150, 3300, "medium",      "N–H stretch",                 True,  False, "Nakamoto Vol.2 p.219"),
    ("EDTA", "chelate", "any", "IR",   1280, 1350, "medium",      "C–N stretch",                 True,  False, "Nakamoto Vol.2 p.220"),
    ("EDTA", "chelate", "any", "IR",    500,  600, "medium",      "M–O stretch (carboxylate)",   True,  False, "Nakamoto Vol.2 p.221"),
    ("EDTA", "chelate", "any", "IR",    350,  450, "medium",      "M–N stretch",                 True,  False, "Nakamoto Vol.2 p.221"),
    # Raman
    ("EDTA", "chelate", "any", "Raman",1380, 1420, "strong",      "COO⁻ sym. stretch",           False, True,  "Nakamoto Vol.2 p.218"),
    ("EDTA", "chelate", "any", "Raman",1280, 1350, "medium",      "C–N stretch",                 False, True,  "Nakamoto Vol.2 p.220"),
    ("EDTA", "chelate", "any", "Raman", 500,  600, "strong",      "M–O stretch (carboxylate)",   False, True,  "Nakamoto Vol.2 p.221"),
    ("EDTA", "chelate", "any", "Raman", 350,  450, "medium",      "M–N stretch",                 False, True,  "Nakamoto Vol.2 p.221"),
    # edta lowercase alias
    ("edta", "chelate", "any", "IR",   1580, 1650, "very strong", "COO⁻ asym. stretch",          True,  False, "Nakamoto Vol.2 p.218"),
    ("edta", "chelate", "any", "IR",   1380, 1420, "strong",      "COO⁻ sym. stretch",           True,  False, "Nakamoto Vol.2 p.218"),
    ("edta", "chelate", "any", "IR",    500,  600, "medium",      "M–O stretch (carboxylate)",   True,  False, "Nakamoto Vol.2 p.221"),
    ("edta", "chelate", "any", "IR",    350,  450, "medium",      "M–N stretch",                 True,  False, "Nakamoto Vol.2 p.221"),
    ("edta", "chelate", "any", "Raman",1380, 1420, "strong",      "COO⁻ sym. stretch",           False, True,  "Nakamoto Vol.2 p.218"),
    ("edta", "chelate", "any", "Raman", 500,  600, "strong",      "M–O stretch (carboxylate)",   False, True,  "Nakamoto Vol.2 p.221"),

    # =========================================================================
    # CYCLOPENTADIENYL  Cp   (Nakamoto Vol.2 pp. 334–360)
    # η5-ligand — bonds through all five carbon atoms (hapticity 5).
    # Very characteristic: C-H stretch ~3080 cm⁻¹, ring stretch ~1400 cm⁻¹,
    # C-H out-of-plane bend ~800 cm⁻¹, and M-ring tilt/stretch below 500.
    # Ferrocene [Fe(Cp)2] is the benchmark complex.
    # =========================================================================
    ("Cp",   "chelate", "any", "IR",   3060, 3110, "medium",      "C–H stretch (Cp ring)",       True,  False, "Nakamoto Vol.2 p.336"),
    ("Cp",   "chelate", "any", "IR",   1390, 1440, "strong",      "C=C ring stretch",            True,  False, "Nakamoto Vol.2 p.337"),
    ("Cp",   "chelate", "any", "IR",   1000, 1030, "medium",      "C–H in-plane bend",           True,  False, "Nakamoto Vol.2 p.338"),
    ("Cp",   "chelate", "any", "IR",    780,  830, "very strong", "C–H out-of-plane bend",       True,  False, "Nakamoto Vol.2 p.338"),
    ("Cp",   "chelate", "any", "IR",    440,  500, "medium",      "M–ring tilt",                 True,  False, "Nakamoto Vol.2 p.340"),
    ("Cp",   "chelate", "any", "IR",    300,  400, "medium",      "M–ring stretch",              True,  False, "Nakamoto Vol.2 p.340"),
    # Fe specific — ferrocene benchmark
    ("Cp",   "chelate", "Fe", "IR",     478,  478, "medium",      "M–ring tilt (ferrocene)",     True,  False, "Nakamoto Vol.2 p.341"),
    ("Cp",   "chelate", "Fe", "IR",     309,  309, "medium",      "M–ring stretch (ferrocene)",  True,  False, "Nakamoto Vol.2 p.341"),
    # Raman — ring breathing and M-ring stretch are particularly strong
    ("Cp",   "chelate", "any", "Raman",3060, 3110, "medium",      "C–H stretch (Cp ring)",       False, True,  "Nakamoto Vol.2 p.336"),
    ("Cp",   "chelate", "any", "Raman",1390, 1440, "strong",      "C=C ring stretch",            False, True,  "Nakamoto Vol.2 p.337"),
    ("Cp",   "chelate", "any", "Raman", 900,  940, "very strong", "ring breathing mode",         False, True,  "Nakamoto Vol.2 p.338"),
    ("Cp",   "chelate", "any", "Raman", 300,  400, "very strong", "M–ring stretch",              False, True,  "Nakamoto Vol.2 p.340"),
    ("Cp",   "chelate", "Fe", "Raman",  309,  309, "very strong", "M–ring stretch (ferrocene)",  False, True,  "Nakamoto Vol.2 p.341"),

    # =========================================================================
    # TRIMETHYLPHOSPHINE  PMe3   (Nakamoto Vol.2 pp. 326–333; Socrates p.221)
    # Monodentate P-donor ligand. Key bands: P-C stretch ~700 cm⁻¹,
    # C-H stretches ~2900 cm⁻¹, and M-P stretch (Raman active) ~200–350 cm⁻¹.
    # =========================================================================
    ("PMe3", "terminal", "any", "IR",  2880, 2980, "medium",      "C–H stretch (CH₃)",           True,  False, "Socrates p.221"),
    ("PMe3", "terminal", "any", "IR",  1400, 1460, "medium",      "CH₃ deformation",             True,  False, "Socrates p.221"),
    ("PMe3", "terminal", "any", "IR",   920,  980, "medium",      "CH₃ rock",                    True,  False, "Socrates p.222"),
    ("PMe3", "terminal", "any", "IR",   680,  730, "strong",      "P–C stretch",                 True,  False, "Nakamoto Vol.2 p.328"),
    ("PMe3", "terminal", "any", "IR",   280,  350, "medium",      "M–P stretch",                 True,  False, "Nakamoto Vol.2 p.330"),
    # Raman — M-P stretch very strong
    ("PMe3", "terminal", "any", "Raman",2880, 2980, "medium",     "C–H stretch (CH₃)",           False, True,  "Socrates p.221"),
    ("PMe3", "terminal", "any", "Raman", 680,  730, "strong",     "P–C stretch",                 False, True,  "Nakamoto Vol.2 p.328"),
    ("PMe3", "terminal", "any", "Raman", 280,  350, "very strong","M–P stretch",                 False, True,  "Nakamoto Vol.2 p.330"),

    # =========================================================================
    # TRIETHYLPHOSPHINE  PEt3   (Nakamoto Vol.2 pp. 326–333; Socrates p.222)
    # Similar to PMe3 but with ethyl groups. P-C stretch shifts slightly.
    # Additional CH₂ bands from the ethyl chains.
    # =========================================================================
    ("PEt3", "terminal", "any", "IR",  2850, 2980, "medium",      "C–H stretch (CH₂/CH₃)",      True,  False, "Socrates p.222"),
    ("PEt3", "terminal", "any", "IR",  1440, 1470, "medium",      "CH₂ scissor",                 True,  False, "Socrates p.222"),
    ("PEt3", "terminal", "any", "IR",  1370, 1400, "medium",      "CH₃ deformation",             True,  False, "Socrates p.222"),
    ("PEt3", "terminal", "any", "IR",   700,  760, "strong",      "P–C stretch",                 True,  False, "Nakamoto Vol.2 p.329"),
    ("PEt3", "terminal", "any", "IR",   250,  330, "medium",      "M–P stretch",                 True,  False, "Nakamoto Vol.2 p.330"),
    # Raman
    ("PEt3", "terminal", "any", "Raman",2850, 2980, "medium",     "C–H stretch (CH₂/CH₃)",      False, True,  "Socrates p.222"),
    ("PEt3", "terminal", "any", "Raman", 700,  760, "strong",     "P–C stretch",                 False, True,  "Nakamoto Vol.2 p.329"),
    ("PEt3", "terminal", "any", "Raman", 250,  330, "very strong","M–P stretch",                 False, True,  "Nakamoto Vol.2 p.330"),

    # =========================================================================
    # METHYL  CH3⁻   (Nakamoto Vol.2 pp. 361–375; Socrates p.48)
    # Ligand alkyle sigma-donneur. Bandes diagnostiques :
    # C-H stretch ~2900 cm⁻¹, déformation CH3 ~1200 cm⁻¹,
    # M-C stretch ~400–600 cm⁻¹ (très diagnostic).
    # =========================================================================
    ("CH3", "terminal", "any", "IR",   2900, 2990, "strong",      "C–H stretch (CH₃)",           True,  False, "Nakamoto Vol.2 p.362"),
    ("CH3", "terminal", "any", "IR",   2800, 2870, "medium",      "C–H sym. stretch (CH₃)",      True,  False, "Nakamoto Vol.2 p.362"),
    ("CH3", "terminal", "any", "IR",   1380, 1460, "medium",      "CH₃ sym. deformation",        True,  False, "Nakamoto Vol.2 p.363"),
    ("CH3", "terminal", "any", "IR",   1150, 1250, "medium",      "CH₃ rock",                    True,  False, "Nakamoto Vol.2 p.364"),
    ("CH3", "terminal", "any", "IR",    400,  600, "strong",      "M–C stretch",                 True,  False, "Nakamoto Vol.2 p.365"),
    # Metal specific
    ("CH3", "terminal", "Fe", "IR",     510,  540, "strong",      "Fe–C stretch",                True,  False, "Nakamoto Vol.2 p.367"),
    ("CH3", "terminal", "Pt", "IR",     550,  590, "strong",      "Pt–C stretch",                True,  False, "Nakamoto Vol.2 p.368"),
    ("CH3", "terminal", "Co", "IR",     480,  520, "strong",      "Co–C stretch",                True,  False, "Nakamoto Vol.2 p.367"),
    # Raman
    ("CH3", "terminal", "any", "Raman", 2900, 2990, "strong",     "C–H stretch (CH₃)",           False, True,  "Nakamoto Vol.2 p.362"),
    ("CH3", "terminal", "any", "Raman", 1380, 1460, "medium",     "CH₃ sym. deformation",        False, True,  "Nakamoto Vol.2 p.363"),
    ("CH3", "terminal", "any", "Raman",  400,  600, "very strong","M–C stretch",                 False, True,  "Nakamoto Vol.2 p.365"),
]


# ---------------------------------------------------------------------------
# Database class
# ---------------------------------------------------------------------------

class IRBandDB:
    """
    Interface to the IR/Raman band SQLite database.

    Parameters
    ----------
    db_path : str, optional
        Path to the SQLite file. Defaults to ':memory:' (in-memory).

    """

    def __init__(self, db_path: str = ":memory:"):
        self.db_path = db_path
        self._conn   = sqlite3.connect(db_path)
        self._conn.row_factory = sqlite3.Row
        self._create_table()
        self._seed()

    # ------------------------------------------------------------------
    # Setup
    # ------------------------------------------------------------------

    def _create_table(self) -> None:
        self._conn.execute("""
            CREATE TABLE IF NOT EXISTS ir_ra_bands (
                id              INTEGER PRIMARY KEY AUTOINCREMENT,
                ligand          TEXT    NOT NULL,
                coordination    TEXT    NOT NULL,
                metal           TEXT    NOT NULL DEFAULT 'any',
                spectrum_type   TEXT    NOT NULL,
                wn_min          REAL    NOT NULL,
                wn_max          REAL    NOT NULL,
                intensity       TEXT    NOT NULL,
                assignment      TEXT    NOT NULL,
                ir_active       INTEGER NOT NULL DEFAULT 1,
                raman_active    INTEGER NOT NULL DEFAULT 0,
                source          TEXT
            )
        """)
        self._conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_ligand ON ir_ra_bands(ligand)"
        )
        self._conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_spectrum ON ir_ra_bands(spectrum_type)"
        )
        self._conn.commit()

    def _seed(self) -> None:
        """Insert seed data only if the table is empty."""
        count = self._conn.execute("SELECT COUNT(*) FROM ir_ra_bands").fetchone()[0]
        if count > 0:
            return
        self._conn.executemany("""
            INSERT INTO ir_ra_bands
              (ligand, coordination, metal, spectrum_type,
               wn_min, wn_max, intensity, assignment,
               ir_active, raman_active, source)
            VALUES (?,?,?,?,?,?,?,?,?,?,?)
        """, SEED_BANDS)
        self._conn.commit()

    # ------------------------------------------------------------------
    # Public query interface
    # ------------------------------------------------------------------

    def get_bands(
        self,
        ligand:          str,
        spectrum_type:   str           = "IR",
        coordination:    Optional[str] = None,
        metal:           Optional[str] = None,
        oxidation_state: Optional[int] = None,   # reserved for future use
    ) -> list[BandRecord]:
        """
        Retrieve band records for a ligand.

        Parameters
        ----------
        ligand : str
            Ligand formula as used in the parser, e.g. "CN", "NH3", "en".
        spectrum_type : str
            "IR" or "Raman".
        coordination : str, optional
            Filter by coordination mode: "terminal", "bridging", "chelate",
            "free". If None, returns all modes.
        metal : str, optional
            Metal symbol e.g. "Fe". Returns metal-specific rows first,
            then generic "any" rows, deduplicating by assignment.

        Returns
        -------
        list[BandRecord]
            Sorted by wavenumber center, most specific match first.
        """
        query  = """
            SELECT * FROM ir_ra_bands
            WHERE ligand = ?
              AND spectrum_type = ?
        """
        params: list = [ligand, spectrum_type]

        if coordination:
            query  += " AND (coordination = ? OR coordination = 'any')"
            params.append(coordination)

        query += " ORDER BY wn_min ASC"

        rows = self._conn.execute(query, params).fetchall()

        # If a specific metal is given, prefer metal-specific rows
        # and fall back to "any" rows only when no specific row exists
        if metal:
            specific   = [r for r in rows if r["metal"] == metal]
            generic    = [r for r in rows if r["metal"] == "any"]

            # Keep only generic rows whose band family is not already covered
            # by a metal-specific row. For example, Ni-Cl stretch should
            # replace the generic M-Cl stretch in the app, not duplicate it.
            covered = {
                self._band_family_key(r, metal)
                for r in specific
            }
            supplement = [
                r for r in generic
                if self._band_family_key(r, metal) not in covered
            ]

            rows = specific + supplement

        return [self._row_to_record(r) for r in rows]

    @staticmethod
    def _band_family_key(row: sqlite3.Row, metal: str) -> tuple[str, str]:
        """Return a normalized key for matching specific and generic bands."""
        assignment = row["assignment"]
        assignment = assignment.replace(f"{metal}-", "M-")
        assignment = assignment.replace(f"{metal}–", "M–")
        assignment = re.sub(r"\([^)]*\)", "", assignment)
        assignment = re.sub(r"\s+", " ", assignment).strip()
        return row["coordination"], assignment

    def get_all_ligands(self) -> list[str]:
        """Return all ligand symbols present in the database."""
        rows = self._conn.execute(
            "SELECT DISTINCT ligand FROM ir_ra_bands ORDER BY ligand"
        ).fetchall()
        return [r[0] for r in rows]

    def get_bands_in_range(
        self,
        wn_low:        float,
        wn_high:       float,
        spectrum_type: str = "IR",
    ) -> list[BandRecord]:
        """
        Find all bands whose range overlaps [wn_low, wn_high].
        Useful for identifying what ligand a mystery peak might belong to.
        """
        rows = self._conn.execute("""
            SELECT * FROM ir_ra_bands
            WHERE spectrum_type = ?
              AND wn_max >= ?
              AND wn_min <= ?
            ORDER BY wn_min ASC
        """, (spectrum_type, wn_low, wn_high)).fetchall()
        return [self._row_to_record(r) for r in rows]

    def add_band(
        self,
        ligand:        str,
        coordination:  str,
        metal:         str,
        spectrum_type: str,
        wn_min:        float,
        wn_max:        float,
        intensity:     str,
        assignment:    str,
        ir_active:     bool = True,
        raman_active:  bool = False,
        source:        str  = "user",
    ) -> None:
        """Add a custom band entry to the database."""
        self._conn.execute("""
            INSERT INTO ir_ra_bands
              (ligand, coordination, metal, spectrum_type,
               wn_min, wn_max, intensity, assignment,
               ir_active, raman_active, source)
            VALUES (?,?,?,?,?,?,?,?,?,?,?)
        """, (ligand, coordination, metal, spectrum_type,
              wn_min, wn_max, intensity, assignment,
              int(ir_active), int(raman_active), source))
        self._conn.commit()

    def summary(self, ligand: str) -> None:
        """Pretty-print all bands for a ligand to the console."""
        print(f"\n{'='*75}")
        print(f"  IR/Raman bands for ligand: {ligand}")
        print(f"{'='*75}")
        print(f"{'Coord.':<12} {'Type':<6} {'Range (cm⁻¹)':<18} {'Intensity':<14} Assignment")
        print(f"{'-'*75}")
        for stype in ("IR", "Raman"):
            bands = self.get_bands(ligand, spectrum_type=stype)
            for b in bands:
                print(
                    f"{b.coordination:<12} {stype:<6} "
                    f"{b.wn_min:>6.0f}–{b.wn_max:<6.0f}   "
                    f"{b.intensity:<14} {b.assignment}"
                )
        print()

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _row_to_record(self, row: sqlite3.Row) -> BandRecord:
        return BandRecord(
            ligand        = row["ligand"],
            coordination  = row["coordination"],
            metal         = row["metal"],
            spectrum_type = row["spectrum_type"],
            wn_min        = row["wn_min"],
            wn_max        = row["wn_max"],
            intensity     = row["intensity"],
            assignment    = row["assignment"],
            ir_active     = bool(row["ir_active"]),
            raman_active  = bool(row["raman_active"]),
            source        = row["source"] or "",
        )

    def close(self) -> None:
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
