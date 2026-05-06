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

Usage
-----
    from data.database.ir_ra_bands import IRBandDB

    db = IRBandDB()                      # uses in-memory DB (no file needed)
    bands = db.get_bands("CN")           # all IR bands for cyanide
    bands = db.get_bands("NH3", spectrum_type="Raman")
    bands = db.get_bands("CO", metal="Cr", oxidation_state=0)
    db.summary("CN")                     # pretty-print to console
"""

import sqlite3
from dataclasses import dataclass
from typing import Optional
import os


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

            # Keep only generic rows whose assignment isn't already covered
            covered    = {r["assignment"] for r in specific}
            supplement = [r for r in generic if r["assignment"] not in covered]

            rows = specific + supplement

        return [self._row_to_record(r) for r in rows]

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
