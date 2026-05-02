# Function Reference — streamlit_app1.py

Quick reference to understand what each function does, written to help with the notebook.

---

## `gaussian(x, center, height, sigma)`

Computes a single Gaussian curve evaluated over the array `x`.

- `center` — wavenumber position of the peak (cm⁻¹)
- `height` — peak amplitude, between 0 and 1
- `sigma` — controls the width (standard deviation in cm⁻¹)

This is the core mathematical primitive of the spectrum renderer.
Every literature band becomes one Gaussian curve; overlapping curves are summed to form the full spectrum.

---

## `build_spectrum(bands, sigma, wn_range)`

Takes a list of `BandRecord` objects from the database and builds a composite spectrum.

Steps:
1. Creates a wavenumber axis `x` from `wn_range[0]` to `wn_range[1]` (3000 points).
2. For each band, reads its center (`(wn_min + wn_max) / 2`) and maps its intensity label
   (`"very strong"`, `"strong"`, …) to a numeric height using `INTENSITY_SCALE`.
3. Calls `gaussian()` and adds the result to the running sum `y`.
4. Normalizes `y` to [0, 1] so all spectra share the same vertical scale.

Returns `(x, y)` — two NumPy arrays ready to plot.

---

## `plot_spectrum(bands, title, sigma)`

Wraps `build_spectrum()` into a Plotly figure.

- Draws the composite spectrum as a filled line trace.
- Adds a vertical dashed marker at each band center so the reader can see
  which literature position each Gaussian corresponds to.
- Sets the x-axis in reverse order (convention in IR spectroscopy: high wavenumber on the left).

Returns a `plotly.graph_objects.Figure`.

---

## `bands_to_df(bands)`

Converts a list of `BandRecord` objects into a pandas `DataFrame` for display in Streamlit.

Each row contains: ligand formula, coordination mode, center wavenumber, full range, intensity label, and vibrational assignment.
Rows are sorted by center wavenumber (ascending) so the table reads like a spectrum from low to high frequency.

---

## App flow (main script, no function)

| Step | What happens |
|------|-------------|
| User types a formula and clicks **Analyze** | `st.stop()` prevents any output before the button is pressed |
| `parse_formula(formula)` | Calls the `coordchem` parser; extracts metal, ligands, charges, denticity, donor atoms |
| `geometry_report(parsed)` | Calls the geometry module; returns oxidation state, d-count, predicted geometry |
| `IRBandDB().get_bands(lig, spectrum_type, metal)` | Queries the SQLite database for every ligand in the complex; returns metal-specific rows first, then generic rows |
| `plot_spectrum(all_bands, ...)` | Renders the combined Gaussian spectrum with Plotly |
| `bands_to_df(all_bands)` | Displays the raw band data as a sortable table |

---

## Key objects from `coordchem`

### `ParsedComplex` (from `parser.py`)
Dataclass returned by `parse_formula()`. Useful fields:

| Field | Type | Content |
|-------|------|---------|
| `metal` | `str` | Metal symbol, e.g. `"Fe"` |
| `ligands` | `dict[str, int]` | `{"CN": 6}` |
| `complex_charge` | `int` | Charge on the complex ion |
| `oxidation_state` | `int` | Computed from charge balance |
| `coordination_number` | `int` | Sum of all ligand denticities |
| `ligand_names` | `dict[str, str]` | IUPAC names |
| `ligand_charges` | `dict[str, int]` | Per-ligand charge |
| `ligand_denticity` | `dict[str, int]` | Bite number |
| `donor_atoms` | `dict[str, str]` | Donor element symbol |
| `warnings` / `errors` | `list[str]` | Parsing issues |

### `BandRecord` (from `database/ir_bands.py`)
Dataclass returned by `IRBandDB.get_bands()`. Useful fields:

| Field | Type | Content |
|-------|------|---------|
| `ligand` | `str` | e.g. `"CN"` |
| `coordination` | `str` | `"terminal"`, `"bridging"`, `"chelate"`, `"free"` |
| `metal` | `str` | `"any"` or specific symbol |
| `spectrum_type` | `str` | `"IR"` or `"Raman"` |
| `wn_min` / `wn_max` | `float` | Wavenumber range (cm⁻¹) |
| `center` | `float` | `(wn_min + wn_max) / 2` (property) |
| `intensity` | `str` | `"very strong"`, `"strong"`, `"medium"`, `"weak"` |
| `assignment` | `str` | e.g. `"C≡N stretch"` |
| `source` | `str` | Bibliographic reference |
