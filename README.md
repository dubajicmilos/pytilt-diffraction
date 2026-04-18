# pytilt-diffraction

Interactive single-crystal X-ray diffraction simulator for perovskites with
tunable Glazer octahedral-tilt systems. A matplotlib GUI that lets you pick a
tilt pattern (`a+b-b-`, `a0a0c-`, ...) and watch the reciprocal-space pattern
update as you drag the tilt-angle sliders, switch zone axes, or step through
HK_L layers.

Built on top of the [pytilting](https://gitlab.com/pyseries/pytilting)
tilt-generator (vendored under `vendor/pytilting/`, GPL v2). The structure
factor / `F(hkl)` engine and the GUI are original.

## Install

```bash
git clone <this-repo-url>
cd pytilt-diffraction
pip install -e .
```

The pytilting tilt generator is already bundled in `vendor/pytilting/`.

## Run

```bash
python -m pytilt_diffraction.simulator
# or, after install:
pytilt-gui
```

Optional composition override:

```bash
python -m pytilt_diffraction.simulator Cs Pb Br 5.874
```

## Web app (Streamlit)

A browser-hosted version of the simulator lives in `streamlit_app.py`. It
reuses the same physics core (vendored pytilting + `DiffractionCalculator`);
only the widget layer is swapped for Streamlit controls.

```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```

Deploy: push to GitHub, connect the repo at
[share.streamlit.io](https://share.streamlit.io), and it auto-redeploys on
every commit.

## Controls

- **Glazer tilt system** (radio, two columns): 18 of the 23 Howard-Stokes
  systems mapped to their aristotype space group.
- **omega_x / omega_y / omega_z** (sliders): tilt magnitudes in degrees.
  Magnitudes that the Glazer letters demand equal are kept tied automatically.
- **HK_L layer** (slider): integer layer in the supercell reciprocal lattice.
  For a 2x2x2 supercell, layer 1 == 0.5 in the parent pseudocubic cell
  (shows R-point / superlattice reflections).
- **Zone axis** (radio): view direction of the reciprocal-space slice.
- **Material** (radio, two columns): 16 ABX3 perovskite presets covering
  halides (CsPbCl3, CsPbBr3, CsPbI3, CsSnBr3, CsGeBr3, RbPbBr3), oxides
  (SrTiO3, BaTiO3, CaTiO3, PbTiO3, BaZrO3, LaAlO3, LaMnO3, LaFeO3), and
  fluorides (KMgF3, KNiF3). Picking a preset swaps the basis and snaps the
  lattice constant to the high-temperature cubic aristotype value. The
  tilt, zone, and layer selections are preserved.
- **a0** (slider, 3.0 - 7.5 A): fine-tune the cubic lattice constant without
  changing the composition. Useful for temperature / composition sweeps
  around a given preset.
- **d_min / h_max / spot size / label threshold** (sliders): viewing
  parameters.
- **Powder** (button): open a side window with a kinematic powder pattern
  for the current composition, lattice, and tilt. Picks Cu / Mo / Co / Cr /
  Ag K-alpha wavelengths, with sliders for 2theta range, h_max, and peak
  FWHM. Updates automatically when you touch the main window.
- **Save PNG / Export hkl / Reset**: export buttons flash a confirmation
  banner at the top of the window; files land next to the script.

## Package layout

```
pytilt-diffraction/
    pytilt_diffraction/
        __init__.py
        calculator.py   # CIFParser + DiffractionCalculator
        simulator.py    # interactive matplotlib GUI
    vendor/
        pytilting/      # upstream Glazer-tilt generator, vendored (GPL v2)
    tests/
        test_physics.py
        test_diffsims_parity.py
    examples/
    pyproject.toml
    README.md
    LICENSE
```

## Licensing

This project depends on and redistributes **pytilting**, which is GPL v2.
Accordingly, this project is distributed under **GPL v2 or later** (see
`LICENSE`). Any redistribution must preserve `vendor/pytilting/LICENSE`.

## Validation

The `DiffractionCalculator` has been cross-checked against
[diffsims](https://github.com/pyxem/diffsims) (pyxem's X-ray / electron
diffraction simulation package). We compute F(hkl) with our Cromer-Mann
parameterisation and compare against diffsims'
`get_kinematical_atomic_scattering_factor` (Doyle-Turner 1968 table) with an
explicit atomic sum, over a 7-material / 11-hkl grid:

- **Halides**: CsPbCl3, CsPbBr3, CsPbI3
- **Oxides**: SrTiO3, BaTiO3, LaAlO3
- **Fluorides**: KMgF3

Across 77 reflections the magnitude agreement is typically under 1% and
always under 5% (absolute slack 0.15 for near-zero structure factors), and
the top-5 strongest reflections match ordering for every material. The
residual comes from the two different form-factor parameterisations, not
from the summation. Our calculator is also ~3x faster per reflection than
diffsims' form-factor path with the same summation (37 us vs 117 us on
cubic CsPbBr3).

Tests and timings: `tests/test_diffsims_parity.py`, `tests/test_physics.py`,
and `tests/BENCHMARK.md`.

```bash
pip install diffsims pytest
python -m pytest tests/ -q
```

A footnote on diffsims 0.7.0: its top-level `get_kinematical_structure_factor`
has a bug in `find_asymmetric_positions` that silently drops all but the
first asymmetric-unit atom. We therefore compare form-factor
parameterisations by summing over atoms ourselves, which is the fair
comparison. See `tests/BENCHMARK.md` for details.
