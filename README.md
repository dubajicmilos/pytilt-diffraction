
# pytilt-diffraction

Interactive single-crystal X-ray diffraction simulator for perovskites with
tunable Glazer octahedral-tilt systems. In the matplotlib GUI you choose a
tilt system, and the reciprocal-space pattern updates as you drag the
tilt-angle sliders, switch zone axes, or step through HKL layers.

The GUI offers the 15 distinct tilt systems in the group-subgroup tree of
Howard & Stokes (1998), a subset of Glazer's (1972) 23 tilt systems:

<img width="1429" height="736" alt="tilts_1" src="https://github.com/user-attachments/assets/e15043ec-38ca-4d6b-aefc-cc3599684918" />

The simulator is built on the [pytilting](https://gitlab.com/pyseries/pytilting)
tilt generator, which is vendored under `vendor/pytilting/` (GPL v2). The
structure-factor engine (`F(hkl)`) and the GUI are original.

**References**

- Glazer, A. M. (1972). *The classification of tilted octahedra in
  perovskites*. Acta Cryst. **B28**, 3384-3392.
- Howard, C. J. & Stokes, H. T. (1998). *Group-Theoretical Analysis of
  Octahedral Tilting in Perovskites*. Acta Cryst. **B54**, 782-789.

## Install

```bash
git clone <this-repo-url>
cd pytilt-diffraction
pip install -e .
```

The pytilting tilt generator is bundled in `vendor/pytilting/`, so it does
not need a separate install.

## Run

```bash
python -m pytilt_diffraction.simulator
# or, after install:
pytilt-gui
```


## Screenshots

The two screenshots below show CsPbI3 with `a+a+a+`, alpha = 6 deg,
zone [001], L = 0 (parent rlu).

| log intensity | linear intensity |
|---|---|
| ![log mode](docs/screenshots/log_after.png) | ![linear mode](docs/screenshots/linear_after.png) |

The next screenshot shows a half-integer parent layer: `L_super = 1`, which
is L = 0.5 in the parent pseudocubic cell. The 2x2x2 supercell halves the
reciprocal-lattice spacing along each axis, so odd `L_super` slices fall
*between* the parent Bragg peaks. Every visible peak on such a slice is a
superlattice reflection produced by the octahedral tilts:

![log mode, L_super=1](docs/screenshots/log_after_L1.png)

Marker colour (white to dark blue, `Blues` colormap) and marker size both
scale with reflection intensity, so darker and larger markers mean stronger
reflections. Linear mode applies a mild gamma compression so that weak peaks
stay visible next to the dominant Bragg reflections; log mode uses the
log-stretched intensities directly.

## Generic CIF viewer

For single-crystal diffraction from an arbitrary CIF, without the Glazer-tilt
machinery or the perovskite assumption, there is a separate GUI built on the
same `DiffractionCalculator`:

```bash
python -m pytilt_diffraction.cif_viewer path/to/structure.cif
# or, with no path -> file picker
python -m pytilt_diffraction.cif_viewer
```

It has the same controls as the Glazer simulator (zone axis, HKL layer
slider, d_min, h_max, log/linear, twin (3 domains), labels).

The HKL layer slider works for every zone axis. It sets the integer
constant `L` in the zone law `h*u + k*v + l*w = L`, so:

  - zone [001]: layer L  =>  hk-plane at l = L
  - zone [100]: layer L  =>  kl-plane at h = L
  - zone [110]: layer L  =>  diagonal slice  h + k = L
  - zone [111]: layer L  =>  diagonal slice  h + k + l = L

The plot title prints the zone-law constraint explicitly, so it is clear
which slice is shown.

The viewer has three buttons for saving results:

- **Save PNG**: the current pattern.
- **Export hkl**: `(h, k, l, d, |F|, I, I_norm)` table as `.txt`.
- **Export 2D matrix**: rasterises the visible slice as a regular
  `n_pix x n_pix` array of summed Gaussians, written as both
  `.npy` and `.csv` (plus a `*_meta.txt` file describing the extent and
  shape). Use it for side-by-side comparison with experimental detector
  images, or load it into MATLAB / Origin / Igor.

![CIF viewer](docs/screenshots/cif_viewer_log.png)

The `export grid (n_pix)` slider at the bottom of the window sets the
resolution of the exported matrix (64 to 1024 pixels per side).

## Web app (Streamlit)

Live demo (no install required):

- Glazer simulator: <https://pytilt-diffraction-milos.streamlit.app/>
- CIF viewer (upload a `.cif`):
  <https://pytilt-diffraction-milos.streamlit.app/CIF_viewer>

The browser versions of both modes are in `streamlit_app.py` (Glazer
simulator) and `pages/2_CIF_viewer.py` (CIF upload). They use the same
`DiffractionCalculator` as the desktop GUIs; only the widget layer is
Streamlit.

```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```

Streamlit detects the `pages/` folder automatically, so the left-hand
sidebar lists two pages: the Glazer simulator (default) and the CIF
viewer, where you can drag and drop a `.cif` file and download the hkl
table and the 2D matrix, as in the desktop GUI.

The live demo is deployed on [share.streamlit.io](https://share.streamlit.io)
and redeploys automatically on every push to `main`.

## Controls

The desktop Glazer simulator has these controls:

- **Glazer tilt system** (radio, two columns): the 15 distinct tilt
  systems of Howard & Stokes (1998), labelled with the standard
  hettotype space-group symbol. Glazer (1972) listed 23 tilt systems;
  Howard & Stokes showed that 8 of these are crystallographically
  equivalent to others, which leaves 15.
- **omega_x / omega_y / omega_z** (sliders): tilt magnitudes in degrees.
  Magnitudes that the Glazer letters require to be equal move together
  automatically.
- **HKL layer** (slider): integer L index in the supercell reciprocal
  lattice, which selects the constant-L slice through reciprocal space
  that the pattern shows. For the default 2x2x2 supercell, `L_super = 1`
  corresponds to L = 0.5 in the parent pseudocubic cell (the R-point /
  superlattice layer).
- **Zone axis** (radio): view direction of the reciprocal-space slice.
- **Material** (radio, two columns): 16 ABX3 perovskite presets covering
  halides (CsPbCl3, CsPbBr3, CsPbI3, CsSnBr3, CsGeBr3, RbPbBr3), oxides
  (SrTiO3, BaTiO3, CaTiO3, PbTiO3, BaZrO3, LaAlO3, LaMnO3, LaFeO3), and
  fluorides (KMgF3, KNiF3). Choosing a preset swaps the basis and sets the
  lattice constant to the high-temperature cubic aristotype value. The
  tilt, zone, and layer selections are preserved.
- **a0** (slider, 3.0 to 7.5 A): fine-tunes the cubic lattice constant
  without changing the composition. Useful for temperature or composition
  sweeps around a given preset.
- **d_min / h_max / spot size / label threshold** (sliders): viewing
  parameters.
- **Powder** (button): opens a side window with a kinematic powder pattern
  for the current composition, lattice constant, and tilt. You can choose
  Cu / Mo / Co / Cr / Ag K-alpha radiation, and sliders set the 2theta
  range, h_max, and peak FWHM. The powder pattern updates automatically
  when the composition, lattice constant, or tilt changes in the main
  window.
- **Save PNG / Export hkl / Reset**: the two export buttons show a
  confirmation banner at the top of the window and write their files to
  the `pytilt_diffraction/` package folder, next to `simulator.py`.

## Package layout

```
pytilt-diffraction/
    pytilt_diffraction/
        __init__.py
        calculator.py   # CIFParser + DiffractionCalculator + rasterize_plane
        simulator.py    # Glazer GUI (perovskite + tilt sliders)
        cif_viewer.py   # generic CIF GUI (any structure)
    streamlit_app.py    # Streamlit: Glazer mode (default page)
    pages/
        2_CIF_viewer.py # Streamlit: CIF upload mode (sidebar nav)
    vendor/
        pytilting/      # upstream Glazer-tilt generator, vendored (GPL v2)
    tests/
        test_physics.py
        test_diffsims_parity.py
    examples/
    docs/
        screenshots/    # PNGs referenced by README
    pyproject.toml
    README.md
    LICENSE
```

## Licensing

This project depends on and redistributes pytilting, which is licensed
under GPL v2. This project is therefore distributed under GPL v2 or later
(see `LICENSE`). Any redistribution must preserve `vendor/pytilting/LICENSE`.

pytilt-diffraction builds on pytilting; its authors ask users to cite:
N. Xie, J. Zhang, S. Raza, N. Zhang, X. Chen and D. Wang, "Generation of
low-symmetry perovskite structures for ab initio computation", J. Phys.:
Condens. Matter 32, 315901 (2020). https://doi.org/10.1088/1361-648X/ab7f6a

## Validation

The `DiffractionCalculator` has been cross-checked against
[diffsims](https://github.com/pyxem/diffsims), pyxem's X-ray and electron
diffraction simulation package. We compute F(hkl) with our Cromer-Mann
parameterisation and compare it with F(hkl) built from diffsims'
`get_kinematical_atomic_scattering_factor` (Doyle-Turner 1968 table) and an
explicit atomic sum, over a grid of 7 materials and 11 hkl:

- Halides: CsPbCl3, CsPbBr3, CsPbI3
- Oxides: SrTiO3, BaTiO3, LaAlO3
- Fluorides: KMgF3

Across the 77 reflections, the difference in |F| is typically under 1% and
always under 5% (with an absolute tolerance of 0.15 for near-zero structure
factors), and the 5 strongest reflections come out in the same order for
every material. The residual difference comes from the two form-factor
parameterisations, not from the summation. Per reflection, our calculator
is also ~3x faster than diffsims' form-factor function with the same
summation (37 us vs 117 us on cubic CsPbBr3).

The tests are in `tests/test_diffsims_parity.py` and `tests/test_physics.py`;
the timings are in `tests/BENCHMARK.md`.

```bash
pip install diffsims pytest
python -m pytest tests/ -q
```

Note on diffsims 0.7.0: its top-level `get_kinematical_structure_factor`
has a bug in `find_asymmetric_positions` that silently drops all but the
first asymmetric-unit atom. We therefore compare the form-factor
parameterisations by summing over atoms ourselves. See
`tests/BENCHMARK.md` for details.
