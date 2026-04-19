"""
Streamlit front-end for pytilt-diffraction.

Reuses the physics core:
    - `pytilting.Distortion` (vendored)   -> applies Glazer tilts
    - `pytilt_diffraction.calculator.DiffractionCalculator`
          -> structure factor + powder pattern

Only the UI layer is replaced: matplotlib widgets are swapped for Streamlit
controls in the sidebar, and the figures are rendered with `st.pyplot`.

Run locally:
    pip install -r requirements.txt
    streamlit run streamlit_app.py

Deploy: push to GitHub, connect the repo at https://share.streamlit.io.
"""

import math
import os
import sys
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

# Make the vendored pytilting package importable (its modules use top-level
# imports like `from puc import Puc`, so the `src/` folder itself must be on
# sys.path rather than being treated as a subpackage).
HERE = os.path.dirname(os.path.abspath(__file__))
_VENDOR_SRC = os.path.normpath(os.path.join(HERE, "vendor", "pytilting", "src"))
if _VENDOR_SRC not in sys.path:
    sys.path.insert(0, _VENDOR_SRC)

from distortion import Distortion  # noqa: E402

from pytilt_diffraction.calculator import CIFParser, DiffractionCalculator  # noqa: E402
from pytilt_diffraction.simulator import (  # noqa: E402
    GLAZER_PRESETS,
    MATERIALS,
    _fmt_hkl,
    glazer_equality_groups,
)


# ---------------------------------------------------------------------------
# Page setup
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="pytilt-diffraction",
    page_icon=":atom_symbol:",
    layout="wide",
)

st.title("Glazer-tilt diffraction simulator")
st.caption(
    "Pick an ABX3 perovskite, choose a Howard-Stokes tilt system, and see "
    "the single-crystal zone slice and powder pattern update live."
)


# ---------------------------------------------------------------------------
# Core pipeline, cached on the subset of parameters that affects the crystal.
# ---------------------------------------------------------------------------
@st.cache_resource(show_spinner=False)
def build_calculator(symbols_tuple, a0, glazer, omega_tuple, grid=(2, 2, 2)):
    """Tilt the supercell, write a CIF, return a ready DiffractionCalculator.

    Cached on (composition, a0, glazer, omega) so that scanning unrelated
    widgets (h_max, d_min, spot size, wavelength, ...) does not rebuild the
    structure.
    """
    distortion = Distortion(
        system={
            "symbols": list(symbols_tuple),
            "lattice_constant": float(a0),
            "grid": tuple(grid),
            "covera": 1.0,
        }
    )
    distortion.distort = {
        "glazer": glazer,
        "omega": tuple(math.radians(abs(w)) for w in omega_tuple),
        "u": (0.0, 0.0, 0.0),
        "k_u": 2 * math.pi * np.zeros((3, 3)),
        "local_mode": [0.0] * 5,
        "modes": [],
    }
    atoms = distortion.get_atoms()
    cif_path = os.path.join(
        tempfile.gettempdir(),
        f"pytilt_streamlit_{abs(hash((symbols_tuple, a0, glazer, omega_tuple)))}.cif",
    )
    atoms.write(cif_path)
    return DiffractionCalculator(CIFParser(cif_path))


def apply_glazer_consistency(glazer, omega_deg):
    """Enforce equal-letter -> equal-magnitude and '0' -> 0. Same rule the
    matplotlib simulator uses; duplicated here so we don't need a GUI
    object."""
    groups = glazer_equality_groups(glazer)
    signs = [glazer[1], glazer[3], glazer[5]]
    mags = [abs(w) for w in omega_deg]
    for grp in groups:
        master = mags[grp[0]]
        for i in grp:
            mags[i] = master
    for i in range(3):
        if signs[i] == "0":
            mags[i] = 0.0
    return mags


# ---------------------------------------------------------------------------
# Sidebar controls
# ---------------------------------------------------------------------------
with st.sidebar:
    st.header("Structure")

    mat_labels = [f"{m[0]}  (a0 = {m[4]:.3f} A)" for m in MATERIALS]
    mat_idx = st.selectbox(
        "Material preset",
        range(len(MATERIALS)),
        index=1,  # CsPbBr3
        format_func=lambda i: mat_labels[i],
    )
    _, A, B, X, a0_default = MATERIALS[mat_idx]

    a0 = st.slider(
        "Lattice constant a0 (A)", 3.0, 7.5, float(a0_default), step=0.001
    )

    glazer_labels = [f"{g}  ({sg})  -  {note}" for g, sg, note in GLAZER_PRESETS]
    glazer_idx = st.selectbox(
        "Glazer tilt system",
        range(len(GLAZER_PRESETS)),
        index=0,
        format_func=lambda i: glazer_labels[i],
    )
    glazer = GLAZER_PRESETS[glazer_idx][0]

    st.caption(
        "Tilt magnitudes. Equal-letter axes are tied automatically; '0' "
        "axes are forced to 0."
    )
    omega_x = st.slider("omega_x (deg)", 0.0, 20.0, 6.0, step=0.1)
    omega_y = st.slider("omega_y (deg)", 0.0, 20.0, 6.0, step=0.1)
    omega_z = st.slider("omega_z (deg)", 0.0, 20.0, 6.0, step=0.1)
    omega_deg = apply_glazer_consistency(glazer, [omega_x, omega_y, omega_z])

    st.divider()
    st.header("Single-crystal view")
    zone_opts = {
        "[001]  hk0": (0, 0, 1),
        "[100]  0kl": (1, 0, 0),
        "[010]  h0l": (0, 1, 0),
        "[110]  hh-l": (1, 1, 0),
        "[111]": (1, 1, 1),
    }
    zone_label = st.selectbox("Zone axis", list(zone_opts.keys()), index=0)
    zone_axis = zone_opts[zone_label]

    layer = st.slider("HK_L layer (supercell)", -4, 4, 0, step=1)
    d_min = st.slider("d_min (A)", 0.3, 2.5, 0.8, step=0.05)
    h_max_sc = st.slider("h_max (single crystal)", 3, 14, 8, step=1)
    spot_scale = st.slider("spot size", 30.0, 800.0, 250.0)
    show_labels = st.checkbox("show hkl labels", value=True)
    log_scale = st.checkbox("log intensity", value=True)
    log_floor_exp = st.slider("log floor (10^x)", -6, -2, -4, step=1)

    st.divider()
    st.header("Powder pattern")
    sources = {
        "Cu Ka1  1.5406 A": 1.5406,
        "Mo Ka1  0.7093 A": 0.7093,
        "Co Ka1  1.7890 A": 1.7890,
        "Cr Ka1  2.2897 A": 2.2897,
        "Ag Ka1  0.5594 A": 0.5594,
    }
    src_label = st.selectbox("Source", list(sources.keys()), index=0)
    wavelength = sources[src_label]
    two_theta_max = st.slider("2theta max (deg)", 20.0, 140.0, 90.0)
    h_max_pow = st.slider("h_max (powder)", 4, 16, 10, step=1)
    fwhm = st.slider("peak FWHM (deg)", 0.02, 1.0, 0.10, step=0.01)


# ---------------------------------------------------------------------------
# Build the calculator once per distinct structure and run both scans.
# ---------------------------------------------------------------------------
symbols_tuple = (A, B, X)
omega_tuple = tuple(round(w, 4) for w in omega_deg)
calc = build_calculator(symbols_tuple, round(a0, 6), glazer, omega_tuple)

# Global intensity normaliser (strongest layer-0 peak), matches the desktop
# simulator so that odd-layer superlattice peaks don't get rescaled to 1.
ref_layer0 = calc.get_plane_reflections(zone_axis, h_max_sc, d_min, I_min=0.0, layer=0)
I_max_ref = max((r["I"] for r in ref_layer0), default=1.0)

refl = calc.get_plane_reflections(
    zone_axis,
    h_max_sc,
    d_min,
    I_min=10.0 ** (log_floor_exp - 1),
    layer=layer,
    I_max_ref=I_max_ref,
)
refl = calc.get_2d_coordinates(refl, zone_axis, grid=(2, 2, 2))

peaks = calc.powder_pattern(
    wavelength=wavelength,
    two_theta_max=two_theta_max,
    h_max=h_max_pow,
    I_min=1e-4,
)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
BG, PANEL, TEXT, TEXT_DIM = "#0d1117", "#161b22", "#e6edf3", "#8b949e"
ACCENT, ACCENT_2, GRID = "#ff6b35", "#4a90e2", "#21262d"


def draw_zone(refl, zone_axis, layer, spot_scale, show_labels,
              log_scale=True, log_floor=1e-4):
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor("#000000")
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=TEXT_DIM)
    ax.set_aspect("equal")
    ax.grid(alpha=0.15, color=GRID)

    if not refl:
        ax.text(
            0.5,
            0.5,
            "(no reflections on this layer)",
            color=TEXT_DIM,
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
    else:
        xs = np.array([r["x"] for r in refl])
        ys = np.array([r["y"] for r in refl])
        Is = np.array([r["I_norm"] for r in refl])
        if log_scale:
            floor = max(log_floor, 1e-10)
            Iplot = (np.log10(np.clip(Is, floor, 1.0)) - np.log10(floor)) / (-np.log10(floor))
            cbar_label = f"log10(I / I_max), floor = 10^{int(np.log10(floor))}"
        else:
            Iplot = Is
            cbar_label = "I / I_max"
        sizes = spot_scale * np.sqrt(np.clip(Iplot, 0.04, 1.0)) + 6
        sc = ax.scatter(xs, ys, s=sizes, c=Iplot, cmap="hot", alpha=0.95,
                        vmin=0, vmax=1)
        plt.colorbar(sc, ax=ax, label=cbar_label, shrink=0.8)
        if show_labels:
            thr = 0.05
            for r, I in zip(refl, Is):
                if I < thr:
                    continue
                ax.annotate(
                    f"({_fmt_hkl(r['h_p'])} {_fmt_hkl(r['k_p'])} {_fmt_hkl(r['l_p'])})",
                    xy=(r["x"], r["y"]),
                    xytext=(4, 4),
                    textcoords="offset points",
                    color=TEXT,
                    fontsize=8,
                )

    u, v, w = zone_axis
    ax.set_title(
        f"Zone [{u}{v}{w}]   layer = {layer}   ({len(refl)} reflections)",
        color=TEXT,
    )
    ax.set_xlabel("reciprocal lattice (parent pseudocubic)", color=TEXT)
    ax.set_ylabel("reciprocal lattice (parent pseudocubic)", color=TEXT)
    return fig


def draw_powder(peaks, wavelength, fwhm, two_theta_max, n_labels=8):
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(9.5, 4.5))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor("#000000")
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=TEXT_DIM)
    ax.set_xlim(0, two_theta_max)
    ax.set_ylim(0, 1.08)
    ax.grid(alpha=0.15, color=GRID)
    ax.set_xlabel("2theta (deg)", color=TEXT)
    ax.set_ylabel("Intensity (a.u.)", color=TEXT)

    if not peaks:
        ax.text(
            0.5,
            0.5,
            "(no peaks in range)",
            color=TEXT_DIM,
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return fig

    two_thetas = np.array([pk["two_theta"] for pk in peaks])
    I_norms = np.array([pk["I_norm"] for pk in peaks])

    # Gaussian-convolved profile
    x = np.linspace(0.0, two_theta_max, 3000)
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    two_sig2 = 2.0 * sigma * sigma
    y = np.zeros_like(x)
    for tt, I in zip(two_thetas, I_norms):
        y += I * np.exp(-((x - tt) ** 2) / two_sig2)
    y /= y.max() if y.max() > 0 else 1.0

    ax.vlines(two_thetas, 0, I_norms, color=ACCENT_2, alpha=0.6, lw=1.0)
    ax.plot(x, y, color=ACCENT, lw=1.3)

    order = np.argsort(-I_norms)
    labelled = set()
    drawn = 0
    for idx in order:
        if drawn >= n_labels:
            break
        key = round(two_thetas[idx], 2)
        if key in labelled:
            continue
        labelled.add(key)
        h, k, l = peaks[idx]["hkl"]
        ax.annotate(
            f"({h}{k}{l})",
            xy=(two_thetas[idx], I_norms[idx]),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            color=TEXT,
            fontsize=8,
        )
        drawn += 1

    ax.set_title(
        f"Powder pattern   lambda = {wavelength:.4f} A   sigma = {sigma:.3f} deg",
        color=TEXT,
    )
    return fig


# ---------------------------------------------------------------------------
# Top banner (composition + tilt readout)
# ---------------------------------------------------------------------------
sg = GLAZER_PRESETS[glazer_idx][1]
banner = (
    f"**{A}{B}{X}3**  |  a0 = {a0:.4f} A  |  Glazer **{glazer}** ({sg})  |  "
    f"omega = ({omega_deg[0]:.1f}, {omega_deg[1]:.1f}, {omega_deg[2]:.1f}) deg"
)
st.markdown(banner)


col_a, col_b = st.columns([1, 1.3])
with col_a:
    st.subheader("Reciprocal-space zone slice")
    st.pyplot(
        draw_zone(refl, zone_axis, layer, spot_scale, show_labels,
                  log_scale=log_scale, log_floor=10.0 ** log_floor_exp),
        clear_figure=True,
    )
with col_b:
    st.subheader("Powder pattern")
    st.pyplot(
        draw_powder(peaks, wavelength, fwhm, two_theta_max),
        clear_figure=True,
    )


# ---------------------------------------------------------------------------
# Tables + downloads
# ---------------------------------------------------------------------------
tab_zone, tab_pow = st.tabs(["Single-crystal reflections", "Powder peaks"])

with tab_zone:
    if refl:
        rows = [
            {
                "h": round(r["h_p"], 4),
                "k": round(r["k_p"], 4),
                "l": round(r["l_p"], 4),
                "d (A)": round(r["d"], 4),
                "I / I_max": round(r["I_norm"], 4),
            }
            for r in sorted(refl, key=lambda r: -r["I_norm"])
        ]
        st.dataframe(rows, use_container_width=True, hide_index=True)
        csv = "h,k,l,d_A,I_norm\n" + "\n".join(
            f"{r['h']},{r['k']},{r['l']},{r['d (A)']},{r['I / I_max']}" for r in rows
        )
        st.download_button(
            "Download reflection CSV",
            csv,
            file_name=f"{A}{B}{X}3_{glazer}_zone{''.join(map(str, zone_axis))}_L{layer}.csv",
            mime="text/csv",
        )
    else:
        st.info("No reflections on this layer at the current cutoffs.")

with tab_pow:
    if peaks:
        rows = [
            {
                "h": pk["hkl"][0],
                "k": pk["hkl"][1],
                "l": pk["hkl"][2],
                "2theta (deg)": round(pk["two_theta"], 3),
                "d (A)": round(pk["d"], 4),
                "I / I_max": round(pk["I_norm"], 4),
                "multiplicity": pk["mult"],
            }
            for pk in sorted(peaks, key=lambda pk: pk["two_theta"])
        ]
        st.dataframe(rows, use_container_width=True, hide_index=True)
        csv = (
            "h,k,l,two_theta_deg,d_A,I_norm,multiplicity\n"
            + "\n".join(
                f"{r['h']},{r['k']},{r['l']},{r['2theta (deg)']},{r['d (A)']},"
                f"{r['I / I_max']},{r['multiplicity']}"
                for r in rows
            )
        )
        st.download_button(
            "Download powder CSV",
            csv,
            file_name=f"{A}{B}{X}3_{glazer}_powder.csv",
            mime="text/csv",
        )
    else:
        st.info("No powder peaks in the selected 2theta range.")


with st.expander("About"):
    st.markdown(
        """
This web app reuses the same physics core as the desktop
`pytilt-diffraction` simulator: the vendored
[pytilting](https://gitlab.com/pyseries/pytilting) Glazer-tilt generator
(GPL v2) drives the `Distortion` object, which emits a tilted ABX3
supercell. Our `DiffractionCalculator` then computes structure factors
with Cromer-Mann X-ray form factors and bins them into a kinematic
powder pattern (Lorentz-polarisation applied).

The `F(hkl)` engine is cross-checked against
[diffsims](https://github.com/pyxem/diffsims) on 7 perovskites (77
reflections): agreement is typically <1 %, always <5 %. See
`tests/test_diffsims_parity.py` and `tests/BENCHMARK.md` in the
[GitHub repo](https://github.com/dubajicmilos/pytilt-diffraction).
"""
    )
