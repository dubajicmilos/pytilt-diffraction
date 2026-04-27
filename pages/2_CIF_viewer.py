r"""
Streamlit page: generic-CIF single-crystal diffraction viewer.

Browser equivalent of `python -m pytilt_diffraction.cif_viewer`.
The user uploads a CIF; we run the same DiffractionCalculator pipeline
and offer downloads of the hkl table and the 2D rasterised matrix.
"""
from __future__ import annotations

import os
import sys
import io
import tempfile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import streamlit as st

# --- vendored pytilting on sys.path (matches streamlit_app.py) -------------
HERE = os.path.dirname(os.path.abspath(__file__))         # pages/
ROOT = os.path.dirname(HERE)                              # repo root
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'vendor', 'pytilting', 'src'))

from pytilt_diffraction.calculator import (              # noqa: E402
    CIFParser, DiffractionCalculator,
)


# ---------------------------------------------------------------------------
st.set_page_config(page_title='CIF viewer  -  pytilt-diffraction',
                   layout='wide')
st.title('Single-crystal diffraction from any CIF')
st.caption(
    'Upload any single-crystal `.cif` and explore the kinematic diffraction '
    'pattern. Same physics core (`DiffractionCalculator`) as the desktop GUI; '
    'no Glazer / tilt machinery -- the structure is whatever the CIF says.'
)

# ---------------------------------------------------------------------------
# Sidebar: upload + controls
# ---------------------------------------------------------------------------
with st.sidebar:
    st.header('CIF file')
    uploaded = st.file_uploader('Upload .cif', type=['cif'])

    st.divider()
    st.header('Pattern')
    zone_opts = {
        '[001]': (0, 0, 1),
        '[100]': (1, 0, 0),
        '[010]': (0, 1, 0),
        '[110]': (1, 1, 0),
        '[111]': (1, 1, 1),
    }
    zone_label = st.selectbox('Zone axis', list(zone_opts.keys()), index=0)
    zone_axis = zone_opts[zone_label]

    layer = st.slider('HKL layer', -8, 8, 0, step=1)
    d_min = st.slider('d_min (A)', 0.3, 2.5, 0.5, step=0.05)
    h_max = st.slider('h_max', 3, 25, 12, step=1)
    spot_scale = st.slider('spot size', 30.0, 800.0, 250.0)
    log_scale = st.checkbox('log intensity', value=True)
    twin_3 = st.checkbox(
        'twin (3 cubic-parent domains)',
        value=False,
        help=('Incoherent sum over 3 domains: '
              'I = |F(h,k,l)|^2 + |F(l,k,h)|^2 + |F(h,l,k)|^2. '
              'Physical for cubic-parent twinning (e.g. a0a0c+/-).'),
    )
    log_floor_exp = st.slider('log floor (10^x)', -6, -2, -4, step=1)

    st.divider()
    st.header('2D matrix export')
    export_n_pix = st.select_slider(
        'grid (n_pix)', options=[128, 256, 384, 512, 768, 1024], value=512,
    )
    export_extent = st.slider('extent (rlu, +/-)', 1.0, 5.0, 2.5, step=0.5)
    export_sigma = st.slider('Gaussian sigma (rlu)', 0.005, 0.05,
                             0.015, step=0.005)


# ---------------------------------------------------------------------------
# No file -> show a hint and stop
# ---------------------------------------------------------------------------
if uploaded is None:
    st.info('Upload a `.cif` in the sidebar to begin.')
    st.stop()

# ---------------------------------------------------------------------------
# Parse + compute (cached on the file bytes so repeated control changes
# are fast; only the slider-driven steps recompute)
# ---------------------------------------------------------------------------
@st.cache_resource(show_spinner=False)
def _build_calc(file_bytes: bytes, name: str):
    """Write the upload to a temp file (CIFParser opens by path) and
    return the populated calculator."""
    tmp = os.path.join(tempfile.gettempdir(), f'pytilt_streamlit_{name}')
    with open(tmp, 'wb') as f:
        f.write(file_bytes)
    cif = CIFParser(tmp)
    calc = DiffractionCalculator(cif)
    return calc, cif, tmp


calc, cif, tmp_path = _build_calc(uploaded.getvalue(), uploaded.name)

# Reference intensity for cross-layer normalisation (matches desktop GUI)
ref0 = calc.get_plane_reflections(zone_axis, h_max, d_min, I_min=0.0,
                                  layer=0, twin=twin_3)
I_max_ref = max((r['I'] for r in ref0), default=1.0)

refl = calc.get_plane_reflections(
    zone_axis, h_max, d_min,
    I_min=10.0 ** (log_floor_exp - 1),
    layer=layer, I_max_ref=I_max_ref, twin=twin_3,
)
refl = calc.get_2d_coordinates(refl, zone_axis, grid=None)

# ---------------------------------------------------------------------------
# Layout: pattern + info / table / matrix
# ---------------------------------------------------------------------------
col_pattern, col_panel = st.columns([2, 1])

with col_pattern:
    if not refl:
        st.warning('No reflections passed the d_min / I_min filters.')
    else:
        x = np.array([r['x'] for r in refl])
        y = np.array([r['y'] for r in refl])
        I = np.array([r['I_norm'] for r in refl])
        if log_scale:
            floor = 10.0 ** log_floor_exp
            Iv = (np.log10(np.clip(I, floor, 1.0)) - np.log10(floor)) \
                 / (-np.log10(floor))
        else:
            Iv = np.clip(I, 0.0, 1.0) ** 0.4
        Icolor = np.clip(Iv, 0.18, 1.0)
        sizes  = spot_scale * np.clip(Iv, 0.02, 1.0) + 6.0

        fig, ax = plt.subplots(figsize=(6.5, 6.5))
        ax.set_facecolor('#fafafa')
        ax.scatter(x, y, s=sizes, c=Icolor, cmap='Blues',
                   alpha=0.95, edgecolors='#1f2328', linewidths=0.5,
                   vmin=0.0, vmax=1.0, zorder=5)
        ax.scatter([0], [0], s=110, c='#c8102e', edgecolors='#1f2328',
                   linewidths=1.0, zorder=10)
        ax.set_aspect('equal')
        ax.grid(True, color='#d0d7de', linewidth=0.5)
        ax.axhline(0, color='#d0d7de', linewidth=0.8)
        ax.axvline(0, color='#d0d7de', linewidth=0.8)
        ax.set_xlabel('reciprocal lattice (rlu)')
        ax.set_ylabel('reciprocal lattice (rlu)')
        ax.set_title(
            f'zone {zone_label}   L = {layer:+d}   .   '
            f'{uploaded.name}   .   {cif.space_group}'
            + ('   .   3-twin' if twin_3 else ''),
            fontsize=10,
        )
        st.pyplot(fig, clear_figure=True)

with col_panel:
    st.subheader('Structure')
    cp = cif.cell_params
    st.text(
        f"file        {uploaded.name}\n"
        f"space grp   {cif.space_group}\n"
        f"atoms       {len(cif.atoms)}\n"
        f"\n"
        f"a = {cp['a']:8.4f} A\n"
        f"b = {cp['b']:8.4f} A\n"
        f"c = {cp['c']:8.4f} A\n"
        f"alpha = {cp['alpha']:6.2f} deg\n"
        f"beta  = {cp['beta']:6.2f} deg\n"
        f"gamma = {cp['gamma']:6.2f} deg\n"
        f"\n"
        f"reflections   {len(refl)}\n"
        f"d_min (A)     {d_min:.2f}\n"
        f"h_max         {h_max}\n"
        f"layer L       {layer:+d}\n"
        f"twin          {'on' if twin_3 else 'off'}"
    )

# ---------------------------------------------------------------------------
# Strongest reflections + downloads
# ---------------------------------------------------------------------------
st.divider()
top_left, top_right = st.columns([2, 1])

with top_left:
    st.subheader('Strongest reflections')
    if refl:
        top = sorted(refl, key=lambda r: -r['I_norm'])[:30]
        rows = [
            {'h': r['h'], 'k': r['k'], 'l': r['l'],
             'd (A)': round(r['d'], 4),
             '|F|':   round(r['F'], 3),
             'I':     r['I'],
             'I_norm': round(r['I_norm'], 4)}
            for r in top
        ]
        st.dataframe(rows, use_container_width=True, height=320)

        # Build a hkl text export string for download
        buf = io.StringIO()
        buf.write(f"# CIF: {uploaded.name}\n")
        buf.write(f"# zone axis : {zone_axis}\n")
        buf.write(f"# layer L   : {layer:+d}\n")
        buf.write(f"# d_min (A) : {d_min:.3f}\n")
        buf.write(f"# h_max     : {h_max}\n")
        buf.write(f"# twin (3)  : {twin_3}\n")
        buf.write(f"# columns   : h k l d(A) |F| I I_norm\n")
        for r in sorted(refl, key=lambda x: -x['I']):
            buf.write(f"{r['h']:4d} {r['k']:4d} {r['l']:4d} "
                      f"{r['d']:8.4f} {r['F']:10.3f} "
                      f"{r['I']:12.3e} {r['I_norm']:8.4f}\n")
        st.download_button(
            'Download hkl table (.txt)',
            data=buf.getvalue().encode('utf-8'),
            file_name=f'reflections_{os.path.splitext(uploaded.name)[0]}'
                      f'_zone{zone_label.strip("[]")}_L{layer:+d}.txt'
                      .replace('+', 'p').replace('-', 'm'),
            mime='text/plain',
        )

with top_right:
    st.subheader('2D matrix export')
    st.caption(
        'Detector-style rasterisation: each reflection becomes a Gaussian '
        'of width sigma in rlu, summed onto a regular grid. Use for '
        'side-by-side comparison with experimental detector images.'
    )
    if refl:
        ext = float(export_extent)
        img = calc.rasterize_plane(
            refl, extent=(-ext, ext, -ext, ext),
            n_pix=int(export_n_pix), sigma=float(export_sigma),
            use_norm=True,
        )
        fig2, ax2 = plt.subplots(figsize=(4.5, 4.5))
        if img.max() > 0:
            ax2.imshow(img, origin='lower',
                       extent=(-ext, ext, -ext, ext),
                       cmap='Blues',
                       norm=LogNorm(vmin=max(img.max() * 1e-4, 1e-12),
                                    vmax=img.max()),
                       interpolation='bilinear', aspect='equal')
        ax2.set_xlabel('rlu')
        ax2.set_ylabel('rlu')
        ax2.set_title(f'rasterised, {img.shape[0]}x{img.shape[1]}',
                      fontsize=9)
        st.pyplot(fig2, clear_figure=True)

        # NPY download
        npy_buf = io.BytesIO()
        np.save(npy_buf, img)
        st.download_button(
            'Download matrix (.npy)',
            data=npy_buf.getvalue(),
            file_name=f'matrix_{os.path.splitext(uploaded.name)[0]}'
                      f'_zone{zone_label.strip("[]")}_L{layer:+d}.npy'
                      .replace('+', 'p').replace('-', 'm'),
            mime='application/octet-stream',
        )

        # CSV download
        csv_buf = io.StringIO()
        np.savetxt(csv_buf, img, delimiter=',', fmt='%.6e')
        st.download_button(
            'Download matrix (.csv)',
            data=csv_buf.getvalue().encode('utf-8'),
            file_name=f'matrix_{os.path.splitext(uploaded.name)[0]}'
                      f'_zone{zone_label.strip("[]")}_L{layer:+d}.csv'
                      .replace('+', 'p').replace('-', 'm'),
            mime='text/csv',
        )
