r"""
Generic single-crystal CIF diffraction viewer.

Loads any single-crystal CIF and lets you browse the kinematic
diffraction pattern interactively:

  - radio buttons for zone axis  ([001], [100], [010], [110], [111])
  - slider for HKL layer index (integer L in the CIF's reciprocal lattice;
    no parent / supercell rescaling -- this is "what the CIF says")
  - sliders for d_min, h_max, spot size, label threshold
  - log / linear intensity, hkl labels, twin-3 (cubic 90 deg domains)
  - Save PNG, Export hkl table (.txt), Export 2D matrix (.npy + .csv)

The visual encoding (Blues cmap, intensity-mapped color + size, light plot
panel, dark surrounding GUI) matches the Glazer simulator, but all the
tilt / Glazer / material-preset machinery is dropped -- it's irrelevant
for a generic CIF.

Run:
    python -m pytilt_diffraction.cif_viewer path/to/structure.cif

If no path is given, a Tk file-picker is shown.
"""
from __future__ import annotations

import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons, Slider, CheckButtons, Button

from .calculator import CIFParser, DiffractionCalculator
from .simulator import TiltDiffractionSimulator, _fmt_hkl

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# Helper: fall back to a Tk file picker if no path was given on the CLI
# ---------------------------------------------------------------------------
def _pick_cif_dialog():
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception:
        return None
    root = tk.Tk(); root.withdraw()
    path = filedialog.askopenfilename(
        title='Open CIF',
        filetypes=[('CIF files', '*.cif'), ('All files', '*.*')],
    )
    root.destroy()
    return path or None


# ---------------------------------------------------------------------------
# Viewer
# ---------------------------------------------------------------------------
class CIFViewer:
    """
    Single-crystal kinematic diffraction viewer for any CIF.

    Reuses the rendering primitives from `TiltDiffractionSimulator`
    (`_draw_pattern`, `_draw_info`, `_draw_reflections`) by sharing the
    same attribute names: self.ax, self.reflections, self.log_scale,
    self.twin_3, self.spot_scale, self.show_labels, self.label_thresh,
    self.I_floor, plus the same color constants.
    """

    # Inherit the GUI palette from TiltDiffractionSimulator so the look
    # is consistent (light plot panel, dark surrounding GUI).
    BG               = TiltDiffractionSimulator.BG
    PANEL            = TiltDiffractionSimulator.PANEL
    PLOT_BG          = TiltDiffractionSimulator.PLOT_BG
    ACCENT           = TiltDiffractionSimulator.ACCENT
    ACCENT_2         = TiltDiffractionSimulator.ACCENT_2
    TEXT             = TiltDiffractionSimulator.TEXT
    TEXT_DIM         = TiltDiffractionSimulator.TEXT_DIM
    GRID             = TiltDiffractionSimulator.GRID
    PLOT_TEXT        = TiltDiffractionSimulator.PLOT_TEXT
    PLOT_TEXT_DIM    = TiltDiffractionSimulator.PLOT_TEXT_DIM
    PLOT_GRID        = TiltDiffractionSimulator.PLOT_GRID
    SPOT_EDGE        = TiltDiffractionSimulator.SPOT_EDGE
    ZERO_SPOT_FACE   = TiltDiffractionSimulator.ZERO_SPOT_FACE
    SPOT_CMAP        = TiltDiffractionSimulator.SPOT_CMAP

    ZONE_OPTIONS = [
        ('[001]', (0, 0, 1)),
        ('[100]', (1, 0, 0)),
        ('[010]', (0, 1, 0)),
        ('[110]', (1, 1, 0)),
        ('[111]', (1, 1, 1)),
    ]

    def __init__(self, cif_path):
        self.cif_path = os.path.abspath(cif_path)
        self.cif  = CIFParser(self.cif_path)
        self.calc = DiffractionCalculator(self.cif)

        # Pattern parameters
        self.zone_axis  = (0, 0, 1)
        self.h_max      = 12
        self.d_min      = 0.5
        self.I_min      = 1e-5
        self.spot_scale = 250.0
        self.show_labels   = True
        self.label_thresh  = 0.05
        self.log_scale     = True
        self.I_floor       = 1e-4
        self.layer         = 0
        self.twin_3        = False

        # Matrix-export parameters
        self.export_n_pix = 512
        self.export_sigma_rlu = 0.015
        self.export_extent_rlu = 2.5     # +/- this in rlu

        # Glazer-only attributes that _draw_pattern peeks at -- give them
        # neutral values so the title still renders cleanly.
        self.glazer = ''                 # blank: not a Glazer system
        self.omega_deg = [0.0, 0.0, 0.0]
        self.grid = (1, 1, 1)            # no supercell rescaling
        self.symbols = list(self.cif.atom_types) if hasattr(
            self.cif, 'atom_types') else []

        # Internal state
        self.reflections = []
        self.I_max_ref   = 1.0
        self._suppress   = False

    # ------------------------------------------------------------------
    # Pipeline (matches TiltDiffractionSimulator's API surface)
    # ------------------------------------------------------------------
    def _update_I_max_ref(self):
        ref = self.calc.get_plane_reflections(
            self.zone_axis, self.h_max, self.d_min, I_min=0.0, layer=0,
            twin=self.twin_3,
        )
        self.I_max_ref = max((r['I'] for r in ref), default=1.0)

    def recompute_pattern(self):
        self.reflections = self.calc.get_plane_reflections(
            self.zone_axis, self.h_max, self.d_min, self.I_min,
            layer=self.layer, I_max_ref=self.I_max_ref, twin=self.twin_3,
        )
        # CIF viewer: no supercell rescaling. Pass grid=None so h_p/k_p/l_p
        # equal h/k/l (genuine "what the CIF says" rlu).
        self.reflections = self.calc.get_2d_coordinates(
            self.reflections, self.zone_axis, grid=None,
        )

    def _parent_layer(self):
        # Used by _draw_pattern's title. With grid=None the "parent" layer
        # equals the supercell layer -- they're the same coordinate system.
        return float(self.layer)

    def _composition_banner(self):
        return f"CIF: {os.path.basename(self.cif_path)}"

    # ------------------------------------------------------------------
    # Visual rendering. The scatter / cmap / sizing logic comes for free
    # from TiltDiffractionSimulator._draw_pattern; we then override only
    # the title and axis labels (no Glazer / supercell concepts here).
    # _draw_reflections is reused as-is because its output is generic.
    # _draw_info is rewritten because the Glazer version assumes
    # ABX3 + a Glazer tilt pattern.
    # ------------------------------------------------------------------
    _draw_reflections = TiltDiffractionSimulator._draw_reflections

    def _draw_pattern(self):
        # Reuse the borrowed rendering, then replace the title and axis
        # labels with CIF-appropriate text.
        TiltDiffractionSimulator._draw_pattern(self)
        za = ''.join(str(i) for i in self.zone_axis)
        sg = getattr(self.cif, 'space_group', '?') or '?'
        self.ax.set_title(
            f"zone [{za}]   {self._zone_law_str()}   .   "
            f"{os.path.basename(self.cif_path)}   .   {sg}"
            + ('   .   3-twin domains' if self.twin_3 else ''),
            color=self.PLOT_TEXT, fontsize=11, pad=10,
        )
        self.ax.set_xlabel('reciprocal lattice (rlu)',
                           color=self.PLOT_TEXT_DIM)
        self.ax.set_ylabel('reciprocal lattice (rlu)',
                           color=self.PLOT_TEXT_DIM)

    def _zone_law_str(self):
        """Human-readable form of the zone-law constraint h*u + k*v + l*w = L
        currently in effect, e.g.  'l = +2'  for zone [001] layer 2,
        'h + k = +1'  for zone [110] layer 1,  'h + k + l = -3'  for
        zone [111] layer -3."""
        u, v, w = self.zone_axis
        parts = []
        for c, sym in zip((u, v, w), ('h', 'k', 'l')):
            if c == 0:
                continue
            if c == 1:
                parts.append(sym)
            elif c == -1:
                parts.append('-' + sym)
            else:
                parts.append(f'{c}{sym}')
        return ' + '.join(parts) + f' = {self.layer:+d}'

    def _draw_info(self):
        self.info_ax.clear()
        self.info_ax.axis('off')
        self.info_ax.set_facecolor(self.PANEL)

        cp = self.cif.cell_params
        sg = getattr(self.cif, 'space_group', '?') or '?'
        n_atoms = len(getattr(self.cif, 'atoms', []) or [])
        n_refl = len(self.reflections)

        # Build the unique-element list from whatever shape the CIF parser
        # exposes. Fall back gracefully if a field is missing.
        elems = []
        for atom in getattr(self.cif, 'atoms', []) or []:
            sym = atom.get('element') or atom.get('label') or ''
            sym = ''.join(c for c in sym if c.isalpha())
            if sym and sym not in elems:
                elems.append(sym)
        formula = ' '.join(elems) if elems else '?'

        info = (
            f"CIF\n"
            f"---\n"
            f"file       {os.path.basename(self.cif_path)}\n"
            f"space grp  {sg}\n"
            f"atoms      {n_atoms}\n"
            f"elements   {formula}\n"
            f"\n"
            f"UNIT CELL\n"
            f"---------\n"
            f"a = {cp['a']:8.4f} A\n"
            f"b = {cp['b']:8.4f} A\n"
            f"c = {cp['c']:8.4f} A\n"
            f"alpha = {cp['alpha']:6.2f} deg\n"
            f"beta  = {cp['beta']:6.2f} deg\n"
            f"gamma = {cp['gamma']:6.2f} deg\n"
            f"\n"
            f"PATTERN\n"
            f"-------\n"
            f"zone axis    {self.zone_axis}\n"
            f"layer L      {self.layer:+d}\n"
            f"d_min (A)    {self.d_min:5.2f}\n"
            f"h_max        {self.h_max}\n"
            f"reflections  {n_refl}\n"
            f"twin (3 dom) {'on' if self.twin_3 else 'off'}\n"
        )
        self.info_ax.text(
            0.05, 0.97, info, transform=self.info_ax.transAxes,
            fontsize=8.5, color=self.TEXT, family='monospace',
            va='top',
        )

    # ------------------------------------------------------------------
    # GUI
    # ------------------------------------------------------------------
    def run(self):
        self._update_I_max_ref()
        self.recompute_pattern()

        plt.style.use('dark_background')
        self.fig = plt.figure(figsize=(15, 9))
        self.fig.canvas.manager.set_window_title(
            f'pytilt-diffraction  .  CIF viewer  .  {os.path.basename(self.cif_path)}'
        )
        self.fig.patch.set_facecolor(self.BG)

        # Main reciprocal-space plot
        self.ax = self.fig.add_axes([0.05, 0.30, 0.55, 0.62])
        self.ax.set_facecolor(self.PLOT_BG)

        # Info / readout panels
        self.info_ax = self.fig.add_axes([0.62, 0.56, 0.17, 0.36])
        self.info_ax.set_facecolor(self.PANEL)
        self.info_ax.set_xticks([]); self.info_ax.set_yticks([])
        for s in self.info_ax.spines.values():
            s.set_color(self.GRID)

        self.refl_ax = self.fig.add_axes([0.81, 0.56, 0.17, 0.36])
        self.refl_ax.set_facecolor(self.PANEL)
        self.refl_ax.set_xticks([]); self.refl_ax.set_yticks([])
        for s in self.refl_ax.spines.values():
            s.set_color(self.GRID)

        # Title banner
        self.fig.text(0.05, 0.955,
                      'Single-Crystal Diffraction (any CIF)',
                      color=self.TEXT, fontsize=15, fontweight='bold')
        self.composition_txt = self.fig.text(
            0.05, 0.928, self._composition_banner(),
            color=self.TEXT_DIM, fontsize=10,
        )
        self.status_txt = self.fig.text(
            0.62, 0.955, '', color=self.ACCENT,
            fontsize=9, va='center',
        )

        self._build_controls()
        self._draw_pattern()
        self._draw_info()
        self._draw_reflections()

        plt.show()

    def _build_controls(self):
        # Zone axis radio
        ax_zlbl = self.fig.add_axes([0.62, 0.46, 0.18, 0.03])
        ax_zlbl.axis('off'); ax_zlbl.set_facecolor(self.PANEL)
        ax_zlbl.text(0, 0.5, 'Zone axis (view direction)',
                     color=self.TEXT, fontsize=10, fontweight='bold',
                     va='center')

        ax_zone = self.fig.add_axes([0.62, 0.27, 0.18, 0.18])
        ax_zone.set_facecolor(self.PANEL)
        labels = [opt[0] for opt in self.ZONE_OPTIONS]
        active = labels.index('[001]')
        self.radio_zone = RadioButtons(
            ax_zone, labels, active=active, activecolor=self.ACCENT_2,
        )
        for lbl in self.radio_zone.labels:
            lbl.set_color(self.TEXT); lbl.set_fontsize(9)
        for sp in self.radio_zone.ax.spines.values():
            sp.set_color(self.GRID)
        self.radio_zone.on_clicked(self._on_zone)

        # Sliders (left, below plot)
        sl_color = self.ACCENT_2

        ax_layer = self.fig.add_axes([0.08, 0.21, 0.45, 0.020])
        self.s_layer = Slider(ax_layer, 'HKL layer', -8, 8,
                              valinit=self.layer, valstep=1, color=sl_color)
        ax_dmin = self.fig.add_axes([0.08, 0.17, 0.45, 0.020])
        self.s_dmin = Slider(ax_dmin, 'd_min (A)', 0.3, 2.5,
                             valinit=self.d_min, color=sl_color)
        ax_hmax = self.fig.add_axes([0.08, 0.13, 0.45, 0.020])
        self.s_hmax = Slider(ax_hmax, 'h_max', 3, 25,
                             valinit=self.h_max, valstep=1, color=sl_color)
        ax_spot = self.fig.add_axes([0.08, 0.09, 0.45, 0.020])
        self.s_spot = Slider(ax_spot, 'spot size', 30.0, 800.0,
                             valinit=self.spot_scale, color=sl_color)
        ax_lthr = self.fig.add_axes([0.08, 0.05, 0.45, 0.020])
        self.s_lthr = Slider(ax_lthr, 'label I/I_max >=', 0.0, 0.5,
                             valinit=self.label_thresh, color=sl_color)
        ax_npix = self.fig.add_axes([0.08, 0.01, 0.45, 0.020])
        self.s_npix = Slider(ax_npix, 'export grid (n_pix)', 64, 1024,
                             valinit=self.export_n_pix, valstep=64,
                             color=self.ACCENT)

        for s in (self.s_layer, self.s_dmin, self.s_hmax, self.s_spot,
                  self.s_lthr, self.s_npix):
            s.label.set_color(self.TEXT)
            s.valtext.set_color(self.TEXT)

        self.s_layer.on_changed(self._on_layer)
        self.s_dmin .on_changed(self._on_dmin)
        self.s_hmax .on_changed(self._on_hmax)
        self.s_spot .on_changed(self._on_spot)
        self.s_lthr .on_changed(self._on_lthr)
        self.s_npix .on_changed(self._on_npix)

        # Check boxes
        ax_chk = self.fig.add_axes([0.81, 0.27, 0.17, 0.085])
        ax_chk.set_facecolor(self.PANEL)
        self.check = CheckButtons(
            ax_chk,
            ['Show hkl labels', 'log intensity', 'twin (3 domains)'],
            [self.show_labels, self.log_scale, self.twin_3],
        )
        for lbl in self.check.labels:
            lbl.set_color(self.TEXT); lbl.set_fontsize(8)
        self.check.on_clicked(self._on_check)

        # Buttons
        ax_save = self.fig.add_axes([0.62, 0.20, 0.17, 0.040])
        self.btn_save = Button(ax_save, 'Save PNG',
                               color=self.ACCENT, hovercolor='#ff8c5c')
        self.btn_save.label.set_color('white')
        self.btn_save.on_clicked(self._on_save)

        ax_exhkl = self.fig.add_axes([0.81, 0.20, 0.17, 0.040])
        self.btn_exhkl = Button(ax_exhkl, 'Export hkl',
                                color=self.ACCENT_2, hovercolor='#6ba7e8')
        self.btn_exhkl.label.set_color('white')
        self.btn_exhkl.on_clicked(self._on_export_hkl)

        ax_exmtx = self.fig.add_axes([0.62, 0.15, 0.17, 0.040])
        self.btn_exmtx = Button(ax_exmtx, 'Export 2D matrix',
                                color='#8b5cf6', hovercolor='#a78bfa')
        self.btn_exmtx.label.set_color('white')
        self.btn_exmtx.on_clicked(self._on_export_matrix)

        ax_loadcif = self.fig.add_axes([0.81, 0.15, 0.17, 0.040])
        self.btn_load = Button(ax_loadcif, 'Load CIF...',
                               color=self.GRID, hovercolor='#30363d')
        self.btn_load.label.set_color(self.TEXT)
        self.btn_load.on_clicked(self._on_load_cif)

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------
    def _redraw(self, full=True):
        if full:
            self._update_I_max_ref()
            self.recompute_pattern()
        self._draw_pattern()
        self._draw_info()
        self._draw_reflections()
        self.fig.canvas.draw_idle()

    def _on_zone(self, label):
        for lab, axis in self.ZONE_OPTIONS:
            if lab == label:
                self.zone_axis = axis
                break
        self._redraw(full=True)

    def _on_layer(self, val):
        self.layer = int(val)
        self._redraw(full=True)

    def _on_dmin(self, val):
        self.d_min = float(val)
        self._redraw(full=True)

    def _on_hmax(self, val):
        self.h_max = int(val)
        self._redraw(full=True)

    def _on_spot(self, val):
        self.spot_scale = float(val)
        self._redraw(full=False)

    def _on_lthr(self, val):
        self.label_thresh = float(val)
        self._redraw(full=False)

    def _on_npix(self, val):
        self.export_n_pix = int(val)
        # No redraw needed; only affects matrix export.

    def _on_check(self, label):
        if 'labels' in label:
            self.show_labels = not self.show_labels
            self._redraw(full=False)
        elif 'log' in label:
            self.log_scale = not self.log_scale
            self._redraw(full=False)
        elif 'twin' in label:
            self.twin_3 = not self.twin_3
            self._redraw(full=True)

    def _set_status(self, msg, ok=True):
        self.status_txt.set_text(msg)
        self.status_txt.set_color(self.ACCENT if ok else '#f85149')
        self.fig.canvas.draw_idle()

    def _stub(self):
        return os.path.splitext(os.path.basename(self.cif_path))[0]

    def _zonetag(self):
        return ''.join(str(i) for i in self.zone_axis)

    def _on_save(self, _):
        fname = os.path.join(
            os.path.dirname(self.cif_path),
            f'diffraction_{self._stub()}_zone{self._zonetag()}'
            f'_L{self.layer:+d}.png'.replace('+', 'p').replace('-', 'm'),
        )
        try:
            self.fig.savefig(fname, dpi=200, facecolor=self.BG,
                             bbox_inches='tight')
            self._set_status(f"Saved PNG: {os.path.basename(fname)}", ok=True)
        except Exception as ex:
            self._set_status(f"Save failed: {ex}", ok=False)

    def _on_export_hkl(self, _):
        fname = os.path.join(
            os.path.dirname(self.cif_path),
            f'reflections_{self._stub()}_zone{self._zonetag()}'
            f'_L{self.layer:+d}.txt'.replace('+', 'p').replace('-', 'm'),
        )
        try:
            with open(fname, 'w') as f:
                f.write(f"# CIF: {os.path.basename(self.cif_path)}\n")
                f.write(f"# zone axis : {self.zone_axis}\n")
                f.write(f"# layer L   : {self.layer:+d}\n")
                f.write(f"# d_min (A) : {self.d_min:.3f}\n")
                f.write(f"# h_max     : {self.h_max:d}\n")
                f.write(f"# twin (3)  : {self.twin_3}\n")
                f.write(f"# columns   : h k l d(A) |F| I I_norm\n")
                for r in sorted(self.reflections, key=lambda x: -x['I']):
                    f.write(f"{r['h']:4d} {r['k']:4d} {r['l']:4d} "
                            f"{r['d']:8.4f} {r['F']:10.3f} "
                            f"{r['I']:12.3e} {r['I_norm']:8.4f}\n")
            self._set_status(f"Exported hkl: {os.path.basename(fname)}",
                             ok=True)
        except Exception as ex:
            self._set_status(f"hkl export failed: {ex}", ok=False)

    def _on_export_matrix(self, _):
        ext = self.export_extent_rlu
        extent = (-ext, ext, -ext, ext)
        try:
            img = self.calc.rasterize_plane(
                self.reflections, extent=extent,
                n_pix=self.export_n_pix,
                sigma=self.export_sigma_rlu, use_norm=True,
            )
        except Exception as ex:
            self._set_status(f"raster failed: {ex}", ok=False)
            return

        base = os.path.join(
            os.path.dirname(self.cif_path),
            f'matrix_{self._stub()}_zone{self._zonetag()}'
            f'_L{self.layer:+d}'.replace('+', 'p').replace('-', 'm'),
        )
        npy_path = base + '.npy'
        csv_path = base + '.csv'
        meta_path = base + '_meta.txt'
        try:
            np.save(npy_path, img)
            np.savetxt(csv_path, img, delimiter=',', fmt='%.6e')
            with open(meta_path, 'w') as f:
                f.write(f"# CIF        : {os.path.basename(self.cif_path)}\n")
                f.write(f"# zone axis  : {self.zone_axis}\n")
                f.write(f"# layer L    : {self.layer:+d}\n")
                f.write(f"# extent rlu : {extent}\n")
                f.write(f"# shape      : {img.shape}\n")
                f.write(f"# sigma rlu  : {self.export_sigma_rlu}\n")
                f.write(f"# I_max      : {img.max():.3e}\n")
                f.write(f"# axis x     : k or k_perp1 (rlu, see zone)\n")
                f.write(f"# axis y     : k_perp2 (rlu, see zone)\n")
            self._set_status(
                f"Saved matrix: {os.path.basename(npy_path)} "
                f"({img.shape[0]}x{img.shape[1]})", ok=True,
            )
        except Exception as ex:
            self._set_status(f"matrix save failed: {ex}", ok=False)

    def _on_load_cif(self, _):
        path = _pick_cif_dialog()
        if not path:
            self._set_status('Load cancelled', ok=False)
            return
        try:
            self.cif_path = os.path.abspath(path)
            self.cif  = CIFParser(self.cif_path)
            self.calc = DiffractionCalculator(self.cif)
            self.composition_txt.set_text(self._composition_banner())
            self._redraw(full=True)
            self._set_status(f"Loaded: {os.path.basename(path)}", ok=True)
        except Exception as ex:
            self._set_status(f"Load failed: {ex}", ok=False)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    if len(sys.argv) > 1:
        cif = sys.argv[1]
    else:
        cif = _pick_cif_dialog()
        if not cif:
            print('No CIF given on command line and no file selected '
                  'in the dialog. Exiting.', file=sys.stderr)
            sys.exit(1)
    if not os.path.isfile(cif):
        print(f'CIF file not found: {cif}', file=sys.stderr)
        sys.exit(2)
    CIFViewer(cif).run()


if __name__ == '__main__':
    main()
