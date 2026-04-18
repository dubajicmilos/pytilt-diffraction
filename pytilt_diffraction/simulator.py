"""
INTERACTIVE GLAZER-TILT SINGLE-CRYSTAL DIFFRACTION
==================================================

Couples the pytilting octahedral-tilt generator (Glazer notation) to the
X-ray structure-factor / single-crystal diffraction engine.

Move the tilt-angle sliders or switch Glazer systems and the reciprocal
lattice pattern updates live. UI is modelled on CrystalMaker SingleCrystal.

Usage
-----
    python interactive_tilt_diffraction.py                 # defaults (CsPbBr3-like)
    python interactive_tilt_diffraction.py Cs Pb Br 5.874  # custom ABX3 + a0

Controls
--------
    * Glazer preset radio buttons: pick one of the 23 tilt systems
    * omega_x, omega_y, omega_z sliders: tilt magnitudes (degrees)
      -> magnitudes that the Glazer letters demand to be equal are
         automatically kept equal (driven by omega_x)
    * Zone-axis radio: view direction of reciprocal-space plane
    * d_min, spot-size, h_max, label-threshold sliders
    * Save PNG / Export reflection list buttons
"""

import sys
import os
import tempfile
import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons, TextBox

# Make the vendored pytilting package importable. Its modules use top-level
# imports (e.g. `from puc import Puc`) so the src/ folder itself has to be
# on sys.path rather than packaging it as a subpackage.
HERE = os.path.dirname(os.path.abspath(__file__))
_VENDOR_SRC = os.path.normpath(
    os.path.join(HERE, '..', 'vendor', 'pytilting', 'src')
)
if _VENDOR_SRC not in sys.path:
    sys.path.insert(0, _VENDOR_SRC)

from distortion import Distortion                              # noqa: E402
from pytilt_diffraction.calculator import (                    # noqa: E402
    CIFParser, DiffractionCalculator,
)


# ---------------------------------------------------------------------------
# The 23 Glazer tilt systems, grouped as in Howard & Stokes (1998)
# ---------------------------------------------------------------------------
GLAZER_PRESETS = [
    ('a0a0a0', 'Pm-3m',     'cubic (no tilt)'),
    ('a0a0c+', 'P4/mbm',    'single in-phase'),
    ('a0a0c-', 'I4/mcm',    'single out-of-phase'),
    ('a0b+b+', 'I4/mmm',    'two equal in-phase'),
    ('a0b-b-', 'Imcm',      'two equal out-of-phase'),
    ('a0b+c-', 'Cmcm',      'mixed in/out'),
    ('a0b-c-', 'C2/m',      'two unequal out-of-phase'),
    ('a+a+a+', 'Im-3',      'three equal in-phase'),
    ('a+b+b+', 'Immm',      'three in-phase'),
    ('a+b+c+', 'Immm',      'three unequal in-phase'),
    ('a+a+c-', 'P4_2/nmc',  'pair + one out'),
    ('a+b+c-', 'Pmmn',      'two in + one out'),
    ('a+b-b-', 'Pnma',      '*GdFeO3 / Pbnm*'),
    ('a+a-a-', 'Pnma',      'common perovskite'),
    ('a+b-c-', 'P2_1/m',    'one in + two out'),
    ('a-a-a-', 'R-3c',      'rhombohedral (LaAlO3)'),
    ('a-b-b-', 'I2/a',      'two equal out-of-phase'),
    ('a-b-c-', 'F-1',       'triclinic'),
]

# ---------------------------------------------------------------------------
# Material presets. High-T cubic aristotype lattice constants (pseudocubic
# for compounds that distort at room temperature). a0 can be overridden at
# runtime with the a0 slider.
# ---------------------------------------------------------------------------
MATERIALS = [
    # halide perovskites
    ('CsPbCl3', 'Cs', 'Pb', 'Cl', 5.605),
    ('CsPbBr3', 'Cs', 'Pb', 'Br', 5.874),
    ('CsPbI3',  'Cs', 'Pb', 'I',  6.289),
    ('CsSnBr3', 'Cs', 'Sn', 'Br', 5.800),
    ('CsGeBr3', 'Cs', 'Ge', 'Br', 5.636),
    ('RbPbBr3', 'Rb', 'Pb', 'Br', 5.850),
    # oxide perovskites
    ('SrTiO3',  'Sr', 'Ti', 'O',  3.905),
    ('BaTiO3',  'Ba', 'Ti', 'O',  4.004),
    ('CaTiO3',  'Ca', 'Ti', 'O',  3.838),
    ('PbTiO3',  'Pb', 'Ti', 'O',  3.970),
    ('BaZrO3',  'Ba', 'Zr', 'O',  4.193),
    ('LaAlO3',  'La', 'Al', 'O',  3.791),
    ('LaMnO3',  'La', 'Mn', 'O',  3.935),
    ('LaFeO3',  'La', 'Fe', 'O',  3.926),
    # fluoride perovskites
    ('KMgF3',   'K',  'Mg', 'F',  3.989),
    ('KNiF3',   'K',  'Ni', 'F',  4.012),
]


# Map Glazer sign/letter pattern -> which omega components must be equal
# Returns list of groups; each group's magnitudes are forced equal
def glazer_equality_groups(g):
    """
    Given a glazer string like 'a+b-b-', return [[0], [1,2]]  meaning
    omega[1] and omega[2] must share the same magnitude.
    """
    letters = [g[0], g[2], g[4]]
    groups = []
    seen = {}
    for i, L in enumerate(letters):
        seen.setdefault(L, []).append(i)
    for L in 'abc':
        if L in seen:
            groups.append(seen[L])
    return groups


def glazer_nonzero_axes(g):
    """Return list of axis indices (0/1/2) that carry a tilt (sign != '0')."""
    signs = [g[1], g[3], g[5]]
    return [i for i, s in enumerate(signs) if s != '0']


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------
class TiltDiffractionSimulator:
    """Interactive Glazer-tilt + single-crystal diffraction."""

    # Palette (SingleCrystal-ish: dark with warm spot colormap)
    BG         = '#0d1117'
    PANEL      = '#161b22'
    PLOT_BG    = '#000000'
    ACCENT     = '#ff6b35'
    ACCENT_2   = '#4a90e2'
    TEXT       = '#e6edf3'
    TEXT_DIM   = '#8b949e'
    GRID       = '#21262d'
    SPOT_CMAP  = 'hot'

    def __init__(self, symbols=('Cs', 'Pb', 'Br'), a0=5.874, grid=(2, 2, 2)):
        self.symbols = list(symbols)
        self.a0 = float(a0)
        self.grid = tuple(grid)

        self.distortion = Distortion(system={
            'symbols': self.symbols,
            'lattice_constant': self.a0,
            'grid': self.grid,
            'covera': 1.0,
        })
        self.tmp_cif = os.path.join(tempfile.gettempdir(), 'pytilt_live.cif')

        # Diffraction / view parameters
        self.glazer     = 'a0a0a0'
        self.omega_deg  = [6.0, 6.0, 6.0]    # magnitudes in degrees
        self.zone_axis  = (0, 0, 1)
        self.h_max      = 8
        self.d_min      = 0.8
        self.I_min      = 1e-3
        self.spot_scale = 250.0
        self.show_labels   = True
        self.label_thresh  = 0.05
        self.log_scale     = False
        # HKx layer selector: integer layer index in *supercell* Miller units.
        # For the default (2,2,2) grid this means half-integer steps in the
        # parent pseudocubic cell (layer=1 supercell ~ 0.5 parent).
        self.layer      = 0

        self.cif  = None
        self.calc = None
        self.reflections = []

        # Prevent slider callbacks from firing during programmatic updates
        self._suppress = False

    # ------------------------------------------------------------------
    # Core pipeline: tilt -> CIF -> diffraction
    # ------------------------------------------------------------------
    def _effective_omega_rad(self):
        """Apply Glazer consistency: equal letters => equal magnitudes,
        different letters => distinct magnitudes (small nudge if they
        happen to collide), '0' signs => zero.  Returns radians."""
        groups = glazer_equality_groups(self.glazer)
        signs  = [self.glazer[1], self.glazer[3], self.glazer[5]]
        mags = [abs(a) for a in self.omega_deg]

        # 1) equalise within groups (drive by first member)
        for grp in groups:
            master = mags[grp[0]]
            for i in grp:
                mags[i] = master

        # 2) enforce distinct magnitudes between different-letter groups
        #    nudge the second one down by 0.5 deg if they collide
        for gi in range(len(groups)):
            for gj in range(gi + 1, len(groups)):
                vi, vj = mags[groups[gi][0]], mags[groups[gj][0]]
                if abs(vi - vj) < 1e-3:
                    vj_new = max(0.0, vj - 0.5) if vj > 0.5 else vj + 0.5
                    for i in groups[gj]:
                        mags[i] = vj_new

        # 3) zero-out where sign == '0'
        for i in range(3):
            if signs[i] == '0':
                mags[i] = 0.0

        return [math.radians(m) for m in mags]

    def set_composition(self, A, B, X, a0=None):
        """Swap out the ABX3 basis and lattice constant while keeping the
        current Glazer tilt pattern, zone axis, and layer settings.

        Rebuilds the pytilting Distortion object, regenerates the tilted
        structure, and refreshes the diffraction cache."""
        self.symbols = [A, B, X]
        if a0 is not None:
            self.a0 = float(a0)
        self.distortion = Distortion(system={
            'symbols': self.symbols,
            'lattice_constant': self.a0,
            'grid': self.grid,
            'covera': 1.0,
        })
        self.rebuild_structure()
        self.recompute_pattern()

    def rebuild_structure(self):
        omega = self._effective_omega_rad()
        self.distortion.distort = {
            'glazer':     self.glazer,
            'omega':      tuple(omega),
            'u':          (0.0, 0.0, 0.0),
            'k_u':        2 * math.pi * np.zeros((3, 3)),
            'local_mode': [0.0] * 5,
            'modes':      [],
        }
        atoms = self.distortion.get_atoms()
        atoms.write(self.tmp_cif)

        # Reload via the existing CIF / diffraction machinery
        self.cif  = CIFParser(self.tmp_cif)
        self.calc = DiffractionCalculator(self.cif)

        # Global intensity reference: strongest reflection at the fundamental
        # layer (L=0). Used to keep I_norm absolute across layers so that
        # layers containing only numerical-noise peaks get filtered out.
        self._update_I_max_ref()

    def _update_I_max_ref(self):
        ref = self.calc.get_plane_reflections(
            self.zone_axis, self.h_max, self.d_min, I_min=0.0, layer=0,
        )
        self.I_max_ref = max((r['I'] for r in ref), default=1.0)

    def recompute_pattern(self):
        self.reflections = self.calc.get_plane_reflections(
            self.zone_axis, self.h_max, self.d_min, self.I_min,
            layer=self.layer, I_max_ref=self.I_max_ref,
        )
        self.reflections = self.calc.get_2d_coordinates(
            self.reflections, self.zone_axis,
        )

    def _parent_layer(self):
        """Convert supercell integer layer to parent pseudocubic fractional
        layer. For default (2,2,2) grid, L_parent = L_super / 2."""
        u, v, w = self.zone_axis
        # use the grid component(s) relevant to the zone axis
        comps = []
        if u != 0: comps.append(self.grid[0])
        if v != 0: comps.append(self.grid[1])
        if w != 0: comps.append(self.grid[2])
        n = min(comps) if comps else 1
        return self.layer / n if n else 0.0

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------
    def run(self):
        self.rebuild_structure()
        self.recompute_pattern()

        plt.style.use('dark_background')
        self.fig = plt.figure(figsize=(17, 10))
        self.fig.canvas.manager.set_window_title(
            'pytilting  .  Interactive Glazer-tilt diffraction'
        )
        self.fig.patch.set_facecolor(self.BG)

        # --- main reciprocal-space plot ---------------------------------
        self.ax = self.fig.add_axes([0.05, 0.30, 0.52, 0.63])
        self.ax.set_facecolor(self.PLOT_BG)

        # --- info / readout panels --------------------------------------
        self.info_ax = self.fig.add_axes([0.60, 0.56, 0.17, 0.37])
        self.info_ax.set_facecolor(self.PANEL)
        self.info_ax.set_xticks([]); self.info_ax.set_yticks([])
        for s in self.info_ax.spines.values():
            s.set_color(self.GRID)

        self.refl_ax = self.fig.add_axes([0.79, 0.56, 0.19, 0.37])
        self.refl_ax.set_facecolor(self.PANEL)
        self.refl_ax.set_xticks([]); self.refl_ax.set_yticks([])
        for s in self.refl_ax.spines.values():
            s.set_color(self.GRID)

        # --- title banner ------------------------------------------------
        self.fig.text(0.05, 0.955,
                      'Glazer-tilt Single-Crystal Diffraction Simulator',
                      color=self.TEXT, fontsize=15, fontweight='bold')
        self.composition_txt = self.fig.text(
            0.05, 0.930, self._composition_banner(),
            color=self.TEXT_DIM, fontsize=10,
        )
        # status banner for save/export feedback
        self.status_txt = self.fig.text(
            0.60, 0.955, '', color=self.ACCENT,
            fontsize=9, va='center',
        )

        # Reference to optional side-by-side powder diffraction window; see
        # _on_powder. The main update path keeps this in sync on tilt / a0 /
        # composition changes.
        self.powder_win = None

        self._build_controls()
        self._draw_pattern()
        self._draw_info()
        self._draw_reflections()

        plt.show()

    def _build_controls(self):
        # =======   Glazer preset radio  =======
        ax_lbl = self.fig.add_axes([0.60, 0.49, 0.38, 0.03]); ax_lbl.axis('off')
        ax_lbl.text(0, 0.5, 'Glazer tilt system  (sign pattern  ->  space group)',
                    color=self.TEXT, fontsize=10, fontweight='bold', va='center')

        # Split into two columns of radio buttons to fit 18 entries
        n = len(GLAZER_PRESETS)
        half = (n + 1) // 2
        col1 = GLAZER_PRESETS[:half]
        col2 = GLAZER_PRESETS[half:]
        col1_labels = [f"{g}   {sg}" for g, sg, _ in col1]
        col2_labels = [f"{g}   {sg}" for g, sg, _ in col2]

        self.ax_glazer1 = self.fig.add_axes([0.60, 0.05, 0.09, 0.43])
        self.ax_glazer1.set_facecolor(self.PANEL)
        self.radio_glazer1 = RadioButtons(
            self.ax_glazer1, col1_labels, active=0, activecolor=self.ACCENT,
        )
        self._style_radio(self.radio_glazer1)
        self.radio_glazer1.on_clicked(lambda lab: self._on_glazer(lab, column=1))

        self.ax_glazer2 = self.fig.add_axes([0.70, 0.05, 0.09, 0.43])
        self.ax_glazer2.set_facecolor(self.PANEL)
        self.radio_glazer2 = RadioButtons(
            self.ax_glazer2, col2_labels, active=0, activecolor=self.ACCENT,
        )
        self._style_radio(self.radio_glazer2)
        self.radio_glazer2.on_clicked(lambda lab: self._on_glazer(lab, column=2))
        # At startup col1 is the active column ('a0a0a0'), so clear col2's dot
        self._visual_deselect(self.radio_glazer2)

        # =======   Zone axis radio   =======
        zone_lbl = self.fig.add_axes([0.80, 0.53, 0.18, 0.02]); zone_lbl.axis('off')
        zone_lbl.text(0, 0.5, 'Zone axis  (view direction)',
                      color=self.TEXT, fontsize=10, fontweight='bold', va='center')

        self.ax_zone = self.fig.add_axes([0.80, 0.41, 0.18, 0.12])
        self.ax_zone.set_facecolor(self.PANEL)
        self.zone_labels = [
            '[001]  hk0',
            '[100]  0kl',
            '[010]  h0l',
            '[110]  hh-l',
            '[111]  (h+k+l)=0',
        ]
        self.zone_map = {
            self.zone_labels[0]: (0, 0, 1),
            self.zone_labels[1]: (1, 0, 0),
            self.zone_labels[2]: (0, 1, 0),
            self.zone_labels[3]: (1, 1, 0),
            self.zone_labels[4]: (1, 1, 1),
        }
        self.radio_zone = RadioButtons(
            self.ax_zone, self.zone_labels, active=0, activecolor=self.ACCENT_2,
        )
        self._style_radio(self.radio_zone)
        self.radio_zone.on_clicked(self._on_zone)

        # =======   Material preset radios  =======
        mat_lbl = self.fig.add_axes([0.80, 0.39, 0.18, 0.02]); mat_lbl.axis('off')
        mat_lbl.text(0, 0.5, 'Material (ABX3 preset)',
                     color=self.TEXT, fontsize=10, fontweight='bold', va='center')

        n_mat = len(MATERIALS)
        half = (n_mat + 1) // 2
        col1 = MATERIALS[:half]
        col2 = MATERIALS[half:]
        col1_labels = [f"{m[0]}   a={m[4]:.3f}" for m in col1]
        col2_labels = [f"{m[0]}   a={m[4]:.3f}" for m in col2]

        self.ax_mat1 = self.fig.add_axes([0.80, 0.23, 0.09, 0.16])
        self.ax_mat1.set_facecolor(self.PANEL)
        self.radio_mat1 = RadioButtons(
            self.ax_mat1, col1_labels, active=0, activecolor=self.ACCENT,
        )
        self._style_radio(self.radio_mat1)
        self.radio_mat1.on_clicked(lambda lab: self._on_material(lab, column=1))

        self.ax_mat2 = self.fig.add_axes([0.89, 0.23, 0.09, 0.16])
        self.ax_mat2.set_facecolor(self.PANEL)
        self.radio_mat2 = RadioButtons(
            self.ax_mat2, col2_labels, active=0, activecolor=self.ACCENT,
        )
        self._style_radio(self.radio_mat2)
        self.radio_mat2.on_clicked(lambda lab: self._on_material(lab, column=2))
        # Figure out which column currently holds the active preset; keep the
        # other deselected so only one radio dot shows.
        self._sync_material_active_radio()

        # =======   Lattice a0 slider  =======
        ax_a0 = self.fig.add_axes([0.80, 0.195, 0.18, 0.020])
        self.s_a0 = Slider(ax_a0, 'a0 (A)', 3.0, 7.5,
                           valinit=self.a0, valstep=0.001,
                           color=self.ACCENT)
        self.s_a0.label.set_color(self.TEXT)
        self.s_a0.valtext.set_color(self.TEXT)
        self.s_a0.on_changed(self._on_a0)

        # =======   Tilt-angle sliders (degrees)   =======
        sl_color = self.ACCENT
        ax_wx = self.fig.add_axes([0.08, 0.23, 0.45, 0.020])
        self.s_wx = Slider(ax_wx, 'omega_x  (deg)', 0.0, 20.0,
                           valinit=self.omega_deg[0], color=sl_color)
        ax_wy = self.fig.add_axes([0.08, 0.20, 0.45, 0.020])
        self.s_wy = Slider(ax_wy, 'omega_y  (deg)', 0.0, 20.0,
                           valinit=self.omega_deg[1], color=sl_color)
        ax_wz = self.fig.add_axes([0.08, 0.17, 0.45, 0.020])
        self.s_wz = Slider(ax_wz, 'omega_z  (deg)', 0.0, 20.0,
                           valinit=self.omega_deg[2], color=sl_color)
        for s in (self.s_wx, self.s_wy, self.s_wz):
            s.label.set_color(self.TEXT)
            s.valtext.set_color(self.TEXT)
        self.s_wx.on_changed(lambda v: self._on_omega(0, v))
        self.s_wy.on_changed(lambda v: self._on_omega(1, v))
        self.s_wz.on_changed(lambda v: self._on_omega(2, v))

        # =======   HKx layer slider (integer, supercell units)   =======
        ax_layer = self.fig.add_axes([0.08, 0.14, 0.45, 0.020])
        self.s_layer = Slider(ax_layer, 'HK_L layer (super)', -4, 4,
                              valinit=self.layer, valstep=1,
                              color=self.ACCENT)
        self.s_layer.label.set_color(self.TEXT)
        self.s_layer.valtext.set_color(self.TEXT)
        self.s_layer.on_changed(self._on_layer)

        # Parent-cell annotation (e.g. "parent L = 0.50")
        self.layer_ax = self.fig.add_axes([0.54, 0.14, 0.05, 0.020])
        self.layer_ax.axis('off')
        self.layer_txt = self.layer_ax.text(
            0, 0.5, self._layer_hint_text(),
            color=self.ACCENT_2, fontsize=8, va='center',
        )

        # =======   Viewing sliders   =======
        ax_dmin = self.fig.add_axes([0.08, 0.11, 0.45, 0.020])
        self.s_dmin = Slider(ax_dmin, 'd_min (A)', 0.3, 2.5,
                             valinit=self.d_min, color=self.ACCENT_2)
        ax_hmax = self.fig.add_axes([0.08, 0.08, 0.45, 0.020])
        self.s_hmax = Slider(ax_hmax, 'h_max', 3, 14,
                             valinit=self.h_max, valstep=1, color=self.ACCENT_2)
        ax_spot = self.fig.add_axes([0.08, 0.05, 0.45, 0.020])
        self.s_spot = Slider(ax_spot, 'spot size', 30.0, 800.0,
                             valinit=self.spot_scale, color=self.ACCENT_2)
        ax_lthr = self.fig.add_axes([0.08, 0.02, 0.45, 0.020])
        self.s_lthr = Slider(ax_lthr, 'label I/I_max >=', 0.0, 0.5,
                             valinit=self.label_thresh, color=self.ACCENT_2)
        for s in (self.s_dmin, self.s_hmax, self.s_spot, self.s_lthr):
            s.label.set_color(self.TEXT)
            s.valtext.set_color(self.TEXT)

        self.s_dmin.on_changed(self._on_dmin)
        self.s_hmax.on_changed(self._on_hmax)
        self.s_spot.on_changed(self._on_spot)
        self.s_lthr.on_changed(self._on_lthr)

        # =======   Check boxes + buttons   =======
        ax_chk = self.fig.add_axes([0.80, 0.140, 0.18, 0.045])
        ax_chk.set_facecolor(self.PANEL)
        self.check = CheckButtons(ax_chk, ['Show hkl labels', 'log intensity'],
                                  [self.show_labels, self.log_scale])
        for lbl in self.check.labels:
            lbl.set_color(self.TEXT); lbl.set_fontsize(8)
        self.check.on_clicked(self._on_check)

        ax_save   = self.fig.add_axes([0.80, 0.098, 0.085, 0.032])
        self.btn_save = Button(ax_save, 'Save PNG',
                               color=self.ACCENT, hovercolor='#ff8c5c')
        self.btn_save.label.set_color('white')
        self.btn_save.on_clicked(self._on_save)

        ax_export = self.fig.add_axes([0.895, 0.098, 0.085, 0.032])
        self.btn_export = Button(ax_export, 'Export hkl',
                                 color=self.ACCENT_2, hovercolor='#6ba7e8')
        self.btn_export.label.set_color('white')
        self.btn_export.on_clicked(self._on_export)

        ax_powder = self.fig.add_axes([0.80, 0.058, 0.085, 0.032])
        self.btn_powder = Button(ax_powder, 'Powder',
                                 color='#8b5cf6', hovercolor='#a78bfa')
        self.btn_powder.label.set_color('white')
        self.btn_powder.on_clicked(self._on_powder)

        ax_reset  = self.fig.add_axes([0.895, 0.058, 0.085, 0.032])
        self.btn_reset = Button(ax_reset, 'Reset',
                                color=self.GRID, hovercolor='#30363d')
        self.btn_reset.label.set_color(self.TEXT)
        self.btn_reset.on_clicked(self._on_reset)

    def _style_radio(self, radio):
        for c in radio.circles if hasattr(radio, 'circles') else []:
            c.set_edgecolor(self.TEXT_DIM)
        for lbl in radio.labels:
            lbl.set_color(self.TEXT)
            lbl.set_fontsize(8)
        for s in radio.ax.spines.values():
            s.set_color(self.GRID)

    def _visual_deselect(self, radio):
        """Clear the visual "active" dot on a RadioButtons without firing its
        callback. matplotlib has no built-in deselect, but we can paint every
        button transparent via the PathCollection it uses internally.

        Use float tuples — integer (0,0,0,0) tuples cause matplotlib to store
        the facecolors array as int dtype, which then truncates the next
        activecolor write (e.g. '#ff6b35' -> (1,0,0,1) instead of orange)."""
        if hasattr(radio, '_buttons') and hasattr(radio._buttons, 'set_facecolors'):
            n = len(radio.labels)
            radio._buttons.set_facecolors([(0.0, 0.0, 0.0, 0.0)] * n)
        elif hasattr(radio, 'circles'):
            for c in radio.circles:
                c.set_facecolor('none')

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------
    def _on_glazer(self, label, column):
        g = label.split()[0]
        self.glazer = g

        # Only one radio-dot should be visible at a time, but matplotlib's
        # RadioButtons cannot natively clear its selection. Paint the other
        # column's dots transparent so the active column is unambiguous.
        other = self.radio_glazer2 if column == 1 else self.radio_glazer1
        self._visual_deselect(other)

        # Disable '0' sliders visually, realign tied magnitudes
        self._sync_sliders_to_glazer()

        self._recompute_and_redraw()

    def _sync_sliders_to_glazer(self):
        """Enforce Glazer consistency: equal letters -> equal magnitudes,
        different letters -> distinct magnitudes, '0' -> 0.  Distinct-
        group magnitudes are seeded from canonical ratios so that switching
        systems never produces a decode_glazer consistency failure."""
        self._suppress = True
        groups = glazer_equality_groups(self.glazer)
        signs  = [self.glazer[1], self.glazer[3], self.glazer[5]]

        # Scale base from current omega_x if sensible, else default
        base = self.omega_deg[0] if self.omega_deg[0] > 0.5 else 6.0
        canonical = [base, base * 0.6, base * 0.35]

        mags = [0.0, 0.0, 0.0]
        for gi, grp in enumerate(groups):
            v = canonical[gi] if gi < len(canonical) else canonical[-1]
            for i in grp:
                mags[i] = v
        for i in range(3):
            if signs[i] == '0':
                mags[i] = 0.0

        self.omega_deg = mags
        self.s_wx.set_val(mags[0])
        self.s_wy.set_val(mags[1])
        self.s_wz.set_val(mags[2])
        self._suppress = False

    def _on_omega(self, axis, val):
        if self._suppress:
            return
        self._suppress = True

        # If the user drags an axis forbidden by '0', snap it back to 0
        signs = [self.glazer[1], self.glazer[3], self.glazer[5]]
        if signs[axis] == '0':
            self.omega_deg[axis] = 0.0
            [self.s_wx, self.s_wy, self.s_wz][axis].set_val(0.0)
            self._suppress = False
            return

        self.omega_deg[axis] = val

        # Propagate to tied axes (same Glazer letter)
        groups = glazer_equality_groups(self.glazer)
        for grp in groups:
            if axis in grp:
                for i in grp:
                    if i != axis:
                        self.omega_deg[i] = val
                        [self.s_wx, self.s_wy, self.s_wz][i].set_val(val)

        self._suppress = False
        self._recompute_and_redraw()

    def _composition_banner(self):
        A, B, X = self.symbols
        return (
            f"Composition  {A}{B}{X}3   .   "
            f"a0 = {self.a0:.4f} A   .   "
            f"grid = {self.grid[0]}x{self.grid[1]}x{self.grid[2]}"
        )

    def _update_composition_banner(self):
        self.composition_txt.set_text(self._composition_banner())

    def _sync_material_active_radio(self):
        """Match the radio highlight to the currently-loaded composition and
        clear the opposite column's dot so only one radio shows as active.

        In current matplotlib versions `set_active` fires the observer, which
        would re-trigger `_on_material` and cascade through sliders that may
        not exist yet during `_build_controls`. We suppress that while we
        drive the widget programmatically."""
        tag = ''.join(self.symbols) + '3'
        n_half = (len(MATERIALS) + 1) // 2
        idx = None
        for i, m in enumerate(MATERIALS):
            if m[0] == tag:
                idx = i
                break
        if idx is None:
            # Composition is not in the preset list (custom CLI start).
            # Leave both columns visually empty so nothing is highlighted.
            self._visual_deselect(self.radio_mat1)
            self._visual_deselect(self.radio_mat2)
            return
        self._suppress = True
        try:
            if idx < n_half:
                self.radio_mat1.set_active(idx)
                self._visual_deselect(self.radio_mat2)
            else:
                self.radio_mat2.set_active(idx - n_half)
                self._visual_deselect(self.radio_mat1)
        finally:
            self._suppress = False

    def _on_material(self, label, column):
        """User clicked a preset. Swap composition, snap a0 slider, keep
        everything else (tilt / zone / layer) as-is."""
        if self._suppress:
            return
        name = label.split()[0]
        entry = next((m for m in MATERIALS if m[0] == name), None)
        if entry is None:
            return
        _, A, B, X, a0 = entry

        other = self.radio_mat2 if column == 1 else self.radio_mat1
        self._visual_deselect(other)

        self.set_composition(A, B, X, a0)

        self._suppress = True
        self.s_a0.set_val(a0)
        self._suppress = False

        self._update_composition_banner()
        self._draw_pattern()
        self._draw_info()
        self._draw_reflections()
        self._refresh_powder_if_open()
        self.fig.canvas.draw_idle()

    def _on_a0(self, val):
        if self._suppress:
            return
        a0 = float(val)
        self.set_composition(self.symbols[0], self.symbols[1], self.symbols[2], a0)
        self._update_composition_banner()
        self._draw_pattern()
        self._draw_info()
        self._draw_reflections()
        self._refresh_powder_if_open()
        self.fig.canvas.draw_idle()

    def _on_powder(self, _):
        """Open (or raise) the powder-diffraction side window."""
        if self.powder_win is not None and plt.fignum_exists(self.powder_win.fig.number):
            # Already open: just bring it to front.
            try:
                self.powder_win.fig.canvas.manager.show()
            except Exception:
                pass
            return
        self.powder_win = PowderWindow(self)
        self._set_status('Powder window opened', ok=True)

    def _refresh_powder_if_open(self):
        if self.powder_win is None:
            return
        if not plt.fignum_exists(self.powder_win.fig.number):
            self.powder_win = None
            return
        self.powder_win.refresh()

    def _on_zone(self, label):
        self.zone_axis = self.zone_map[label]
        self._update_I_max_ref()
        self.recompute_pattern()
        self._update_layer_hint()
        self._draw_pattern()
        self._draw_info()
        self._draw_reflections()
        self.fig.canvas.draw_idle()

    def _on_layer(self, val):
        if self._suppress:
            return
        self.layer = int(round(val))
        self.recompute_pattern()
        self._update_layer_hint()
        self._draw_pattern()
        self._draw_info()
        self._draw_reflections()
        self.fig.canvas.draw_idle()

    def _layer_hint_text(self):
        return f"= {self._parent_layer():+.2f} parent"

    def _update_layer_hint(self):
        if hasattr(self, 'layer_txt'):
            self.layer_txt.set_text(self._layer_hint_text())

    def _on_dmin(self, val):
        if self._suppress:
            return
        self.d_min = float(val)
        self._update_I_max_ref()
        self.recompute_pattern()
        self._draw_pattern(); self._draw_info(); self._draw_reflections()
        self.fig.canvas.draw_idle()

    def _on_hmax(self, val):
        if self._suppress:
            return
        self.h_max = int(val)
        self._update_I_max_ref()
        self.recompute_pattern()
        self._draw_pattern(); self._draw_info(); self._draw_reflections()
        self.fig.canvas.draw_idle()

    def _on_spot(self, val):
        self.spot_scale = float(val)
        self._draw_pattern(); self.fig.canvas.draw_idle()

    def _on_lthr(self, val):
        self.label_thresh = float(val)
        self._draw_pattern(); self.fig.canvas.draw_idle()

    def _on_check(self, label):
        if 'labels' in label:
            self.show_labels = not self.show_labels
        elif 'log' in label:
            self.log_scale = not self.log_scale
        self._draw_pattern(); self.fig.canvas.draw_idle()

    def _set_status(self, msg, ok=True):
        """Flash a message in the top banner so the user sees save/export
        feedback even when the terminal isn't visible (e.g. IDE run)."""
        self.status_txt.set_text(msg)
        self.status_txt.set_color(self.ACCENT if ok else '#f85149')
        self.fig.canvas.draw_idle()

    def _on_save(self, _):
        tag = self.glazer.replace('+', 'p').replace('-', 'm')
        zone = ''.join(str(i) for i in self.zone_axis)
        lay = f"L{self.layer:+d}".replace('+', 'p').replace('-', 'm')
        fname = os.path.join(
            HERE, f'diffraction_{tag}_zone{zone}_{lay}.png'
        )
        try:
            self.fig.savefig(fname, dpi=200, facecolor=self.BG,
                             bbox_inches='tight')
            print(f"[save] wrote {fname}", flush=True)
            self._set_status(f"Saved PNG: {os.path.basename(fname)}", ok=True)
        except Exception as ex:
            print(f"[save] FAILED: {ex}", flush=True)
            self._set_status(f"Save failed: {ex}", ok=False)

    def _on_export(self, _):
        tag = self.glazer.replace('+', 'p').replace('-', 'm')
        lay = f"L{self.layer:+d}".replace('+', 'p').replace('-', 'm')
        fname = os.path.join(HERE, f'reflections_{tag}_{lay}.txt')
        try:
            with open(fname, 'w') as f:
                f.write(f"# Glazer: {self.glazer}\n")
                f.write(f"# omega (deg): {self.omega_deg}\n")
                f.write(f"# zone axis : {self.zone_axis}\n")
                f.write(f"# layer     : L_super = {self.layer:+d}   "
                        f"L_parent = {self._parent_layer():+.3f}\n")
                f.write(f"# d_min     : {self.d_min}\n")
                f.write(f"# {'h':>4} {'k':>4} {'l':>4} {'d(A)':>8} "
                        f"{'|F|':>10} {'I':>12} {'I/Imax':>8}\n")
                for r in sorted(self.reflections, key=lambda x: -x['I']):
                    f.write(
                        f"  {r['h']:4d} {r['k']:4d} {r['l']:4d} "
                        f"{r['d']:8.4f} {r['F']:10.2f} {r['I']:12.2f} "
                        f"{r['I_norm']:8.4f}\n"
                    )
            print(f"[export] wrote {fname}", flush=True)
            self._set_status(f"Exported hkl: {os.path.basename(fname)}", ok=True)
        except Exception as ex:
            print(f"[export] FAILED: {ex}", flush=True)
            self._set_status(f"Export failed: {ex}", ok=False)

    def _on_reset(self, _):
        self.glazer = 'a0a0a0'
        self.omega_deg = [0.0, 0.0, 0.0]
        self.layer = 0
        self._sync_sliders_to_glazer()
        self._suppress = True
        self.s_layer.set_val(0)
        self._suppress = False
        self._update_layer_hint()
        self.radio_glazer1.set_active(0)
        self._visual_deselect(self.radio_glazer2)
        self._recompute_and_redraw()

    # ------------------------------------------------------------------
    def _recompute_and_redraw(self):
        try:
            self.rebuild_structure()
            self.recompute_pattern()
        except SystemExit:
            # decode_glazer calls exit() on bad input; catch it
            print("[warn] invalid Glazer/omega combination; ignoring.")
            return
        self._draw_pattern()
        self._draw_info()
        self._draw_reflections()
        self._refresh_powder_if_open()
        self.fig.canvas.draw_idle()

    # ------------------------------------------------------------------
    # Drawing
    # ------------------------------------------------------------------
    def _draw_pattern(self):
        self.ax.clear()
        self.ax.set_facecolor(self.PLOT_BG)

        if not self.reflections:
            self.ax.text(0.5, 0.5, 'no reflections pass d_min / I_min',
                         color=self.TEXT_DIM, ha='center', va='center',
                         transform=self.ax.transAxes)
        else:
            x = np.array([r['x'] for r in self.reflections])
            y = np.array([r['y'] for r in self.reflections])
            I = np.array([r['I_norm'] for r in self.reflections])
            if self.log_scale:
                Iplot = np.log10(I * 999 + 1) / np.log10(1000)
            else:
                Iplot = I
            sizes = self.spot_scale * np.sqrt(Iplot) + 3.0

            self.ax.scatter(
                x, y, s=sizes, c=Iplot, cmap=self.SPOT_CMAP,
                alpha=0.95, edgecolors='white', linewidths=0.25,
                vmin=0, vmax=1,
            )

            if self.show_labels:
                for r in self.reflections:
                    if r['I_norm'] >= self.label_thresh:
                        self.ax.annotate(
                            f"{r['h']},{r['k']},{r['l']}",
                            (r['x'], r['y']),
                            xytext=(3, 3), textcoords='offset points',
                            fontsize=6, color='#c9d1d9', alpha=0.85,
                        )

        # central 000 spot
        self.ax.scatter([0], [0], s=110, c='#58a6ff', marker='o',
                        edgecolors='white', linewidths=2, zorder=10)
        self.ax.annotate('000', (0, 0), xytext=(6, 6),
                         textcoords='offset points',
                         color='#58a6ff', fontsize=8, fontweight='bold')

        za = ''.join(str(i) for i in self.zone_axis)
        self.ax.set_title(
            f"zone [{za}]   L_s={self.layer:+d}  (parent {self._parent_layer():+.2f})"
            f"   .   Glazer {self.glazer}   .   "
            f"omega = ({self.omega_deg[0]:.2f}, {self.omega_deg[1]:.2f}, "
            f"{self.omega_deg[2]:.2f}) deg",
            color=self.TEXT, fontsize=11, pad=10,
        )
        self.ax.set_aspect('equal')
        self.ax.grid(True, color=self.GRID, alpha=0.5, linestyle='-', linewidth=0.4)
        self.ax.axhline(0, color=self.GRID, linewidth=0.6)
        self.ax.axvline(0, color=self.GRID, linewidth=0.6)
        self.ax.tick_params(colors=self.TEXT_DIM, labelsize=8)
        for s in self.ax.spines.values():
            s.set_color(self.GRID)
        self.ax.set_xlabel('reciprocal h*', color=self.TEXT_DIM)
        self.ax.set_ylabel('reciprocal k*', color=self.TEXT_DIM)

    def _draw_info(self):
        self.info_ax.clear()
        self.info_ax.axis('off')
        self.info_ax.set_facecolor(self.PANEL)

        cp = self.cif.cell_params
        n_refl = len(self.reflections)
        info = (
            f"STRUCTURE\n"
            f"---------\n"
            f"formula    {self.symbols[0]}{self.symbols[1]}{self.symbols[2]}3\n"
            f"space grp  {self.cif.space_group}\n"
            f"atoms      {len(self.cif.atoms)}\n"
            f"\n"
            f"UNIT CELL (supercell)\n"
            f"---------------------\n"
            f"a = {cp['a']:8.4f} A\n"
            f"b = {cp['b']:8.4f} A\n"
            f"c = {cp['c']:8.4f} A\n"
            f"alpha = {cp['alpha']:6.2f} deg\n"
            f"beta  = {cp['beta']:6.2f} deg\n"
            f"gamma = {cp['gamma']:6.2f} deg\n"
            f"\n"
            f"TILT\n"
            f"----\n"
            f"glazer   {self.glazer}\n"
            f"omega_x  {self.omega_deg[0]:6.2f} deg\n"
            f"omega_y  {self.omega_deg[1]:6.2f} deg\n"
            f"omega_z  {self.omega_deg[2]:6.2f} deg\n"
            f"\n"
            f"PATTERN\n"
            f"-------\n"
            f"reflections  {n_refl}\n"
            f"d_min (A)   {self.d_min:5.2f}\n"
            f"layer  L_s = {self.layer:+d}\n"
            f"       L_p = {self._parent_layer():+5.2f}\n"
        )
        self.info_ax.text(
            0.05, 0.97, info, transform=self.info_ax.transAxes,
            fontsize=8.5, color=self.TEXT, family='monospace',
            va='top',
        )

    def _draw_reflections(self):
        self.refl_ax.clear()
        self.refl_ax.axis('off')
        self.refl_ax.set_facecolor(self.PANEL)

        title = "STRONGEST REFLECTIONS\n" \
                "---------------------\n" \
                "  h  k  l   d(A)   I/Imax\n"
        lines = [title]
        top = sorted(self.reflections, key=lambda r: -r['I_norm'])[:18]
        for r in top:
            lines.append(
                f" {r['h']:>2} {r['k']:>2} {r['l']:>2}   "
                f"{r['d']:5.3f}   {r['I_norm']:6.3f}"
            )
        self.refl_ax.text(
            0.05, 0.97, '\n'.join(lines),
            transform=self.refl_ax.transAxes,
            fontsize=8.5, color=self.TEXT, family='monospace',
            va='top',
        )


# ---------------------------------------------------------------------------
# Powder diffraction side window
# ---------------------------------------------------------------------------
class PowderWindow:
    """Side window that renders a kinematic powder pattern for the parent
    simulator's current structure.

    Reads the parent's calc / symbols / a0 / glazer state each refresh, so
    the powder pattern stays in lockstep with composition, lattice-constant,
    and tilt changes in the main window.
    """

    # Common laboratory wavelengths (A)
    WAVELENGTHS = [
        ('Cu Ka1', 1.5406),
        ('Mo Ka1', 0.7093),
        ('Co Ka1', 1.7890),
        ('Cr Ka1', 2.2897),
        ('Ag Ka1', 0.5594),
    ]

    def __init__(self, parent):
        self.parent = parent
        self.wavelength = 1.5406
        self.two_theta_max = 90.0
        self.h_max = 10
        self.fwhm = 0.10
        self.n_labels = 8
        self.peaks = []

        p = parent
        plt.style.use('dark_background')
        self.fig = plt.figure(figsize=(11, 6.5))
        self.fig.canvas.manager.set_window_title(
            'pytilting  .  Powder diffraction'
        )
        self.fig.patch.set_facecolor(p.BG)

        self.ax = self.fig.add_axes([0.07, 0.28, 0.68, 0.63])
        self.ax.set_facecolor(p.PLOT_BG)
        for s in self.ax.spines.values():
            s.set_color(p.GRID)
        self.ax.tick_params(colors=p.TEXT_DIM)

        self.title_txt = self.fig.text(
            0.07, 0.945, '', color=p.TEXT, fontsize=12, fontweight='bold',
        )
        self.sub_txt = self.fig.text(
            0.07, 0.918, '', color=p.TEXT_DIM, fontsize=9,
        )

        # peak-list panel
        self.list_ax = self.fig.add_axes([0.77, 0.28, 0.21, 0.63])
        self.list_ax.set_facecolor(p.PANEL)
        self.list_ax.set_xticks([]); self.list_ax.set_yticks([])
        for s in self.list_ax.spines.values():
            s.set_color(p.GRID)

        self._build_controls()
        self.refresh()

        self.fig.canvas.mpl_connect('close_event', self._on_close)
        plt.show(block=False)
        try:
            self.fig.canvas.manager.show()
        except Exception:
            pass

    def _build_controls(self):
        p = self.parent

        # Wavelength radio (X-ray source)
        src_lbl = self.fig.add_axes([0.07, 0.21, 0.25, 0.025]); src_lbl.axis('off')
        src_lbl.text(0, 0.5, 'X-ray source',
                     color=p.TEXT, fontsize=9, fontweight='bold', va='center')
        self.ax_src = self.fig.add_axes([0.07, 0.04, 0.12, 0.17])
        self.ax_src.set_facecolor(p.PANEL)
        self.src_labels = [f"{name}  {wl:.4f} A" for name, wl in self.WAVELENGTHS]
        self.radio_src = RadioButtons(
            self.ax_src, self.src_labels, active=0, activecolor=p.ACCENT,
        )
        for lbl in self.radio_src.labels:
            lbl.set_color(p.TEXT); lbl.set_fontsize(8)
        for s in self.ax_src.spines.values():
            s.set_color(p.GRID)
        self.radio_src.on_clicked(self._on_source)

        # Sliders: 2theta_max, h_max, FWHM
        ax_tt = self.fig.add_axes([0.28, 0.175, 0.45, 0.020])
        self.s_tt = Slider(ax_tt, '2theta max (deg)', 20.0, 140.0,
                           valinit=self.two_theta_max, color=p.ACCENT_2)
        ax_hm = self.fig.add_axes([0.28, 0.125, 0.45, 0.020])
        self.s_hm = Slider(ax_hm, 'h_max', 4, 16,
                           valinit=self.h_max, valstep=1, color=p.ACCENT_2)
        ax_fw = self.fig.add_axes([0.28, 0.075, 0.45, 0.020])
        self.s_fw = Slider(ax_fw, 'peak FWHM (deg)', 0.02, 1.0,
                           valinit=self.fwhm, color=p.ACCENT_2)
        for s in (self.s_tt, self.s_hm, self.s_fw):
            s.label.set_color(p.TEXT)
            s.valtext.set_color(p.TEXT)
        self.s_tt.on_changed(self._on_slider)
        self.s_hm.on_changed(self._on_slider)
        self.s_fw.on_changed(self._on_slider_fwhm)

        # Save / Export buttons
        ax_save = self.fig.add_axes([0.80, 0.175, 0.08, 0.035])
        self.btn_save = Button(ax_save, 'Save PNG',
                               color=p.ACCENT, hovercolor='#ff8c5c')
        self.btn_save.label.set_color('white')
        self.btn_save.on_clicked(self._on_save)

        ax_exp = self.fig.add_axes([0.89, 0.175, 0.08, 0.035])
        self.btn_exp = Button(ax_exp, 'Export',
                              color=p.ACCENT_2, hovercolor='#6ba7e8')
        self.btn_exp.label.set_color('white')
        self.btn_exp.on_clicked(self._on_export)

    def _on_source(self, label):
        for name, wl in self.WAVELENGTHS:
            if label.startswith(name):
                self.wavelength = wl
                break
        self.refresh()

    def _on_slider(self, _val):
        self.two_theta_max = float(self.s_tt.val)
        self.h_max = int(self.s_hm.val)
        self.refresh()

    def _on_slider_fwhm(self, val):
        # FWHM affects only the convolved profile, not the peak table.
        self.fwhm = float(val)
        self._draw()

    def _on_close(self, _evt):
        # Tell the parent we're gone so it stops refreshing us.
        if self.parent.powder_win is self:
            self.parent.powder_win = None

    def refresh(self):
        """Recompute peaks from the parent's current calculator and redraw.

        Called on composition/a0/tilt changes in the main window, and on any
        of our own slider changes that alter the peak list."""
        calc = self.parent.calc
        if calc is None:
            self.peaks = []
        else:
            self.peaks = calc.powder_pattern(
                wavelength=self.wavelength,
                two_theta_max=self.two_theta_max,
                h_max=self.h_max,
                I_min=1e-4,
            )
        self._draw()

    def _draw(self):
        p = self.parent
        ax = self.ax
        ax.clear()
        ax.set_facecolor(p.PLOT_BG)
        for s in ax.spines.values():
            s.set_color(p.GRID)
        ax.tick_params(colors=p.TEXT_DIM)
        ax.set_xlabel('2theta (deg)', color=p.TEXT)
        ax.set_ylabel('Intensity (a.u.)', color=p.TEXT)
        ax.set_xlim(0, self.two_theta_max)
        ax.grid(alpha=0.15, color=p.GRID)

        # Title and subtitle (composition + tilt)
        A, B, X = p.symbols
        self.title_txt.set_text(
            f"Powder pattern  .  {A}{B}{X}3   a0 = {p.a0:.4f} A"
        )
        sig = self.fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        self.sub_txt.set_text(
            f"Glazer {p.glazer}  .  omega = "
            f"({p.omega_deg[0]:.1f}, {p.omega_deg[1]:.1f}, {p.omega_deg[2]:.1f}) deg"
            f"   .   lambda = {self.wavelength:.4f} A"
            f"   .   sigma = {sig:.3f} deg"
        )

        if not self.peaks:
            ax.text(0.5, 0.5, '(no peaks in range)',
                    color=p.TEXT_DIM, ha='center', va='center',
                    transform=ax.transAxes, fontsize=11)
            self._draw_peak_list()
            self.fig.canvas.draw_idle()
            return

        two_thetas = np.array([pk['two_theta'] for pk in self.peaks])
        I_norms    = np.array([pk['I_norm']   for pk in self.peaks])

        # Gaussian-convolved continuous profile
        x = np.linspace(0.0, self.two_theta_max, 4000)
        y = np.zeros_like(x)
        sigma = self.fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        two_sig2 = 2.0 * sigma * sigma
        for tt, I in zip(two_thetas, I_norms):
            y += I * np.exp(-((x - tt) ** 2) / two_sig2)
        y_max = y.max() if y.max() > 0 else 1.0
        y /= y_max

        # Stick markers (normalised to same scale as the profile peak max)
        ax.vlines(two_thetas, 0, I_norms, color=p.ACCENT_2, alpha=0.6, lw=1.0)
        ax.plot(x, y, color=p.ACCENT, lw=1.4)
        ax.set_ylim(0, 1.08)

        # Label the top-N strongest peaks with their hkl
        order = np.argsort(-I_norms)
        labelled = set()
        labels_drawn = 0
        for idx in order:
            if labels_drawn >= self.n_labels:
                break
            key = round(two_thetas[idx], 2)
            if key in labelled:
                continue
            labelled.add(key)
            h, k, l = self.peaks[idx]['hkl']
            ax.annotate(
                f"({h}{k}{l})",
                xy=(two_thetas[idx], I_norms[idx]),
                xytext=(0, 6), textcoords='offset points',
                ha='center', va='bottom',
                color=p.TEXT, fontsize=8,
            )
            labels_drawn += 1

        self._draw_peak_list()
        self.fig.canvas.draw_idle()

    def _draw_peak_list(self):
        p = self.parent
        ax = self.list_ax
        ax.clear()
        ax.set_facecolor(p.PANEL)
        for s in ax.spines.values():
            s.set_color(p.GRID)
        ax.set_xticks([]); ax.set_yticks([])

        ax.text(0.04, 0.97, 'Top peaks  (hkl  2th  d  I/I0  m)',
                transform=ax.transAxes,
                color=p.TEXT, fontweight='bold', fontsize=9, va='top')

        if not self.peaks:
            return

        top = sorted(self.peaks, key=lambda pk: -pk['I_norm'])[:22]
        lines = []
        for pk in top:
            h, k, l = pk['hkl']
            lines.append(
                f" {h:>2}{k:>2}{l:>2}  "
                f"{pk['two_theta']:6.2f}  "
                f"{pk['d']:5.3f}  "
                f"{pk['I_norm']:5.3f}  "
                f"{pk['mult']:>2d}"
            )
        ax.text(0.02, 0.91, '\n'.join(lines),
                transform=ax.transAxes,
                color=p.TEXT, family='monospace', fontsize=8, va='top')

    def _on_save(self, _):
        A, B, X = self.parent.symbols
        tag = self.parent.glazer.replace('+', 'p').replace('-', 'm')
        fname = os.path.join(
            HERE, f"powder_{A}{B}{X}3_{tag}_{self.wavelength:.4f}A.png"
        )
        self.fig.savefig(fname, dpi=180, facecolor=self.parent.BG)
        self.parent._set_status(f'Saved {os.path.basename(fname)}', ok=True)

    def _on_export(self, _):
        A, B, X = self.parent.symbols
        tag = self.parent.glazer.replace('+', 'p').replace('-', 'm')
        fname = os.path.join(
            HERE, f"powder_{A}{B}{X}3_{tag}.csv"
        )
        with open(fname, 'w', encoding='utf-8') as fh:
            fh.write('h,k,l,two_theta_deg,d_A,I_norm,multiplicity\n')
            for pk in sorted(self.peaks, key=lambda pk: pk['two_theta']):
                h, k, l = pk['hkl']
                fh.write(
                    f"{h},{k},{l},{pk['two_theta']:.4f},{pk['d']:.5f},"
                    f"{pk['I_norm']:.6f},{pk['mult']}\n"
                )
        self.parent._set_status(f'Exported {os.path.basename(fname)}', ok=True)


# ---------------------------------------------------------------------------
def main():
    """Entry point for `pytilt-gui` / `python -m pytilt_diffraction.simulator`.

    Optional CLI: A B X a0  ->  elements and lattice constant.
    """
    symbols = ('Cs', 'Pb', 'Br')
    a0 = 5.874
    if len(sys.argv) >= 4:
        symbols = (sys.argv[1], sys.argv[2], sys.argv[3])
    if len(sys.argv) >= 5:
        a0 = float(sys.argv[4])

    sim = TiltDiffractionSimulator(symbols=symbols, a0=a0)
    sim.run()


if __name__ == '__main__':
    main()
