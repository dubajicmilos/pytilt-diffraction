r"""
Headless screenshot for the CIF viewer. Generates a CIF on the fly via
the Glazer simulator's tilt machinery (CsPbI3, a+a+a+, alpha = 6 deg),
then opens it with CIFViewer and savefig()s the result.

Run from the repo root:
    python docs/screenshots/_cif_viewer_screenshot.py
"""
import os
import sys
import math
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))                  # docs/screenshots/
ROOT = os.path.dirname(os.path.dirname(HERE))                      # repo root
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'vendor', 'pytilting', 'src'))

from distortion import Distortion                                  # type: ignore
from pytilt_diffraction.cif_viewer import CIFViewer


def make_cspbi3_cif(out_path):
    """Build a CIF for CsPbI3 with a+a+a+ alpha=6deg, 2x2x2 supercell."""
    s = Distortion(system={
        'symbols': ['Cs', 'Pb', 'I'],
        'lattice_constant': 6.32,
        'grid': (2, 2, 2),
        'covera': 1.0,
    })
    a_rad = math.radians(6.0)
    s.distort = {
        'glazer':     'a+a+a+',
        'omega':      (a_rad, a_rad, a_rad),
        'u':          (0.0, 0.0, 0.0),
        'k_u':        2 * math.pi * np.zeros((3, 3)),
        'local_mode': [0.0] * 5,
        'modes':      [],
    }
    s.get_atoms().write(out_path)
    return out_path


def shoot(out_name, **overrides):
    cif = os.path.join(tempfile.gettempdir(), 'pytilt_cif_view_demo.cif')
    make_cspbi3_cif(cif)
    v = CIFViewer(cif)
    for k, val in overrides.items():
        setattr(v, k, val)
    v.run()
    out = os.path.join(HERE, out_name)
    v.fig.savefig(out, dpi=130, facecolor=v.BG, bbox_inches='tight')
    plt.close(v.fig)
    print(f'saved: {out}')


if __name__ == '__main__':
    # Default zone [001], L = 0, log mode.
    shoot('cif_viewer_log.png', log_scale=True, layer=0)
    # Linear mode at the same state.
    shoot('cif_viewer_linear.png', log_scale=False, layer=0)
    # The half-integer (supercell L=1) layer to show the superlattice peaks.
    shoot('cif_viewer_log_L1.png', log_scale=True, layer=1)
    # Demonstrate that the HKL layer slider works for all zone axes,
    # including [110] (h + k = L) and [111] (h + k + l = L).
    shoot('cif_viewer_zone110_L1.png', log_scale=True, layer=1,
          zone_axis=(1, 1, 0))
    shoot('cif_viewer_zone111_L0.png', log_scale=True, layer=0,
          zone_axis=(1, 1, 1))
    shoot('cif_viewer_zone100_L2.png', log_scale=True, layer=2,
          zone_axis=(1, 0, 0))
