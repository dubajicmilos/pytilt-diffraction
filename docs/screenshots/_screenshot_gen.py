r"""
Headless screenshot generator for the README.

Run from the repo root:
    python docs/screenshots/_screenshot_gen.py

Instantiates TiltDiffractionSimulator under the Agg backend, drives it
through `run()`, and savefig()s the result -- no display required.
"""
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))           # docs/screenshots/
ROOT = os.path.dirname(os.path.dirname(HERE))               # repo root
sys.path.insert(0, ROOT)

from pytilt_diffraction.simulator import TiltDiffractionSimulator


def shoot(out_name, glazer='a+a+a+', omega=(6.0, 6.0, 6.0), zone=(0,0,1),
          layer=0, log_scale=True, twin_3=False, show_labels=False,
          d_min=0.8, h_max=8, spot_scale=250.0, label_thresh=0.05):
    sim = TiltDiffractionSimulator(symbols=('Cs','Pb','I'), a0=6.32)
    sim.glazer = glazer
    sim.omega_deg = list(omega)
    sim.zone_axis = zone
    sim.layer = layer
    sim.log_scale = log_scale
    sim.twin_3 = twin_3
    sim.show_labels = show_labels
    sim.d_min = d_min
    sim.h_max = h_max
    sim.spot_scale = spot_scale
    sim.label_thresh = label_thresh
    # run() builds the full figure (axes, widgets, info panel) then calls
    # plt.show(); under the Agg backend plt.show is a no-op so it returns
    # immediately and we can savefig.
    sim.run()
    out = os.path.join(HERE, out_name)
    sim.fig.savefig(out, dpi=130, facecolor=sim.BG, bbox_inches='tight')
    plt.close(sim.fig)
    print(f"saved: {out}")


if __name__ == '__main__':
    # Same physical state (CsPbI3, a+a+a+, alpha=6 deg, [001], L=0) under
    # linear vs log; plus the half-integer parent layer (L_super=1) where
    # only superlattice peaks survive.
    shoot('linear_after.png',  log_scale=False)
    shoot('log_after.png',     log_scale=True)
    shoot('log_after_L1.png',  log_scale=True,  layer=1)
