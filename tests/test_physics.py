"""
Physics sanity tests for the pytilt-diffraction simulator.

- Cubic aristotype (a0a0a0) must show no superlattice peaks at odd supercell
  L (which correspond to half-integer L in the parent pseudocubic cell).
- Non-trivial tilts (Pnma a+b-b-) must produce R-point superlattice peaks on
  those odd layers.
- Zone-law filter: every reported reflection must satisfy h*u+k*v+l*w = L.
"""

import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, os.pardir))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from pytilt_diffraction.simulator import TiltDiffractionSimulator


@pytest.fixture
def sim():
    return TiltDiffractionSimulator()


def test_cubic_has_no_odd_layer_peaks(sim):
    """Pm-3m aristotype: L=1 in the 2x2x2 supercell is L=0.5 in the parent
    pseudocubic. There are no R-point reflections, so the layer must be
    empty (down to I_min)."""
    sim.glazer = "a0a0a0"
    sim.omega_deg = [0.0, 0.0, 0.0]
    sim.rebuild_structure()

    sim.layer = 0
    sim.recompute_pattern()
    n_zero = len(sim.reflections)
    assert n_zero > 0, "cubic L=0 should have reflections"

    sim.layer = 1
    sim.recompute_pattern()
    assert len(sim.reflections) == 0, (
        f"cubic L=1 should have no peaks, got {len(sim.reflections)}"
    )


def test_pnma_has_odd_layer_peaks(sim):
    """Pnma (a+b-b-) with non-zero tilt must produce superlattice peaks on
    the L=1 layer of the 2x2x2 supercell (the R-point of the parent)."""
    sim.glazer = "a+b-b-"
    sim.omega_deg = [8.0, 6.0, 6.0]
    sim.rebuild_structure()

    sim.layer = 1
    sim.recompute_pattern()
    assert len(sim.reflections) > 0, (
        "Pnma with tilt should produce L=1 superlattice peaks"
    )


def test_zone_law_respected(sim):
    """Every reported reflection must lie on the requested layer:
    h*u + k*v + l*w = layer."""
    sim.glazer = "a+b-b-"
    sim.omega_deg = [8.0, 6.0, 6.0]
    sim.rebuild_structure()

    for layer in (0, 1, 2):
        sim.layer = layer
        sim.recompute_pattern()
        u, v, w = sim.zone_axis
        for r in sim.reflections:
            assert r["h"] * u + r["k"] * v + r["l"] * w == layer, (
                f"zone law violated on layer={layer}: {r}"
            )


def test_odd_layer_intensity_is_weak_fraction_of_main(sim):
    """Superlattice peaks should be a small fraction of the fundamental
    (Bragg) peaks. Typical tilt-induced R-point intensity is a few percent
    of the strongest fundamental, not comparable to it."""
    sim.glazer = "a+b-b-"
    sim.omega_deg = [8.0, 6.0, 6.0]
    sim.rebuild_structure()

    sim.layer = 0
    sim.recompute_pattern()
    I_main = max(r["I"] for r in sim.reflections)

    sim.layer = 1
    sim.recompute_pattern()
    I_odd = max(r["I"] for r in sim.reflections)

    ratio = I_odd / I_main
    assert ratio < 0.30, (
        f"odd-layer peak suspiciously strong: ratio={ratio:.3f}"
    )
    assert ratio > 1e-5, (
        f"odd-layer peak should not be vanishing at 8/6/6 deg tilt: "
        f"ratio={ratio:.3e}"
    )
