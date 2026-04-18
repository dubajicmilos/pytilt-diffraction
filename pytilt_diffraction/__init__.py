"""pytilt_diffraction

Interactive single-crystal diffraction simulator for Glazer-tilted
perovskites. Wraps the vendored pytilting tilt generator with a
numpy/matplotlib structure-factor calculator.
"""

from pytilt_diffraction.calculator import (
    CIFParser,
    DiffractionCalculator,
)

__all__ = ["CIFParser", "DiffractionCalculator"]
__version__ = "0.1.0"
