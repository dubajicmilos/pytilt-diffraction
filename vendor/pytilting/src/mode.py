# -*- coding: utf-8 -*-
from __future__ import print_function
"""
.. module:: class to define a mode.

.. moduleauthor:: Dawei Wang <dwang5@zoho.com>

This module provides the definition of a particular mode.

"""
import numpy as np
from ase.spacegroup import *
from ase.build import bulk
from ase.io import read,write
import numpy as np
import unittest

import sys
import os
sys.path.append("../src/")

class Mode(object):
    """
    A mode is defined as:
    - A vector to specify ion displacement in each PUC.
    - A Brillouin zone vector to specify its space variation.

    There are some overlap of this module with the class "distortion".

    While the oxygen octaherdon tilting can be treated in the same way as a spceical mode,
    but for now (May 1, 2021), this module in only intended for ion displacements.
    """

    def __init__(self, q_symbol=None, q=None, disp=None):
        if q_symbol is None:
            self.q_symbol = 'Gamma'
        else:
            self.q_symbol = q_symbol

        if q is None:
            self.q = np.array([
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0]
            ])
        else:
            self.q = q

        # Displacements, five atoms, three directions.
        if disp is None:
            self.disp = np.array([
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0]
            ])
        else:
            self.disp = disp


if __name__ == '__main__':
    pass