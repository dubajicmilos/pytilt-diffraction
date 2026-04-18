# -*- coding: utf-8 -*-
from __future__ import print_function

"""
.. module:: Common_structures_for_Perovskites.

.. moduleauthor:: Dawei Wang <dwang5@zoho.com>

This module provides the functions for primitive unit cells of perovskites.

"""
import numpy as np


class Puc:
    """
    This module provides basic class for primitive unit cells (PUC) of perovskites and some basic operations on the PUC level
    to actuate the distortion specified in the "Distortion" class.

    :param symbols: specifies the symbols of atoms in one primitive unit cell.

    :param shift_index: an integer vector indicating the shift from the origin. It is used to indiacte the position of this unit cell in the supercell. 

    """

    def __init__(self, symbols=['A', 'B', 'O'], shift_index=[0, 0, 0]):
        self.symbols = symbols

        self.atoms = [
            {'tag': symbols[0], 'pos': [0.5, 0.5, 0.5]},
            {'tag': symbols[1], 'pos': [0.0, 0.0, 0.0]},
            {'tag': symbols[2], 'pos': [0.5, 0.0, 0.0]},
            {'tag': symbols[2], 'pos': [0.0, 0.5, 0.0]},
            {'tag': symbols[2], 'pos': [0.0, 0.0, 0.5]}
        ]
        self.shift_index = shift_index

    def bravais_basis(self):
        # return absolute position
        p = list(map(lambda x: np.add(x['pos'], self.shift_index),
                     self.atoms))
        return list(map(lambda x: x.tolist(), p))

    def element_basis(self):
        # print atoms symbols
        return list(map(lambda x: x['tag'], self.atoms))

    def rotate_atom(self, atom, omega, covera):
        # define the rotation of the oxygen octahedron
        if atom['tag'] == self.symbols[2]:
            pos = atom['pos']
            a = np.array(omega)
            b = np.array(pos)
            pos1 = np.cross(a, b * np.array([1., 1., covera])) + b
            atom['pos'] = pos1.tolist()
        return atom

    def rotate(self, omega=[0, 0, 0], covera=1.0):
        """
        Only O atoms are rotated.

        :param omega: Rotation angle of oxygen octahedron in **a**-, **b**-, **c**-axis.
        :param covera: The ratio of the length of unit cell along **c**-axis to **a**-axis, which indicates a small distortion of the unit cell .

        """
        self.atoms = list(
            map(lambda x: self.rotate_atom(x, omega, covera), self.atoms)
        )

    def shift(self, u=[0, 0, 0], local_mode=[[0.0, 0.0, 0.0, 0.0, 0.0]]):
        """
        :param u: Magnitudes of the vibration in **a**-, **b**-, **c**-axis.
        :param local_mode: Order is in delta_A, delta_B, delta_O-perp, delta_O-perp, delta_O-parallel

        """
        # dumb but effective method.
        dA = local_mode[0]
        dB = local_mode[1]
        dOp = local_mode[2]
        dOq = local_mode[4]
        # A little math show the following works.

        dpq = np.array([
            [dA, dA, dA],
            [dB, dB, dB],
            [dOq, dOp, dOp],
            [dOp, dOq, dOp],
            [dOp, dOp, dOq]
        ])
        du = np.array([
            [u[0], 0., 0.],
            [0., u[1], 0.],
            [0., 0., u[2]]
        ])

        rslt = np.matmul(dpq, du)
        self.real_shift(rslt)

    def real_shift(self, disp=np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0]])
                   ):
        #disp is a 5*3 matrix
        for i in range(5):
            self.atoms[i]['pos'] += disp[i]
            # print atoms positions

    def print_atoms(self):
        list(map(lambda x: print(x), self.atoms))
