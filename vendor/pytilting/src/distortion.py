# -*- coding: utf-8 -*-
from __future__ import print_function

"""
.. module:: Common_structures_for_Perovskites.

.. moduleauthor:: Dawei Wang <dwang5@zoho.com>

This file takes care of R-3 related structures, which are even more comlpex
than Pnma related phases.

"""
from ase import Atoms
import numpy as np
import math
from ase.data import reference_states as _refstate
from ase.utils import basestring
from functools import reduce

from puc import Puc
from glazer import decode_glazer, my_equal

from copy import deepcopy


# class Pero_factory:
class Distortion:
    """
    The class sets distortion to perovskite structure.

    Basic Information:

    :param symbols: Symbols of ABO3.
    :param lattice_constant: Lattice constant of cubic ABO3.
    :param grid: Number of unit cells along the three axes.
    :param covera: The ratio of the length of unit cell along **c-** axis to **a-** axis, not the supercell.

    Parameters to set up the distortion:

    :param glazer: Glazer notation of distortion.
    :param omega: The angles of the rotation along the three axes.
    :param u: Magnitudes of the vibration in three axes.
    :param k_u: The pattern of shift.
    :param local_mode: Eigenvalues of vibration.
    """

    def __init__(self,
                 system={
                     'symbols': ['Ba', 'Ti', 'O'],
                     'lattice_constant': 4.0,
                     'grid': (2, 2, 2),
                     'covera': 1.0
                 },
                 distort=None,
                 modes=None
                 ):
        """
        Input parameters to set up the desired distortion.
        Default omega is around z axis, with a magnitude of 0.1
        a0 is the lattice constant in the pseudocubic strucutre.
        For now, the compitiblity check falls on the user.

        """

        self.system = system
        self.int_basis = np.diag(self.system['grid'])
        self.basis_factor = 1.0
        self.symbols = self.system['symbols']
        self.covera = self.system['covera']

        # Converts the natural basis back to the crystallographic basis
        self.inverse_matrix = np.linalg.inv(np.transpose(self.int_basis))

        self.pucs = self.generate_pucs()

        a0 = self.system['lattice_constant']
        self.lattice = [[a0, 0, 0], [0, a0, 0], [0, 0, a0 * self.covera]]
        self.cell = np.dot(self.int_basis, self.lattice)

        # set default distort if it is not set.
        if distort is None:
            self._distort = {
                'glazer': 'a0a0a0',
                'omega': (0.0, 0.0, 0.0),
                'u': (0.00, 0.0, 0.0),
                'k_u': [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
                'modes': []
            }

    @property
    def distort(self):
        # Do something if you want
        return self._distort

    @distort.setter
    def distort(self, val):
        # Do something if you want
        self._distort = val
        # Since the distortion is changed, 
        # I need to reset the pucs and the parameters.
        self.pucs = self.generate_pucs()
        self.set_parameters()

    def generate_pucs(self):
        a = self.int_basis[0][0]
        b = self.int_basis[1][1]
        c = self.int_basis[2][2]

        # set up the primitive unit cells.
        lg = []  # undistorted puc
        for i in range(a):
            for j in range(b):
                for k in range(c):
                    puc = Puc(symbols=self.symbols, shift_index=[i, j, k])
                    lg.append(puc)
        return lg

    def get_atoms(self):
        # return atomic information to ASE.
        self.actuate_distort()
        return Atoms(
            symbols=self.element_names(),
            scaled_positions=self.bravais_basis(),
            cell=self.cell,
            pbc=True
        )

    def set_parameters(self):
        # Deal with oxygen tilting first.
        self.omega = self._distort['omega']
        relat2 = [
            my_equal(abs(self.omega[0]), abs(self.omega[1])),
            my_equal(abs(self.omega[1]), abs(self.omega[2])),
            my_equal(abs(self.omega[2]), abs(self.omega[0]))
        ]
        self.glazer = decode_glazer(self._distort['glazer'])
        self.k_omega = self.glazer[0]
        # the angles of rotation should be consistent with the letters of glazer notation
        if (self.glazer[1] != relat2):
            print("Given omega values are not consistent with the given Glazer notation.")
            exit()
        # Deal with ions displacements.
        self.k_u = self._distort['k_u']
        self.local_mode = self._distort['local_mode']
        self.u = self._distort['u']

    def actuate_distort(self):
        if self._distort is not None:
            # Carry out the oxygen octahedron rotation and the displacements.
            self.rotate()
            self.shift()

    def rotate_puc(self, puc):
        """
        Rotate the ions in one primitive unit cell.
        """
        omega1 = [0.0, 0.0, 0.0]
        k_omega2 = self.k_omega / math.pi
        for i in range(3):
            sign = (-1) ** round(np.dot(k_omega2[i], puc.shift_index))
            omega1[i] = sign * self.omega[i]
        puc.rotate(omega1, self.covera)

    def rotate(self):
        # Rotate the ions in all unit cell.
        list(map(lambda x: self.rotate_puc(x), self.pucs))

    def shift_puc(self, puc):
        """
        Move the ions in one primitive unit cell.
        """
        u1 = [0.0, 0.0, 0.0]
        for i in range(3):
            sign = (-1) ** math.floor(np.dot(self.k_u[i]/math.pi, puc.shift_index))
            u1[i] = sign * self.u[i]
        puc.shift(u1, self.local_mode)

        # In addition to such changes, extra mode shall also be applied.
        disp = [0.0, 0.0, 0.0]
        for mode in self._distort['modes']:
            modeq2=mode.q/math.pi
            for i in range(3):
                sign = (-1) ** math.floor(np.dot(modeq2, puc.shift_index))
                disp[i] = sign * mode.disp[i]
            puc.real_shift(np.transpose(disp))

    def shift(self):
        # Move the ions in all unit cell.
        list(map(lambda x: self.shift_puc(x), self.pucs))

    def print(self):
        list(map(lambda x: x.print_atoms(), self.pucs))
        print("---")

    def bravais_basis(self):
        # return absolute position
        p = list(map(lambda x: x.bravais_basis(), self.pucs))
        q = reduce(lambda x, y: x + y, p)
        return list(map(lambda x: np.matmul(self.inverse_matrix, x), q))

    def element_basis(self):
        # match atoms symbols
        p = list(map(lambda x: x.element_basis(), self.pucs))
        p = reduce(lambda x, y: x + y, p)
        d = dict()
        count = 0
        for i in p:
            if i in d:
                continue
            else:
                d[i] = count
                count += 1
        return list(map(lambda x: d[x], p))

    def element_names(self):
        # return atoms symbols
        p = list(map(lambda x: x.element_basis(), self.pucs))
        p = reduce(lambda x, y: x + y, p)
        return p
