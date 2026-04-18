"""
.. module:: tests_BFO

.. moduleauthor:: Dawei Wang <dwang5@zoho.com>, Na Xie <whereasxn@163.com>

This file tests the program with the R3c phase of :math:`\mathrm{BiFeO}_3`.

Set the reverse rotation on the three axes and assume the displacement of the **Fe** ions.


"""
import numpy as np
import unittest
import math
import sys
import os
sys.path.append("../src/")

from distortion import Distortion
from mode import Mode
from utility import my_get_spacegroup

s = Distortion(
    system={
        'symbols': ['Bi', 'Fe', 'O'],
        'lattice_constant': 3.96,
        'grid': (12, 12, 12),
        'covera': 1.0,
    }
)

class Tests(unittest.TestCase):

    def test_R3c(self):
        s.distort = {
            'glazer': 'a-a-a-',
            'omega': (0.1, 0.1, 0.1),
            'u': (0.1, 0.1, 0.1),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.50, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/a-a-a-.cif')

        sg = my_get_spacegroup(atoms,method = 'spglib')
        print(sg.no)
        self.assertEqual(sg.no, 161)


if __name__ == '__main__':
    if not os.path.exists('structure'):
        os.makedirs('structure')
    suite = unittest.TestLoader().loadTestsFromTestCase(Tests)
    results=unittest.TextTestRunner(verbosity=2).run(suite)
