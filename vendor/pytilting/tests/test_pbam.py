"""
.. module:: try-pbam for PZO.

.. moduleauthor:: Dawei Wang <dwang5@zoho.com>, Na Xie <whereasxn@163.com>


Fri Sep 17 09:31:13 CST 2021
For PZO's complex phases.
This file tests the program with a 4x4x2 supercell, smaller one may not be possible for all the phases listed below.
Note this one also uses the newer version of the program, which provides more distortion modes, but may not necessary.

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
        'symbols': ['Pb', 'Zr', 'O'],
        'lattice_constant': 4.11,
        'grid': (12, 12, 12),
        'covera': 1.0
    }
)

class Tests(unittest.TestCase):

    def test_Pbam(self):
        q_sigma = 2*math.pi*np.array([1.0/4, 1.0/4, 0])
        s.distort = {
            'glazer': 'a-a-c0',
            'omega': (-0.10, 0.10, 0.00),
            'u': (0.1, -0.1, 0.0),
            'k_u': [q_sigma, q_sigma, q_sigma],
            'local_mode': [0.10, 0.00, 0.00, 0.00, 0.00],
            'modes': []

        }

        atoms = s.get_atoms()
        atoms.write('./pbam.cif')
        
        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. Rev. B(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 55)
    

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Tests)
    results = unittest.TextTestRunner(verbosity=2).run(suite)
