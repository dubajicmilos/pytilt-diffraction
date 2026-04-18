"""
.. module:: tests_1

.. moduleauthor:: Dawei Wang <dwang5@zoho.com>, Na Xie <whereasxn@163.com>


Fri Sep 17 09:31:13 CST 2021
For PZO's complex phases.
This file tests the program with a 4x4x2 supercell, smaller one may not be possible for all the phases listed below.
Note this one also uses the newer version of the program, which provides more distortion modes, but may not necessary.

"""
from ase.spacegroup import *
from ase.build import bulk
from ase.io import read, write
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
        'grid': (8, 8, 8),
        'covera': 1.0
    }
)

dirName = 'structure'
try:
    # Create target Directory
    os.mkdir(dirName)
    print("Directory ", dirName, " Created ")
except FileExistsError:
    print("Directory ", dirName, " already exists")

"I thought the s4 mode is necessary, but it turns out not."
'''
    s4 = Mode(q_symbol='s4',
          q=2 * math.pi * np.array([
              [1.0 / 2, 1.0 / 2, 1.0 / 2],
              [1.0 / 2, 1.0 / 2, 1.0 / 2],
              [1.0 / 2, 1.0 / 2, 1.0 / 2]
          ]),
          disp=np.array([
              [0.01, 0.0, 0.0, 0.0, 0.0],
              [-0.01, 0.0, 0.0, 0.0, 0.0],
              [0.00, 0.0, 0.0, 0.0, 0.0]
          ])
          )
'''

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
        
        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. RevB(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 55)

    def test_Pm_3m(self):
        q_sigma = 2*math.pi*np.array([0.0, 0.0, 0.0])
        s.distort = {
            'glazer': 'a0a0a0',
            'omega': (0.00, 0.00, 0.00),
            'u': (0.0, 0.0, 0.0),
            'k_u': [q_sigma, q_sigma, q_sigma],
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. RevB(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 221)

    def test_I4_mcm(self):
        q_sigma = 2*math.pi*np.array([0.0, 0.0, 0.0])
        s.distort = {
            'glazer': 'a-b0b0',
            'omega': (0.10, 0.00, 0.00),
            'u': (0.0, 0.0, 0.0),
            'k_u': [q_sigma, q_sigma, q_sigma],
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. RevB(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 140)

    def test_Imcm(self):
        q_sigma = 2*math.pi*np.array([1.0/2, 1.0/2, 1.0/2])
        s.distort = {
            'glazer': 'a-a-c0',
            'omega': (0.10, 0.10, 0.00),
            'u': (0.1, -0.1, 0.0),
            'k_u': [q_sigma, q_sigma, q_sigma],
            'local_mode': [0.10, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. RevB(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 74)

    def test_R_3c(self):
        q_sigma = 2*math.pi*np.array([0.0, 0.0, 0.0])
        s.distort = {
            'glazer': 'a-a-a-',
            'omega': (0.10, 0.10, 0.10),
            'u': (0.0, 0.0, 0.0),
            'k_u': [q_sigma, q_sigma, q_sigma],
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. RevB(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 167)
 
    def test_Pnma(self):
        q_sigma = 2*math.pi*np.array([0.0, 0.0, 1.0/2])
        s.distort = {
            'glazer': 'a-a-c+',
            'omega': (0.10, 0.10, 0.12),
            'u': (0.1, 0.1, 0.0),
            'k_u': [q_sigma, q_sigma, q_sigma],
            'local_mode': [0.10, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. RevB(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 62)

    def test_P4mm(self):
        q_sigma = 2*math.pi*np.array([0.0, 0.0, 0.0])
        s.distort = {
            'glazer': 'a0a0a0',
            'omega': (0.00, 0.00, 0.00),
            'u': (0.1, 0.0, 0.0),
            'k_u': [q_sigma, q_sigma, q_sigma],
            'local_mode': [0.10, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. RevB(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 99)

    def test_Amm2(self):
        q_sigma = 2*math.pi*np.array([0.0, 0.0, 0.0])
        s.distort = {
            'glazer': 'a0a0a0',
            'omega': (0.00, 0.00, 0.00),
            'u': (0.1, 0.1, 0.0),
            'k_u': [q_sigma, q_sigma, q_sigma],
            'local_mode': [0.10, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. RevB(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 38)

    def test_R3m(self):
        q_sigma = 2*math.pi*np.array([0.0, 0.0, 0.0])
        s.distort = {
            'glazer': 'a0a0a0',
            'omega': (0.00, 0.00, 0.00),
            'u': (0.1, 0.1, 0.1),
            'k_u': [q_sigma, q_sigma, q_sigma],
            'local_mode': [0.10, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. RevB(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 160)
                                                                                                
    def test_Ima2(self):
        q_sigma = 2*math.pi*np.array([0.0, 0.0, 0.0])
        s.distort = {
            'glazer': 'a-a-c0',
            'omega': (0.10, 0.10, 0.00),
            'u': (0.1, 0.1, 0.0),
            'k_u': [q_sigma, q_sigma, q_sigma],
            'local_mode': [0.10, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. RevB(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 46)

    def test_R3c(self):
        q_sigma = 2*math.pi*np.array([0.0, 0.0, 0.0])
        s.distort = {
            'glazer': 'a-a-a-',
            'omega': (0.10, 0.10, 0.10),
            'u': (0.1, 0.1, 0.1),
            'k_u': [q_sigma, q_sigma, q_sigma],
            'local_mode': [0.10, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. RevB(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 161)

    def test_Pbmm(self):
        q_sigma = 2*math.pi*np.array([1.0/2, 1.0/2, 0.0])
        s.distort = {
            'glazer': 'a0a0a0',
            'omega': (0.00, 0.00, 0.00),
            'u': (0.1, -0.1, 0.0),
            'k_u': [q_sigma, q_sigma, q_sigma],
            'local_mode': [0.10, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: Phys. RevB(2014) 90 220103: Tab. I.
        self.assertEqual(sg.no, 51)

                                                                                            
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Tests)
    results = unittest.TextTestRunner(verbosity=2).run(suite)
