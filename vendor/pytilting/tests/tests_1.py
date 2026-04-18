"""
.. module:: tests_1

.. moduleauthor:: Dawei Wang <dwang5@zoho.com>, Na Xie <whereasxn@163.com>

This file tests the program by setting up proper distortion for a supercell (2x2x2).
The lattice constant of CsSnI3  6.219 angstrom, A  2x2x2 supercell is sufficient to show the structural change.

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
        'symbols': ['Cs', 'Sn', 'I'],
        'lattice_constant': 6.219,
        'grid': (12, 12, 12),
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


class Tests(unittest.TestCase):
    def test_Immm(self):
        s.distort = {
            'glazer': 'a+b+c+',
            'omega': (0.1, 0.12, 0.17),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/1.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 71)

    def test_Immm2(self):
        s.distort = {
            'glazer': 'a+b+b+',
            'omega': (0.1, 0.12, 0.12),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/2.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 71)

    def test_Im_3(self):
        s.distort = {
            'glazer': 'a+a+a+',
            'omega': (0.1, 0.1, 0.1),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/3.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 204)

    def test_Pmmn(self):
        s.distort = {
            'glazer': 'a+b+c-',
            'omega': (0.1, 0.12, 0.17),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/4.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 59)

    def test_Pmmn2(self):
        s.distort = {
            'glazer': 'a+a+c-',
            'omega': (0.1, 0.1, 0.17),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/5.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 137)

    def test_P42nmc(self):
        s.distort = {
            'glazer': 'a+b+b-',
            'omega': (0.1, 0.12, 0.12),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/6.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 59)

    def test_Pmmn4(self):
        s.distort = {
            'glazer': 'a+a+a-',
            'omega': (0.1, 0.1, 0.1),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/7.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 137)

    def test_P21_2(self):
        s.distort = {
            'glazer': 'a+b-c-',
            'omega': (0.1, 0.12, 0.17),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/8.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 11)

    def test_P21(self):
        s.distort = {
            'glazer': 'a+a-c-',
            'omega': (0.1, 0.1, 0.17),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/9.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 11)

    def test_Pmnb(self):
        s.distort = {
            'glazer': 'a+b-b-',
            'omega': (0.1, 0.12, 0.12),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/10.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 62)

    def test_Pmnb2(self):
        s.distort = {
            'glazer': 'a+a-a-',
            'omega': (0.1, 0.10, 0.10),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/11.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 62)

    def test_F1(self):
        s.distort = {
            'glazer': 'a-b-c-',
            'omega': (0.1, 0.12, 0.17),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/12.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 2)

    def test_I2(self):
        s.distort = {
            'glazer': 'a-b-b-',
            'omega': (0.1, 0.12, 0.12),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/13.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 15)

    def test_R_3c(self):
        s.distort = {
            'glazer': 'a-a-a-',
            'omega': (0.1, 0.1, 0.1),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/14.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 167)

    def test_Immm3(self):
        s.distort = {
            'glazer': 'a0b+c+',
            'omega': (0.0, 0.12, 0.17),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/15.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 71)

    def test_I4_mmm(self):
        s.distort = {
            'glazer': 'a0b+b+',
            'omega': (0.0, 0.12, 0.12),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/16.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 139)

    def test_Bmmb(self):
        s.distort = {
            'glazer': 'a0b+c-',
            'omega': (0.0, 0.12, 0.17),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/17.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 63)

    def test_Bmmb2(self):
        s.distort = {
            'glazer': 'a0b+b-',
            'omega': (0.0, 0.12, 0.12),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/18.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 63)

    def test_F2_m11(self):
        s.distort = {
            'glazer': 'a0b-c-',
            'omega': (0.0, 0.12, 0.17),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/19.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 12)

    def test_Imcm(self):
        s.distort = {
            'glazer': 'a0b-b-',
            'omega': (0.0, 0.12, 0.12),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/19.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 74)

    def test_C4_mmb(self):
        s.distort = {
            'glazer': 'a0a0c+',
            'omega': (0.0, 0.0, 0.12),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/21.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 127)

    def test_I4_mmc(self):
        s.distort = {
            'glazer': 'a0a0c-',
            'omega': (0.0, 0.0, 0.10),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/16.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 140)

    def test_Pm_3m(self):
        s.distort = {
            'glazer': 'a0a0a0',
            'omega': (0.0, 0.0, 0.0),
            'u': (0.00, 0.0, 0.0),
            'k_u': 2*math.pi*np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
            'local_mode': [0.00, 0.00, 0.00, 0.00, 0.00],
            'modes': []
        }

        atoms = s.get_atoms()
        atoms.write('./structure/23.cif')

        sg = my_get_spacegroup(atoms, method='spglib')
        print(sg.no)
        # Ref: J. Phys.: Condens. Matter 25 (2013) 175902: Tab. I.
        self.assertEqual(sg.no, 221)


if __name__ == '__main__':
    if not os.path.exists('structure'):
        os.makedirs('structure')
    suite = unittest.TestLoader().loadTestsFromTestCase(Tests)
    results = unittest.TextTestRunner(verbosity=2).run(suite)
