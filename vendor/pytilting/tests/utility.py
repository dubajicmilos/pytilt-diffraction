#June 2,  2025: Since ASE's method to find the spacegroup again changed, I have
#   to make changes here accordingly. .

'''
from ase.spacegroup import *

def my_get_spacegroup(*args,**kwargs):
        """
        This function provides a thin wrapper to ASE's **get_spacegroup** function because this function had changed before.
		"""
        # For older version of ase.
        # return get_spacegroup(*args,**kwargs)
        return get_spacegroup(*args)
'''


from ase.spacegroup.symmetrize import check_symmetry
from dataclasses import dataclass

@dataclass
class SG:
    no: int

def my_get_spacegroup(*args,**kwargs):
    symprec = 1e-5
    sg_object = check_symmetry(*args, symprec=symprec, verbose=False)
    return SG(sg_object.number)

