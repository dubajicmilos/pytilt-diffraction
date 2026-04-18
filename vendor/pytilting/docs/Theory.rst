Theory
*******

Background
----------

Please see the following papers regardning the theoretical background.

* "Generation of low-symmetry perovskite structures for firsts-principles computation”, by N. Xie, *et al* https://arxiv.org/abs/1910.04507

* A.M. Glazer. "The Classification of Tilted Octahedra in Perovskite", Acta Cryst.: Sect. B, 28, 3384 (1972).

* A.M. Glazer. "Simple ways of determining perovskite structures", Acta Cryst.: Sect. A, 31, 756 (1975).

* W. Zhong, David Vanderbilt, and K. M. Rabe. “First-principles theory of ferroelectric phase transitions for perovskites: The case of BaTiO3”, Phys. Rev. B, 52, 6301 (1995).


Program
---------

The program has two main classes: the *Puc* class and the *Distortion* class stored in two files named "puc.py" and "distortion.py", rsepectively.

The "Puc" class sets the essential information to construct one primitive unit cell while the "Distortion" class can create arbitrary supercells and 
make proper distortion to each of them according to the inputs.

Since processing Glazer notation is a key to the operation of this program, it is separately stored in "glazer.py". 
