"""
SINGLE CRYSTAL DIFFRACTION SIMULATOR
=====================================

Interactive GUI for simulating single crystal X-ray diffraction patterns.

Features:
- Load any CIF file
- Select diffraction plane (zone axis)
- Adjustable d-spacing limits
- Spot sizes proportional to |F|² (intensity)
- Miller index labels
- Wavelength selection
- Export pattern as image

Similar to CrystalMaker's SingleCrystal software.

Requirements:
    pip install matplotlib numpy

Usage:
    python single_crystal_diffraction.py [cif_file]
    
Or in Python:
    from single_crystal_diffraction import SingleCrystalSimulator
    sim = SingleCrystalSimulator("structure.cif")
    sim.run()
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons
from matplotlib.patches import Circle
import re
import os
import sys

# =============================================================================
# ATOMIC SCATTERING FACTORS (Cromer-Mann coefficients)
# =============================================================================

SCATTERING_FACTORS = {
    'H':  ([0.493, 0.323, 0.140, 0.041], [10.511, 26.126, 3.142, 57.800], 0.003),
    'He': ([0.873, 0.631, 0.311, 0.178], [9.104, 3.357, 22.928, 0.982], 0.006),
    'Li': ([1.128, 0.751, 0.618, 0.465], [3.955, 1.052, 85.391, 168.261], 0.038),
    'Be': ([1.592, 1.128, 0.539, 0.703], [43.643, 1.862, 103.483, 0.542], 0.038),
    'B':  ([2.055, 1.333, 1.098, 0.707], [23.219, 1.021, 60.350, 0.140], -0.193),
    'C':  ([2.310, 1.020, 1.589, 0.865], [20.844, 10.208, 0.569, 51.651], 0.216),
    'N':  ([12.213, 3.132, 2.013, 1.166], [0.006, 9.893, 28.997, 0.583], -11.529),
    'O':  ([3.049, 2.287, 1.546, 0.867], [13.277, 5.701, 0.324, 32.909], 0.251),
    'F':  ([3.539, 2.641, 1.517, 1.024], [10.283, 4.294, 0.262, 26.148], 0.278),
    'Ne': ([3.955, 3.112, 1.455, 1.125], [8.404, 3.426, 0.231, 21.718], 0.352),
    'Na': ([4.763, 3.174, 1.267, 1.113], [3.285, 8.842, 0.314, 129.424], 0.676),
    'Mg': ([5.420, 2.174, 1.227, 2.307], [2.828, 79.261, 0.381, 7.194], 0.858),
    'Al': ([6.420, 1.900, 1.594, 1.965], [3.039, 0.743, 31.547, 85.089], 1.115),
    'Si': ([6.292, 3.035, 1.989, 1.541], [2.439, 32.334, 0.678, 81.694], 1.141),
    'P':  ([6.435, 4.179, 1.780, 1.491], [1.907, 27.157, 0.526, 68.164], 1.115),
    'S':  ([6.905, 5.203, 1.438, 1.586], [1.468, 22.215, 0.254, 56.172], 0.867),
    'Cl': ([11.460, 7.196, 6.256, 1.645], [0.010, 1.166, 18.519, 47.778], -9.557),
    'Ar': ([7.484, 6.772, 0.654, 1.644], [0.907, 14.841, 43.898, 33.393], 1.444),
    'K':  ([8.219, 7.440, 1.052, 0.866], [12.795, 0.775, 213.187, 41.684], 1.423),
    'Ca': ([8.627, 7.387, 1.590, 1.021], [10.442, 0.660, 85.748, 178.437], 1.375),
    'Sc': ([9.189, 7.368, 1.641, 1.468], [9.021, 0.573, 136.108, 51.353], 1.333),
    'Ti': ([9.759, 7.356, 1.699, 1.902], [7.851, 0.500, 35.634, 116.105], 1.281),
    'V':  ([10.297, 7.351, 2.070, 2.057], [6.866, 0.438, 26.894, 102.478], 1.220),
    'Cr': ([10.640, 7.354, 3.324, 1.492], [6.104, 0.392, 20.263, 98.740], 1.183),
    'Mn': ([11.282, 7.357, 3.019, 2.244], [5.341, 0.343, 17.867, 83.754], 1.090),
    'Fe': ([11.770, 7.357, 3.522, 2.305], [4.761, 0.307, 15.354, 76.881], 1.037),
    'Co': ([12.284, 7.341, 4.003, 2.349], [4.279, 0.278, 13.536, 71.169], 1.012),
    'Ni': ([12.838, 7.292, 4.444, 2.380], [3.878, 0.257, 12.176, 66.342], 1.034),
    'Cu': ([13.338, 7.168, 5.616, 1.673], [3.583, 0.247, 11.397, 64.813], 1.191),
    'Zn': ([14.074, 7.032, 5.165, 2.410], [3.266, 0.233, 10.316, 58.710], 1.304),
    'Ga': ([15.235, 6.701, 4.359, 2.962], [3.067, 0.241, 10.781, 61.414], 1.719),
    'Ge': ([16.082, 6.375, 3.707, 3.683], [2.851, 0.252, 11.447, 54.763], 2.131),
    'As': ([16.672, 6.070, 3.431, 4.278], [2.635, 0.265, 12.948, 47.797], 2.531),
    'Se': ([17.001, 5.820, 3.973, 4.354], [2.410, 0.273, 15.237, 43.816], 2.841),
    'Br': ([17.179, 5.236, 5.638, 3.985], [2.172, 16.580, 0.261, 41.433], 2.956),
    'Kr': ([17.355, 6.729, 5.549, 3.537], [1.938, 16.562, 0.226, 39.397], 2.825),
    'Rb': ([17.178, 9.644, 5.140, 1.529], [1.789, 17.315, 0.275, 164.934], 3.487),
    'Sr': ([17.566, 9.818, 5.422, 2.669], [1.556, 14.099, 0.166, 132.376], 2.506),
    'Y':  ([17.776, 10.295, 5.726, 3.266], [1.403, 12.801, 0.261, 104.354], 1.913),
    'Zr': ([17.876, 10.948, 5.417, 3.657], [1.276, 11.916, 0.118, 87.663], 2.069),
    'Nb': ([17.614, 12.014, 4.042, 3.533], [1.189, 11.766, 0.205, 69.796], 3.756),
    'Mo': ([3.703, 17.236, 12.888, 3.743], [0.277, 1.096, 11.004, 61.658], 4.387),
    'Ru': ([19.267, 12.918, 4.863, 1.568], [0.809, 8.435, 24.800, 94.293], 5.379),
    'Rh': ([19.296, 14.350, 4.734, 1.289], [0.752, 8.218, 25.875, 98.606], 5.328),
    'Pd': ([19.332, 15.502, 5.295, 0.606], [0.699, 7.989, 25.205, 76.899], 5.266),
    'Ag': ([19.281, 17.266, 4.689, 1.091], [0.645, 7.474, 24.661, 99.816], 5.179),
    'Cd': ([19.221, 17.644, 4.461, 1.603], [0.595, 6.909, 24.701, 87.482], 5.069),
    'In': ([19.162, 18.560, 4.295, 2.040], [0.548, 6.378, 25.850, 92.803], 4.939),
    'Sn': ([19.189, 19.101, 4.458, 2.466], [5.830, 0.503, 26.891, 83.957], 4.782),
    'Sb': ([19.642, 19.045, 5.037, 2.683], [5.303, 0.461, 27.907, 75.283], 4.591),
    'Te': ([19.964, 19.014, 6.145, 2.524], [4.817, 0.421, 28.528, 70.840], 4.352),
    'I':  ([20.147, 18.995, 7.514, 2.273], [4.347, 0.381, 27.766, 66.878], 4.071),
    'Xe': ([20.293, 19.030, 8.977, 1.990], [3.928, 0.344, 26.466, 64.266], 3.712),
    'Cs': ([20.389, 19.106, 10.662, 1.495], [3.569, 0.311, 24.388, 213.904], 3.335),
    'Ba': ([20.336, 19.297, 10.888, 2.696], [3.216, 0.276, 20.207, 167.202], 2.773),
    'La': ([20.578, 19.599, 11.373, 3.287], [2.948, 0.244, 18.773, 133.124], 2.147),
    'Ce': ([21.167, 19.769, 11.851, 3.330], [2.812, 0.227, 17.608, 127.113], 1.863),
    'Pr': ([22.044, 19.670, 12.386, 2.824], [2.774, 0.222, 16.767, 143.644], 2.058),
    'Nd': ([22.684, 19.685, 12.774, 2.851], [2.662, 0.211, 15.885, 137.903], 1.991),
    'Sm': ([24.004, 19.426, 13.440, 2.896], [2.473, 0.196, 14.400, 128.007], 2.209),
    'Eu': ([24.627, 19.089, 13.760, 2.923], [2.388, 0.190, 13.754, 123.174], 2.575),
    'Gd': ([25.071, 19.080, 13.852, 3.545], [2.253, 0.181, 12.933, 101.398], 2.419),
    'Tb': ([25.898, 18.219, 14.317, 2.954], [2.243, 0.196, 12.665, 115.362], 3.583),
    'Dy': ([26.507, 17.638, 14.560, 2.965], [2.180, 0.202, 12.190, 111.874], 4.297),
    'Ho': ([26.905, 17.294, 14.558, 3.638], [2.071, 0.197, 11.441, 92.657], 4.567),
    'Er': ([27.656, 16.428, 14.978, 2.982], [2.074, 0.223, 11.361, 105.703], 5.920),
    'Yb': ([28.664, 15.434, 15.309, 2.990], [1.989, 0.257, 10.665, 100.417], 7.567),
    'Lu': ([28.948, 15.221, 15.100, 3.716], [1.902, 0.261, 9.985, 84.330], 7.976),
    'Hf': ([29.144, 15.173, 14.759, 4.300], [1.833, 0.275, 9.371, 72.029], 8.582),
    'Ta': ([29.202, 15.229, 14.514, 4.765], [1.773, 0.295, 9.370, 63.364], 9.244),
    'W':  ([29.082, 15.430, 14.433, 5.120], [1.720, 0.321, 9.243, 57.056], 9.888),
    'Re': ([28.762, 15.719, 14.556, 5.442], [1.672, 0.351, 9.092, 52.086], 10.472),
    'Os': ([28.189, 16.155, 14.931, 5.675], [1.629, 0.387, 8.979, 48.165], 11.000),
    'Ir': ([27.305, 16.730, 15.611, 5.834], [1.593, 0.417, 8.866, 45.001], 11.472),
    'Pt': ([27.006, 17.764, 15.713, 5.784], [1.513, 0.445, 8.812, 38.610], 11.688),
    'Au': ([16.882, 18.591, 25.558, 5.860], [0.461, 8.622, 1.483, 36.396], 12.066),
    'Hg': ([20.681, 19.042, 21.657, 5.968], [0.545, 8.448, 1.573, 38.325], 12.609),
    'Tl': ([27.545, 19.158, 15.538, 5.525], [0.655, 8.707, 1.963, 45.815], 13.174),
    'Pb': ([31.062, 13.064, 18.442, 5.970], [0.690, 2.358, 8.618, 47.258], 13.412),
    'Bi': ([33.369, 12.951, 16.588, 6.469], [0.704, 2.924, 8.794, 48.009], 13.578),
    'U':  ([36.029, 23.054, 15.143, 4.303], [0.529, 3.263, 16.092, 100.613], 13.446),
}


def get_scattering_factor(element, s):
    """Calculate atomic scattering factor f(s) where s = sin(θ)/λ."""
    base_element = re.sub(r'[0-9+-]', '', element)
    if base_element not in SCATTERING_FACTORS:
        base_element = 'C'
    a, b, c = SCATTERING_FACTORS[base_element]
    a, b = np.array(a), np.array(b)
    return c + sum(a[i] * np.exp(-b[i] * s**2) for i in range(4))


# =============================================================================
# CIF PARSER
# =============================================================================

class CIFParser:
    """Parse CIF files to extract crystal structure."""
    
    def __init__(self, filename):
        self.filename = filename
        self.cell_params = {}
        self.atoms = []
        self.space_group = "P 1"
        self._parse()
    
    def _parse(self):
        with open(self.filename, 'r') as f:
            content = f.read()
        
        # Cell parameters
        patterns = {
            'a': r'_cell_length_a\s+([0-9.]+)',
            'b': r'_cell_length_b\s+([0-9.]+)',
            'c': r'_cell_length_c\s+([0-9.]+)',
            'alpha': r'_cell_angle_alpha\s+([0-9.]+)',
            'beta': r'_cell_angle_beta\s+([0-9.]+)',
            'gamma': r'_cell_angle_gamma\s+([0-9.]+)',
        }
        for param, pattern in patterns.items():
            match = re.search(pattern, content)
            if match:
                self.cell_params[param] = float(match.group(1))
        
        # Space group
        sg_match = re.search(r'_space_group_name_H-M_alt\s+"([^"]+)"', content)
        if sg_match:
            self.space_group = sg_match.group(1)
        
        # Atoms
        self._parse_atoms(content)
    
    def _parse_atoms(self, content):
        atom_loop = re.search(
            r'loop_\s*\n((?:\s*_atom_site_\w+\s*\n)+)((?:(?!\s*loop_|\s*_).*\n?)*)',
            content
        )
        if not atom_loop:
            return
        
        headers = re.findall(r'_atom_site_(\w+)', atom_loop.group(1))
        data = atom_loop.group(2).strip()
        
        for line in data.split('\n'):
            line = line.strip()
            if not line or line.startswith('_') or line.startswith('loop_'):
                continue
            values = line.split()
            if len(values) >= len(headers):
                atom = dict(zip(headers, values))
                self.atoms.append({
                    'symbol': atom.get('type_symbol', atom.get('label', 'X')[:2].rstrip('0123456789')),
                    'label': atom.get('label', ''),
                    'x': float(atom.get('fract_x', 0)),
                    'y': float(atom.get('fract_y', 0)),
                    'z': float(atom.get('fract_z', 0)),
                    'occupancy': float(atom.get('occupancy', 1.0)),
                    'B_iso': float(atom.get('B_iso_or_equiv', 0)) if 'B_iso_or_equiv' in atom else 0,
                })
    
    def get_cell_matrix(self):
        a = self.cell_params['a']
        b = self.cell_params['b']
        c = self.cell_params['c']
        alpha = np.radians(self.cell_params['alpha'])
        beta = np.radians(self.cell_params['beta'])
        gamma = np.radians(self.cell_params['gamma'])
        
        cos_a, cos_b, cos_g = np.cos(alpha), np.cos(beta), np.cos(gamma)
        sin_g = np.sin(gamma)
        v = np.sqrt(1 - cos_a**2 - cos_b**2 - cos_g**2 + 2*cos_a*cos_b*cos_g)
        
        return np.array([
            [a, b*cos_g, c*cos_b],
            [0, b*sin_g, c*(cos_a - cos_b*cos_g)/sin_g],
            [0, 0, c*v/sin_g]
        ])
    
    def get_reciprocal_matrix(self):
        M = self.get_cell_matrix()
        return 2 * np.pi * np.linalg.inv(M).T


# =============================================================================
# DIFFRACTION CALCULATOR
# =============================================================================

class DiffractionCalculator:
    """Calculate diffraction pattern data."""
    
    def __init__(self, cif_parser):
        self.cif = cif_parser
        self.G = self.cif.get_reciprocal_matrix()
        self.M = self.cif.get_cell_matrix()
    
    def calculate_structure_factor(self, h, k, l):
        """Calculate F(hkl) and d-spacing."""
        hkl = np.array([h, k, l])
        Q = self.G @ hkl
        Q_mag = np.linalg.norm(Q)
        
        if Q_mag < 1e-10:
            return 0.0 + 0.0j, np.inf
        
        d = 2 * np.pi / Q_mag
        s = Q_mag / (4 * np.pi)
        
        F = 0.0 + 0.0j
        for atom in self.cif.atoms:
            f = get_scattering_factor(atom['symbol'], s)
            T = np.exp(-atom['B_iso'] * s**2)
            phase = 2 * np.pi * (h*atom['x'] + k*atom['y'] + l*atom['z'])
            F += atom['occupancy'] * f * T * np.exp(1j * phase)
        
        return F, d
    
    def get_plane_reflections(self, zone_axis, h_max=20, d_min=0.5, I_min=0.01,
                              layer=0, I_max_ref=None):
        """
        Get reflections for a specific zone/plane at a given layer.

        zone_axis: tuple like (0, 0, 1) for hk0 plane, (1, 0, 0) for 0kl, etc.
        layer: integer layer offset along zone axis (supercell Miller units).
               layer=0 is the zero-layer (HK0); layer=1 samples the next plane
               above it (e.g. HK1 for zone [001]).
        I_max_ref: optional absolute intensity to normalise against. If given,
               I_norm = I / I_max_ref so intensities on weak (e.g. off-zero)
               layers stay small instead of being rescaled to 1. Pass the
               strongest reflection of the full structure so that layers with
               only numerical-noise peaks (as in cubic at odd supercell L)
               get filtered out by I_min.
        """
        u, v, w = zone_axis
        layer = int(layer)
        reflections = []

        for h in range(-h_max, h_max + 1):
            for k in range(-h_max, h_max + 1):
                for l in range(-h_max, h_max + 1):
                    # Zone law generalised to layer L: hu + kv + lw = L
                    if h*u + k*v + l*w != layer:
                        continue
                    if h == 0 and k == 0 and l == 0:
                        continue

                    F, d = self.calculate_structure_factor(h, k, l)
                    F_mag = np.abs(F)
                    I = F_mag ** 2

                    if d < d_min:
                        continue

                    reflections.append({
                        'h': h, 'k': k, 'l': l,
                        'd': d, 'F': F_mag, 'I': I,
                        'phase': np.degrees(np.angle(F))
                    })

        # Normalize intensities
        if reflections:
            if I_max_ref is not None and I_max_ref > 0:
                I_max = I_max_ref
            else:
                I_max = max(r['I'] for r in reflections)
            for r in reflections:
                r['I_norm'] = r['I'] / I_max if I_max > 0 else 0

            # Filter by minimum intensity
            reflections = [r for r in reflections if r['I_norm'] >= I_min]

        return reflections
    
    def powder_pattern(self, wavelength=1.5406, two_theta_max=90.0,
                       h_max=12, I_min=1e-4):
        """
        Compute a kinematic powder diffraction pattern.

        Enumerates every (hkl) up to h_max, computes |F|^2 and d-spacing,
        applies the unpolarised-X-ray Lorentz-polarisation factor, and bins
        reflections with the same d into single peaks (so Laue-equivalent
        hkl are summed automatically).

        Returns a list of dicts, sorted by 2theta, each containing:
            d          d-spacing (A)
            two_theta  diffraction angle (deg)
            I          LP-corrected intensity (not normalised)
            I_norm     I / max(I) across the whole pattern
            hkl        one representative (h,k,l) tuple
            mult       number of hkl grouped into this peak

        Parameters
        ----------
        wavelength : float
            X-ray wavelength in A. Default 1.5406 (Cu K-alpha_1).
        two_theta_max : float
            Upper 2theta cutoff, in degrees.
        h_max : int
            Miller-index sweep limit.
        I_min : float
            Drop peaks with I_norm < I_min after normalisation.
        """
        by_d = {}
        for h in range(-h_max, h_max + 1):
            for k in range(-h_max, h_max + 1):
                for l in range(-h_max, h_max + 1):
                    if h == 0 and k == 0 and l == 0:
                        continue
                    F, d = self.calculate_structure_factor(h, k, l)
                    if d <= wavelength / 2.0:
                        continue
                    sin_theta = wavelength / (2.0 * d)
                    if sin_theta <= 0 or sin_theta >= 1.0:
                        continue
                    theta = np.arcsin(sin_theta)
                    two_theta_deg = 2.0 * np.degrees(theta)
                    if two_theta_deg > two_theta_max:
                        continue
                    I_F = abs(F) ** 2
                    LP = (1.0 + np.cos(2.0 * theta) ** 2) / (
                        np.sin(theta) ** 2 * np.cos(theta)
                    )
                    I_lp = I_F * LP

                    key = round(d, 4)
                    if key in by_d:
                        by_d[key]['I'] += I_lp
                        by_d[key]['mult'] += 1
                    else:
                        by_d[key] = {
                            'd': d,
                            'two_theta': two_theta_deg,
                            'I': I_lp,
                            'hkl': (h, k, l),
                            'mult': 1,
                        }

        peaks = sorted(by_d.values(), key=lambda p: p['two_theta'])
        I_max = max((p['I'] for p in peaks), default=1.0)
        for p in peaks:
            p['I_norm'] = p['I'] / I_max if I_max > 0 else 0.0
        return [p for p in peaks if p['I_norm'] >= I_min]

    def get_2d_coordinates(self, reflections, zone_axis):
        """
        Convert 3D hkl to 2D plot coordinates based on zone axis.
        """
        u, v, w = zone_axis
        zone = np.array([u, v, w], dtype=float)
        
        # Find two perpendicular vectors
        if abs(u) <= abs(v) and abs(u) <= abs(w):
            perp1 = np.array([0, -w, v], dtype=float)
        elif abs(v) <= abs(w):
            perp1 = np.array([-w, 0, u], dtype=float)
        else:
            perp1 = np.array([-v, u, 0], dtype=float)
        
        if np.linalg.norm(perp1) < 1e-10:
            perp1 = np.array([1, 0, 0], dtype=float)
        
        perp1 = perp1 / np.linalg.norm(perp1)
        perp2 = np.cross(zone, perp1)
        perp2 = perp2 / np.linalg.norm(perp2)
        
        for r in reflections:
            hkl = np.array([r['h'], r['k'], r['l']], dtype=float)
            r['x'] = np.dot(hkl, perp1)
            r['y'] = np.dot(hkl, perp2)
        
        return reflections


# =============================================================================
# INTERACTIVE SIMULATOR
# =============================================================================

class SingleCrystalSimulator:
    """Interactive single crystal diffraction simulator."""
    
    def __init__(self, cif_filename=None):
        self.cif = None
        self.calc = None
        self.reflections = []
        
        # Default parameters
        self.zone_axis = (0, 0, 1)  # hk0 plane
        self.h_max = 16
        self.d_min = 0.5
        self.I_min = 0.001
        self.spot_scale = 300
        self.show_labels = True
        self.label_threshold = 0.05
        
        if cif_filename:
            self.load_cif(cif_filename)
    
    def load_cif(self, filename):
        """Load a CIF file."""
        self.cif = CIFParser(filename)
        self.calc = DiffractionCalculator(self.cif)
        self.filename = filename
        print(f"Loaded: {filename}")
        print(f"  Cell: a={self.cif.cell_params['a']:.3f}, "
              f"b={self.cif.cell_params['b']:.3f}, "
              f"c={self.cif.cell_params['c']:.3f} Å")
        print(f"  Space group: {self.cif.space_group}")
        print(f"  Atoms: {len(self.cif.atoms)}")
    
    def calculate_pattern(self):
        """Calculate diffraction pattern for current settings."""
        if self.calc is None:
            return
        
        self.reflections = self.calc.get_plane_reflections(
            self.zone_axis, self.h_max, self.d_min, self.I_min
        )
        self.reflections = self.calc.get_2d_coordinates(
            self.reflections, self.zone_axis
        )
    
    def run(self):
        """Run the interactive simulator."""
        if self.cif is None:
            print("No CIF file loaded!")
            return
        
        # Calculate initial pattern
        self.calculate_pattern()
        
        # Create figure with dark theme
        plt.style.use('dark_background')
        self.fig = plt.figure(figsize=(14, 9))
        self.fig.patch.set_facecolor('#1a1a2e')
        
        # Main diffraction plot
        self.ax = self.fig.add_axes([0.08, 0.18, 0.55, 0.72])
        self.ax.set_facecolor('#0a0a15')
        
        # Info panel
        self.info_ax = self.fig.add_axes([0.66, 0.50, 0.32, 0.42])
        self.info_ax.set_facecolor('#1a1a2e')
        self.info_ax.axis('off')
        
        # Title
        self.fig.suptitle('Single Crystal X-ray Diffraction Simulator', 
                         fontsize=16, color='white', fontweight='bold', y=0.96)
        
        # Create widgets
        self._create_widgets()
        
        # Initial plot
        self._update_plot()
        self._update_info()
        
        plt.show()
    
    def _create_widgets(self):
        """Create control widgets."""
        # Zone axis selector
        ax_zone_label = self.fig.add_axes([0.68, 0.43, 0.15, 0.03])
        ax_zone_label.axis('off')
        ax_zone_label.text(0, 0.5, 'Zone Axis (View Direction):', color='white', 
                          fontsize=10, fontweight='bold', va='center')
        
        ax_zone = self.fig.add_axes([0.68, 0.32, 0.28, 0.11])
        self.zone_radio = RadioButtons(
            ax_zone, 
            ['[001] → hk0 plane', '[100] → 0kl plane', '[010] → h0l plane', 
             '[110] → hh̄l plane', '[111] → special'],
            active=0,
            activecolor='#e94560'
        )
        ax_zone.set_facecolor('#16213e')
        for label in self.zone_radio.labels:
            label.set_color('white')
            label.set_fontsize(9)
        self.zone_radio.on_clicked(self._on_zone_change)
        
        # Sliders
        slider_color = '#e94560'
        
        # Max hkl slider
        ax_hmax = self.fig.add_axes([0.20, 0.08, 0.35, 0.025])
        self.slider_hmax = Slider(ax_hmax, 'Max hkl index', 5, 25, 
                                  valinit=self.h_max, valstep=1, color=slider_color)
        self.slider_hmax.label.set_color('white')
        self.slider_hmax.valtext.set_color('white')
        self.slider_hmax.on_changed(self._on_hmax_change)
        
        # d_min slider  
        ax_dmin = self.fig.add_axes([0.20, 0.05, 0.35, 0.025])
        self.slider_dmin = Slider(ax_dmin, 'd_min (Å)', 0.3, 2.0, 
                                  valinit=self.d_min, color=slider_color)
        self.slider_dmin.label.set_color('white')
        self.slider_dmin.valtext.set_color('white')
        self.slider_dmin.on_changed(self._on_dmin_change)
        
        # Spot size slider
        ax_spot = self.fig.add_axes([0.20, 0.02, 0.35, 0.025])
        self.slider_spot = Slider(ax_spot, 'Spot size', 50, 800, 
                                  valinit=self.spot_scale, color=slider_color)
        self.slider_spot.label.set_color('white')
        self.slider_spot.valtext.set_color('white')
        self.slider_spot.on_changed(self._on_spot_change)
        
        # Checkboxes
        ax_check = self.fig.add_axes([0.68, 0.24, 0.12, 0.06])
        ax_check.set_facecolor('#16213e')
        self.check = CheckButtons(ax_check, ['Show Labels'], [self.show_labels])
        for label in self.check.labels:
            label.set_color('white')
            label.set_fontsize(9)
        self.check.on_clicked(self._on_check)
        
        # Save button
        ax_save = self.fig.add_axes([0.82, 0.24, 0.12, 0.04])
        self.btn_save = Button(ax_save, 'Save Image', color='#e94560', hovercolor='#ff6b6b')
        self.btn_save.label.set_color('white')
        self.btn_save.on_clicked(self._on_save)
        
        # Export data button
        ax_export = self.fig.add_axes([0.68, 0.18, 0.12, 0.04])
        self.btn_export = Button(ax_export, 'Export Data', color='#4a6fa5', hovercolor='#6b8fc5')
        self.btn_export.label.set_color('white')
        self.btn_export.on_clicked(self._on_export)
    
    def _on_zone_change(self, label):
        """Handle zone axis change."""
        zones = {
            '[001] → hk0 plane': (0, 0, 1),
            '[100] → 0kl plane': (1, 0, 0),
            '[010] → h0l plane': (0, 1, 0),
            '[110] → hh̄l plane': (1, 1, 0),
            '[111] → special': (1, 1, 1),
        }
        self.zone_axis = zones.get(label, (0, 0, 1))
        self.calculate_pattern()
        self._update_plot()
        self._update_info()
    
    def _on_hmax_change(self, val):
        self.h_max = int(val)
        self.calculate_pattern()
        self._update_plot()
        self._update_info()
    
    def _on_dmin_change(self, val):
        self.d_min = val
        self.calculate_pattern()
        self._update_plot()
        self._update_info()
    
    def _on_spot_change(self, val):
        self.spot_scale = val
        self._update_plot()
    
    def _on_check(self, label):
        if 'Label' in label:
            self.show_labels = not self.show_labels
        self._update_plot()
    
    def _on_save(self, event):
        """Save the current pattern as PNG."""
        filename = os.path.splitext(self.filename)[0] + '_diffraction.png'
        
        # Create a clean figure for saving
        fig_save, ax_save = plt.subplots(figsize=(10, 10))
        ax_save.set_facecolor('#0a0a15')
        fig_save.patch.set_facecolor('#1a1a2e')
        
        self._draw_pattern(ax_save, for_save=True)
        
        zone_str = f"[{self.zone_axis[0]}{self.zone_axis[1]}{self.zone_axis[2]}]"
        ax_save.set_title(f'{os.path.basename(self.filename)}\nZone Axis: {zone_str} | {len(self.reflections)} reflections',
                         color='white', fontsize=14, pad=15)
        
        fig_save.tight_layout()
        fig_save.savefig(filename, dpi=200, facecolor=fig_save.get_facecolor(),
                        bbox_inches='tight')
        plt.close(fig_save)
        print(f"✓ Saved: {filename}")
    
    def _on_export(self, event):
        """Export reflection data to file."""
        filename = os.path.splitext(self.filename)[0] + '_reflections.txt'
        
        with open(filename, 'w') as f:
            f.write(f"# Single Crystal Diffraction Data\n")
            f.write(f"# Source: {os.path.basename(self.filename)}\n")
            f.write(f"# Zone axis: [{self.zone_axis[0]},{self.zone_axis[1]},{self.zone_axis[2]}]\n")
            f.write(f"# d_min: {self.d_min:.2f} Å\n")
            f.write(f"#\n")
            f.write(f"# {'h':>4} {'k':>4} {'l':>4} {'d(Å)':>8} {'|F|':>10} {'I':>12} {'I_norm':>8}\n")
            
            for r in sorted(self.reflections, key=lambda x: -x['I']):
                f.write(f"  {r['h']:4d} {r['k']:4d} {r['l']:4d} "
                       f"{r['d']:8.4f} {r['F']:10.2f} {r['I']:12.2f} {r['I_norm']:8.4f}\n")
        
        print(f"✓ Exported: {filename}")
    
    def _draw_pattern(self, ax, for_save=False):
        """Draw the diffraction pattern on given axes."""
        ax.clear()
        ax.set_facecolor('#0a0a15')
        
        if not self.reflections:
            ax.text(0.5, 0.5, 'No reflections found', ha='center', va='center',
                   color='white', fontsize=14, transform=ax.transAxes)
            return
        
        # Extract data
        x = np.array([r['x'] for r in self.reflections])
        y = np.array([r['y'] for r in self.reflections])
        I = np.array([r['I_norm'] for r in self.reflections])
        
        # Spot sizes (sqrt for better visualization of dynamic range)
        sizes = self.spot_scale * np.sqrt(I) + 3
        
        # Plot with colormap
        scatter = ax.scatter(x, y, s=sizes, c=I, cmap='hot', 
                            alpha=0.9, edgecolors='white', linewidths=0.2,
                            vmin=0, vmax=1)
        
        # Central beam (000)
        ax.scatter([0], [0], s=120, c='cyan', marker='o', 
                  edgecolors='white', linewidths=2, zorder=10)
        ax.annotate('000', (0, 0), xytext=(5, 5), textcoords='offset points',
                   fontsize=8, color='cyan', fontweight='bold')
        
        # Labels for reflections
        if self.show_labels:
            for r in self.reflections:
                if r['I_norm'] >= self.label_threshold:
                    label = f"{r['h']},{r['k']},{r['l']}"
                    ax.annotate(label, (r['x'], r['y']), 
                               xytext=(3, 3), textcoords='offset points',
                               fontsize=6, color='#cccccc', alpha=0.8)
        
        # Styling
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.15, color='gray', linestyle='-')
        ax.axhline(y=0, color='#333355', alpha=0.5, linestyle='-', linewidth=0.5)
        ax.axvline(x=0, color='#333355', alpha=0.5, linestyle='-', linewidth=0.5)
        
        # Axis labels
        ax.set_xlabel('← Reciprocal Space →', color='white', fontsize=10)
        ax.set_ylabel('← Reciprocal Space →', color='white', fontsize=10)
        ax.tick_params(colors='white', labelsize=8)
        
        for spine in ax.spines.values():
            spine.set_color('#444466')
        
        # Add colorbar for saved image
        if for_save:
            cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label('Relative Intensity', color='white')
            cbar.ax.yaxis.set_tick_params(color='white')
            plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
    
    def _update_plot(self):
        """Update the main diffraction plot."""
        self._draw_pattern(self.ax)
        
        zone_str = f"[{self.zone_axis[0]}{self.zone_axis[1]}{self.zone_axis[2]}]"
        self.ax.set_title(f'Zone Axis: {zone_str}', color='white', fontsize=11, pad=10)
        
        self.fig.canvas.draw_idle()
    
    def _update_info(self):
        """Update the information panel."""
        self.info_ax.clear()
        self.info_ax.axis('off')
        self.info_ax.set_facecolor('#1a1a2e')
        
        # Build info text
        zone_names = {
            (0,0,1): 'hk0', (1,0,0): '0kl', (0,1,0): 'h0l',
            (1,1,0): 'hh̄l', (1,1,1): 'special'
        }
        plane_name = zone_names.get(self.zone_axis, 'custom')
        
        info = f"""╔══════════════════════════════════╗
║     STRUCTURE INFORMATION        ║
╚══════════════════════════════════╝

  File: {os.path.basename(self.filename)[:25]}
  Space Group: {self.cif.space_group}

  ┌─ Unit Cell ─────────────────┐
  │  a = {self.cif.cell_params['a']:8.4f} Å          │
  │  b = {self.cif.cell_params['b']:8.4f} Å          │
  │  c = {self.cif.cell_params['c']:8.4f} Å          │
  │  α = {self.cif.cell_params['alpha']:6.2f}°            │
  │  β = {self.cif.cell_params['beta']:6.2f}°            │
  │  γ = {self.cif.cell_params['gamma']:6.2f}°            │
  └─────────────────────────────┘

  Atoms in cell: {len(self.cif.atoms)}
  Viewing: {plane_name} plane
  Reflections: {len(self.reflections)}
  Resolution: {self.d_min:.2f} Å"""
        
        self.info_ax.text(0.02, 0.98, info, transform=self.info_ax.transAxes,
                         fontsize=9, color='#88ccff', fontfamily='monospace',
                         verticalalignment='top', linespacing=1.3)
        
        # Top reflections
        if self.reflections:
            sorted_refs = sorted(self.reflections, key=lambda r: -r['I'])[:6]
            ref_text = "\n  ┌─ Strongest Reflections ─────┐"
            for r in sorted_refs:
                ref_text += f"\n  │ ({r['h']:2d},{r['k']:2d},{r['l']:2d})  d={r['d']:.3f}Å  │"
            ref_text += "\n  └─────────────────────────────┘"
            
            self.info_ax.text(0.02, 0.22, ref_text, transform=self.info_ax.transAxes,
                             fontsize=8, color='#ffcc88', fontfamily='monospace',
                             verticalalignment='top', linespacing=1.2)


# =============================================================================
# QUICK PLOT FUNCTION
# =============================================================================

def plot_diffraction(cif_filename, zone_axis=(0,0,1), h_max=16, d_min=0.5, 
                     save=None, show=True):
    """
    Quick function to plot a diffraction pattern.
    
    Parameters:
    -----------
    cif_filename : str
        Path to CIF file
    zone_axis : tuple
        Zone axis (u, v, w) - e.g., (0,0,1) for hk0 plane
    h_max : int
        Maximum Miller index
    d_min : float
        Minimum d-spacing (Å)
    save : str, optional
        Filename to save image
    show : bool
        Display the plot
        
    Returns:
    --------
    reflections : list
        List of reflection dictionaries
    """
    cif = CIFParser(cif_filename)
    calc = DiffractionCalculator(cif)
    
    reflections = calc.get_plane_reflections(zone_axis, h_max, d_min, I_min=0.001)
    reflections = calc.get_2d_coordinates(reflections, zone_axis)
    
    print(f"Structure: {cif_filename}")
    print(f"Zone axis: [{zone_axis[0]},{zone_axis[1]},{zone_axis[2]}]")
    print(f"Reflections: {len(reflections)}")
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_facecolor('#0a0a15')
    fig.patch.set_facecolor('#1a1a2e')
    
    if reflections:
        x = [r['x'] for r in reflections]
        y = [r['y'] for r in reflections]
        I = [r['I_norm'] for r in reflections]
        sizes = [300 * np.sqrt(i) + 5 for i in I]
        
        scatter = ax.scatter(x, y, s=sizes, c=I, cmap='hot', alpha=0.9,
                            edgecolors='white', linewidths=0.3)
        
        # Labels
        for r in reflections:
            if r['I_norm'] >= 0.05:
                ax.annotate(f"{r['h']},{r['k']},{r['l']}", (r['x'], r['y']),
                           xytext=(3, 3), textcoords='offset points',
                           fontsize=6, color='#aaaaaa')
        
        plt.colorbar(scatter, label='Relative Intensity')
    
    ax.scatter([0], [0], s=100, c='cyan', marker='o', edgecolors='white', linewidths=2)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2, color='gray')
    
    zone_str = f"[{zone_axis[0]}{zone_axis[1]}{zone_axis[2]}]"
    ax.set_title(f'Single Crystal Diffraction\n{os.path.basename(cif_filename)} | Zone: {zone_str}',
                color='white', fontsize=12)
    ax.set_xlabel('h', color='white')
    ax.set_ylabel('k', color='white')
    ax.tick_params(colors='white')
    
    if save:
        fig.savefig(save, dpi=150, facecolor=fig.get_facecolor(), bbox_inches='tight')
        print(f"Saved: {save}")
    
    if show:
        plt.show()
    
    return reflections


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    if len(sys.argv) > 1:
        cif_file = sys.argv[1]
    else:
        cif_files = [f for f in os.listdir('.') if f.endswith('.cif')]
        if cif_files:
            cif_file = cif_files[0]
        else:
            print("Usage: python single_crystal_diffraction.py <cif_file>")
            print("\nNo CIF file found in current directory.")
            sys.exit(1)
    
    if not os.path.exists(cif_file):
        print(f"File not found: {cif_file}")
        sys.exit(1)
    
    # Run interactive simulator
    sim = SingleCrystalSimulator(cif_file)
    sim.run()
