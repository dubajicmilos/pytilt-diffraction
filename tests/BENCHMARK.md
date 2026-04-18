# Benchmark: pytilt-diffraction vs diffsims

Accuracy and speed comparison of `DiffractionCalculator.calculate_structure_factor`
against diffsims 0.7.0 (`get_kinematical_atomic_scattering_factor`).

Reference structure: cubic CsPbBr3, a = 5.874 Å, Pm-3m aristotype, 5 atoms
(Cs/Pb/3x Br), B_iso = 0, occupancy = 1.

## Accuracy (`test_diffsims_parity.py`)

F(hkl) computed with our Cromer-Mann + explicit atomic sum vs. diffsims'
X-ray form factor with the same summation. Agreement is within 1% on every
tested reflection, typically 0.1%.

| hkl       | \|F\| ours | \|F\| ref (diffsims) | rel err |
|-----------|-----------|----------------------|---------|
| (1, 0, 0) | 60.297    | 60.214               | 0.14%   |
| (1, 1, 0) | 93.825    | 93.848               | 0.02%   |
| (1, 1, 1) | 64.276    | 64.277               | 0.002%  |
| (2, 0, 0) | 205.220   | 205.255              | 0.02%   |
| (2, 1, 0) | 53.905    | 54.087               | 0.34%   |
| (2, 1, 1) | 84.762    | 84.710               | 0.06%   |
| (2, 2, 0) | 185.559   | 185.539              | 0.01%   |
| (3, 0, 0) | 50.163    | 50.221               | 0.12%   |
| (3, 1, 0) | 78.684    | 78.627               | 0.07%   |
| (3, 1, 1) | 49.195    | 49.293               | 0.20%   |
| (2, 2, 2) | 172.395   | 172.343              | 0.03%   |

Remaining difference is the published form-factor table (four-Gaussian +
constant Cromer-Mann coefficients in ours; EMsoft-derived table in
diffsims). Both are X-ray form factors; both reduce to Z at s = 0.

### Why not call `get_kinematical_structure_factor` directly

In diffsims 0.7.0, `find_asymmetric_positions` indexes only the first entry
of `corepos`:

```python
return [
    np.array([np.allclose(xyz, asym_xyz) for xyz in positions])
    for asym_xyz in asymmetric_positions
][0]   # <-- drops all but the first asymmetric position
```

For CsPbBr3 this means only the Cs contribution reaches the sum; Pb and Br
are silently dropped. The parity test therefore compares form-factor
parameterisations with our own summation, which is the fair comparison.

## Speed

### Per-reflection F(hkl) (9,260 hkl on cubic CsPbBr3, 5 atoms)

| Backend        | time    | per reflection | ratio  |
|----------------|---------|----------------|--------|
| ours           | 0.35 s  | 37 μs          | 1.00x  |
| diffsims form factors + our summation | 1.08 s | 117 μs | 3.1x slower |

### Full GUI update path (rebuild structure + pattern scan)

2x2x2 supercell (40 atoms), `a+b-b-` tilt pattern, `[001]` zone, h_max = 10:

| Stage                                 | time    |
|---------------------------------------|---------|
| `rebuild_structure` (pytilting + CIF) | 80 ms   |
| `recompute_pattern` (F(hkl) scan)     | 70 ms   |
| **total per slider tick**             | **150 ms** |
| **update rate**                       | **6.7 Hz** |

The pattern feels live on small supercells. On larger distortions or bigger
tilt grids the pytilting rebuild dominates; the diffraction scan stays
tractable.

## Conclusion

- **Accuracy**: agreement well under 1% vs. an independent X-ray form-factor
  reference.
- **Speed**: our calculator is ~3x faster per reflection than diffsims'
  form-factor function wrapped with the same summation. (Direct comparison
  against diffsims' top-level structure-factor function is blocked by the
  diffsims bug above.)
- **Recommendation**: keep ours. The Cromer-Mann / explicit-sum path is
  faster, numerically compatible, and has no dependency on diffsims' buggy
  asymmetric-unit resolution.

## Reproduce

```bash
pip install diffsims pytest
cd pytilt-diffraction
python -m pytest tests/test_diffsims_parity.py -v
python -m pytest tests/test_physics.py -v
```

Hardware: measurements on Windows 11, Python 3.14.2. Absolute timings vary
across machines; the 3x speed ratio is reproducible.
