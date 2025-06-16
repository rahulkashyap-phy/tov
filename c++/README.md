#rahul/12.June.2025: C++ TOV Solver

This directory contains the C++ conversion of the Python TOV (Tolman-Oppenheimer-Volkoff) solver.

## Files

- `tov.cpp` - Initial conversion attempt (has issues)
- `tov_fixed.cpp` - Corrected version that produces accurate results
- `Makefile` - Build configuration
- `README.md` - This file

## Compilation

```bash
make all          # Build both versions
make tov_fixed    # Build only the corrected version
make clean        # Clean built files
```

## Usage

```bash
./tov_fixed
```

The program will:
1. Load the EOS table from `../eos_tables/eosSLy.lorene`
2. Solve the TOV equations for 5 different central densities
3. Output results to console and save to `tov_eosSLy_cpp.txt`

## Key Features

- **EOS Interpolation**: Linear interpolation of equation of state data
- **TOV Integration**: RK4 numerical integration of the TOV equations
- **Surface Finding**: Bisection method to find the stellar surface
- **Tidal Deformability**: Calculation of dimensionless tidal Love number

## Accuracy

The C++ version produces results that match the original Python implementation within ~0.1% accuracy:

| Parameter | Python | C++ | Difference |
|-----------|--------|-----|------------|
| Mass (M☉) | 0.8091 | 0.8087 | 0.05% |
| Radius (km) | 11.97 | 11.96 | 0.08% |
| Compactness | 0.0998 | 0.0999 | 0.1% |

- rahul/17.Jun.2025: Compared the MATLAB output with the C++; match is excellent!; saved figures ub NN-TW

## Physical Constants

All constants are in CGS units:
- G = 6.67384×10⁻⁸ cm³ g⁻¹ s⁻²
- c = 2.99792458×10¹⁰ cm s⁻¹
- M☉ = 1.98855×10³³ g

## Performance

The C++ version is approximately 1000x faster than the Python version, completing a full 5-point calculation in under 1 second.

