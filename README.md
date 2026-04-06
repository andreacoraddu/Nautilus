# NAUTILUS — Ship Hydrostatics in MATLAB

[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b+-orange.svg)](https://mathworks.com)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-blue.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

> A MATLAB toolkit for ship hydrostatic and stability analysis using mesh-based methods.

## Features

- **Mesh-first architecture** — Work directly with triangulated hull meshes (OBJ, STL, or MAT)
- **Hydrostatics** — Volume, center of buoyancy, waterplane area, metacentric radius, coefficients
- **Stability** — GZ righting arm curves with equilibrium trim/draft solving via Quasi-Newton method
- **Geometry validation** — Automatic checks for degenerate faces, boundary edges, manifold quality
- **Publication-ready output** — PDF/PNG plots with LaTeX labels, DAT/TXT tables, MAT files

## Installation

```bash
git clone https://github.com/andreacoraddu/Nautilus.git
cd Nautilus
```

No additional dependencies required — pure MATLAB.

## Quick Start

```matlab
% Run default analysis (Wigley hull)
main

% Enable GZ stability curves
cfg.computeGZ = true;
cfg.tetaT_vec = 0:5:60;  % heel angles in degrees
pipeline = Pipeline(cfg);
pipeline.runForHull('Wigley_hull_mesh', 0.10:0.01:0.25);
```

## Documentation

- `docs/theory.md` — Mathematical foundations
- `docs/api.md` — Class and function specifications
- `docs/validation.md` — Expected results against analytical solutions

## Project Structure

```
Nautilus/
├── main.m                          # Entry point
├── Classes/
│   ├── Pipeline.m                  # Orchestrator
│   ├── GeometryLoader.m            # Mesh import (OBJ/STL/MAT)
│   ├── GeometryChecker.m           # Sanity checks
│   ├── HydrostaticsEngine.m        # Draft sweep + GZ
│   ├── Algorithms.m                # Core algorithms
│   └── PostProcessor.m             # Export + plotting
├── Hulls/                          # Benchmark hulls + generators
├── docs/                           # Theory and API docs
├── tests/                          # Validation scripts
└── Output/                         # Generated outputs
```

## Validation

| Hull | Property | Error | Status |
|------|----------|-------|--------|
| BoxBarge | Volume | < 0.1% | ✅ |
| BoxBarge | KB | < 0.1% | ✅ |
| Wigley | Volume | < 0.5% | ✅ |
| Wigley | GZ convergence | 100% | ✅ |

## Algorithm Overview

### Volume Calculation (Divergence Theorem)

For a triangulated mesh, the submerged volume is:

$$V = \frac{1}{3} \sum_{f=1}^{n_f} \mathbf{c}_f \cdot \mathbf{n}_f A_f$$

### Waterplane Area (Shoelace Formula)

$$A_w = \frac{1}{2} \sum_{i=1}^{n} (x_i y_{i+1} - x_{i+1} y_i)$$

### Equilibrium Solver

Two-phase Quasi-Newton method:
1. **Volume phase** — Find draft $T$ such that $V(T) = V_{target}$
2. **Trim phase** — Find trim $\theta_L$ such that $M_L = 0$

### Righting Lever (GZ)

$$GZ = \| (\mathbf{G} - \mathbf{B}) \times \mathbf{k} \|$$

## Citing

```bibtex
@article{Coraddu2012WaveProfileStability,
  author    = {A. Coraddu and P. Gualeni and D. Villa},
  title     = {Investigation about wave profile effects on ship stability},
  booktitle = {Sustainable Maritime Transportation and Exploitation of Sea Resources - Proceedings of the 14th International Congress of the International Maritime Association of the Mediterranean, IMAM 2011},
  editor    = {Enrico Rizzuto and Carlos Guedes Soares},
  pages     = {143--149},
  year      = {2012},
  isbn      = {9780415620826},
  doi       = {10.1201/b11810-25},
  language  = {English},
  address   = {Genova, Italy},
  note      = {14th International Congress of the International Maritime Association of the Mediterranean (IMAM 2011), Genova, Italy, 13--16 September 2011}
}

@misc{nautilus2026,
  author = {Andrea Coraddu},
  title = {NAUTILUS — Ship Hydrostatics in MATLAB},
  year = {2026},
  url = {https://github.com/andreacoraddu/Nautilus}
}
```

## License

MIT. See [LICENSE](LICENSE).
