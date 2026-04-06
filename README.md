# NAUTILUS — Ship Hydrostatics in MATLAB

[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b+-orange.svg)](https://mathworks.com)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.19445212-blue.svg)](https://doi.org/10.5281/zenodo.19445212)

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
├── tests/                          # Validation scripts
└── Output/                        # Generated outputs
```

## Theoretical Background

### 1. Volume and Buoyancy Center (Divergence Theorem)

For a closed triangulated mesh, the submerged volume is computed using the divergence theorem:

$$V = \int_V dV = \frac{1}{3} \sum_{f=1}^{n_f} \mathbf{c}_f \cdot \mathbf{n}_f A_f$$

where:
- $n_f$ = number of triangular faces
- $\mathbf{c}_f$ = centroid of face $f$ (mean of its three vertices)
- $\mathbf{n}_f$ = outward unit normal vector of face $f$
- $A_f$ = area of face $f$

The center of buoyancy $\mathbf{B} = (B_x, B_y, B_z)$ is:

$$B_x = \frac{1}{V} \sum_{f=1}^{n_f} c_{f,x} \, (\mathbf{c}_f \cdot \mathbf{n}_f) \, A_f$$
$$B_y = \frac{1}{V} \sum_{f=1}^{n_f} c_{f,y} \, (\mathbf{c}_f \cdot \mathbf{n}_f) \, A_f$$
$$B_z = \frac{1}{V} \sum_{f=1}^{n_f} c_{f,z} \, (\mathbf{c}_f \cdot \mathbf{n}_f) \, A_f$$

**Implementation note:** Only faces below the waterplane ($z < T$) contribute to the volume integral. Face clipping is performed by intersecting each triangle with the $z = T$ plane.

### 2. Waterplane Area and Moments

The waterplane area $A_w$ is computed using the **shoelace formula** on the clipped polygon formed by the intersection of the mesh with the $z = T$ plane:

$$A_w = \frac{1}{2} \sum_{i=1}^{n_w} (x_i y_{i+1} - x_{i+1} y_i)$$

where $n_w$ is the number of vertices in the waterplane polygon.

The waterplane first moments (for calculating center of flotation):

$$M_x = \frac{1}{6} \sum_{i=1}^{n_w} (y_i + y_{i+1}) \, (x_i y_{i+1} - x_{i+1} y_i)$$
$$M_y = \frac{1}{6} \sum_{i=1}^{n_w} (x_i + x_{i+1}) \, (x_i y_{i+1} - x_{i+1} y_i)$$

The center of flotation $F = (F_x, F_y)$ is:

$$F_x = \frac{M_y}{A_w}, \quad F_y = \frac{M_x}{A_w}$$

### 3. Metacentric Radii

The second moments of the waterplane area about the center of flotation:

$$I_x = \frac{1}{12} \sum_{i=1}^{n_w} (y_i + y_{i+1}) \left[ y_i(y_i + y_{i+1}) + y_{i+1}^2 \right] - A_w F_y^2$$
$$I_y = \frac{1}{12} \sum_{i=1}^{n_w} (x_i + x_{i+1}) \left[ x_i(x_i + x_{i+1}) + x_{i+1}^2 \right] - A_w F_x^2$$

The metacentric radii are:

$$KMM_x = \frac{I_x}{V}, \quad KMM_y = \frac{I_y}{V}$$

### 4. Equilibrium Solver

The equilibrium solver finds the draft $T$ and trim $\theta_L$ that satisfy:

$$V(T, \theta_L) = V_{\text{target}}$$
$$M_L(T, \theta_L) = 0$$

where $M_L$ is the longitudinal moment about the centerline.

**Phase 1 — Volume correction:**

$$dT = \frac{\Delta V}{A_w}, \quad T_{\text{new}} = T_{\text{old}} + \alpha \cdot dT$$

where $\alpha = 0.8$ is a relaxation factor and $A_w$ is the waterplane area.

**Phase 2 — Trim correction:**

The longitudinal metacenter $M_L$ is:

$$M_L = B_z + \frac{I_y}{V} \cdot PI_z$$

where $PI$ is the unit vector in the longitudinal direction projected into world coordinates:

$$PI = \left[ -\cos\theta_T \sin\theta_L, \; -\sin\theta_T \cos\theta_L, \; \cos\theta_T \cos\theta_L \right]$$

The trim correction angle (in degrees):

$$\Delta\theta_L = \frac{180}{\pi} \cdot \text{atan2}\left( \frac{(M_L - B)_x \cdot (G - M_L)_z - (M_L - B)_z \cdot (G - M_L)_x}{(M_L - B) \cdot (G - M_L)} \right)$$

The trim is updated with relaxation:

$$\theta_{L,\text{new}} = \theta_{L,\text{old}} + \alpha \cdot \Delta\theta_L$$

The update is clamped to $\pm 15^\circ$ per iteration to prevent divergence.

### 5. Righting Lever (GZ)

The righting lever is the horizontal distance between the center of gravity $G$ and the line of action of the buoyant force through $B$:

$$\mathbf{GZ} = (\mathbf{B} - \mathbf{G}) \times \mathbf{k}$$

where $\mathbf{k} = [0, 0, 1]$ is the vertical unit vector in the world frame.

In scalar form:

$$GZ = \sqrt{(B_x - G_x)^2 + (B_y - G_y)^2}$$

The sign convention: $GZ > 0$ indicates a righting moment (stable), $GZ < 0$ indicates an capsizing moment (unstable).

### 6. Hydrostatic Coefficients

$$\text{KM} = KB + BM = B_z + \frac{I_y}{V}$$

$$\text{KM}_T = KB + BM_T = B_z + \frac{I_x}{V}$$

$$\text{LCB}_x = \frac{M_y}{V} - L_{pp}/2$$

## Validation

| Hull | Property | Analytical | Computed | Error | Status |
|------|----------|------------|----------|-------|--------|
| BoxBarge | Volume (T=1m) | 40.000 m³ | 40.000 m³ | < 0.1% | ✅ |
| BoxBarge | KB (T=1m) | 0.500 m | 0.500 m | < 0.1% | ✅ |
| Wigley | Volume | ~0.248 m³ | ~0.248 m³ | < 0.5% | ✅ |
| Wigley | GZ convergence | 100% | 100% | — | ✅ |

## Citing

If you use NAUTILUS in your research, please cite:

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
