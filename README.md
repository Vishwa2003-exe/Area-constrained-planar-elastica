# Area-Constrained Planar Elastica

This repository contains Python scripts developed as part of a course project for the **'Geometry and Mechanics of Materials'** course. The project focuses on the mathematical and numerical exploration of the equilibrium configurations of a rigid loop constrained by fixed length and enclosed area. The primary goal is to understand the interplay between geometric constraints and bending energy in determining loop configurations.

### Files Description

#### 1. `o3vsE.py`
- **Description**: Computes and plots the relationship between the total angle $\Theta_3$ and the energy $E$ for different values of $\mu$.

#### 2. `VvsK.py`
- **Description**: Analyzes and visualizes the relationship between the potential function $V(K)$ and the curvature parameter $K$.

#### 3. `KdashvsK.py`
- **Description**: Produces plots of $\dot{K}$ (the derivative of the curvature) against $K$.

#### 4. `o3_3d.py`
- **Description**: Extends the analysis to a 3D plot between total angle $\Theta_3$, energy $E$, and $\mu$.

#### 5. `thetai_vs_mu.py`
- **Description**: Explores the relationship between the maximum loop angle at inflection points $\theta_{\text{max}}$ and the Lagrange multiplier $\mu$.

#### 6. `Area-constrained-planar-elastica.pdf`
- **Description**: The paper used as a reference for the numerical simulations.


### Prerequisite
To run these scripts, ensure you have Python installed along with the following libraries:
- NumPy
- Matplotlib
- SciPy

### How to Use
   ```bash
  git clone https://github.com/Vishwa2003-exe/Area-constrained-planar-elastica.git
  cd Area-constrained-planar-elastica
  python script_name.py

