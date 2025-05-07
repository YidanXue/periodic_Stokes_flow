[![made-with-MATLAB](https://img.shields.io/badge/Made%20with-MATLAB_R2023b-orange)]()
[![GitHub latest commit](https://badgen.net/github/last-commit/YidanXue/periodic_Stokes_flow)](https://GitHub.com/YidanXue/periodic_Stokes_flow/commit/)
[![arXiv](https://img.shields.io/badge/DOI-10.1098%2Frspa.2024.0676-blue)](https://doi.org/10.1098/rspa.2024.0676)
![X (formerly Twitter) Follow](https://img.shields.io/twitter/follow/YidanXue?link=https%3A%2F%2Ftwitter.com%2FYidanXue)

Repository for the manuscript "Computing Stokes flows in periodic channels via rational approximation" by Yidan Xue. The arXiv preprint can be found at https://arxiv.org/abs/2407.14864. The PRSA publication can be found at https://doi.org/10.1098/rspa.2024.0676.

Prerequisites
----------------------

The codes need to be run in MATLAB or Octave, its free alternative. The recommended version is MATLAB_R2023b. The user needs to have aaa.m in the same folder (or add it in the default MATLAB path) for the AAA rational approximation. The easist way to do this is to download Chebfun at https://www.chebfun.org. 

MATLAB codes for the paper
----------------------

1) 'poiseuille_smooth.m' computes Stokes flow in a smooth periodic channel driven by a pressure gradient (Poiseuille problem). This code reproduces Fig. 3.
2) 'poiseuille_sharp.m' computes Stokes flow in a periodic channel with sharp corners driven by a pressure gradient (Poiseuille problem). This code reproduces Fig. 4.
3) 'Pozrikidis.m' computes Stokes flow in periodic channels constricted by a moving flat wall and a steady sinusoidal wall. This code reproduces Fig. 5.
4) 'poincare_section.m' computes the Poincar&#233; sections of the unsteady Couette problem between two sinusoidal walls. This code reproduces Fig. 6.
5) 'poincare_case1.m', 'poincare_case2.m' and 'poincare_case3.m' reproduce the 3 unsteady Couette problems presented in Fig. 9.
