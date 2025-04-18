# Integrability_footprints
Code used in the paper Dissipation-Induced Threshold on Integrability Footprints [arXiv:2504:10255](https://arxiv.org/abs/2504.10255).

## /src
The src directory contains the code responsible for sampling all random ensembles. The file General.fl has functions used commonly in all the .jl files, such as generating a random Kraus map, whereas the others are model-specific.

## /main
The main directory has the code used to generate random matrices (rand_matrices.jl) and compute the respective eigenvalues (eig_from_file.jl). It also contains two notebooks, one where we compute the angular velocity and the other used to plot the data as in the main text. Important: It was not possible to import most of the eigenvalues due to memory restrictions; thus, the Plots notebook will not run successfully (only Figure 2 can be obtained since it does not need eigenvalues explicitly).

## /files
Contains two examples of eigenvalues.

### Organization of the files
We computed the eigenvalues for each case with $\kappa$ ranging from $0.01$ to $1.0$, in steps of 0.01. For each value of $\kappa$ we also computed the eigenvalues for $\kappa \pm d\kappa$, with $d\kappa = \kappa/1000$ (necessary to compute angular velocities and accelerations). Thus, each line corresponds to one of these values. If we only want to know the change in $\kappa$ in steps of $0.01$, we start with the second line and search every three lines.

## /data
Contains the data used in the main text to plot the angular velocities and the schematic.
