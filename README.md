# Inflated-Dynamic-Laplacian
This repository contains Julia code used for the numerical implementation of the trajectory-based version of the inflated Dynamic Laplacian, a mathematical operator used to identify the birth and death of coherent sets within a velocity system over a finite time interval. This trajectory-based method is described in the paper: 

Jason Atnip, Gary Froyland, and Péter Koltai. "An inflated dynamic Laplacian to track the emergence and disappearance of semi-material coherent sets", 2024. https://arxiv.org/abs/2403.10360

The code provided here is for two systems, with definitions and descriptions of both of these systems provided in Sections 6.1 and 6.3 (respectively) of the above paper. The two systems in question are:

1. The Switching Double Gyre, a flow system involving two counter-rotating gyres of unequal size within a rectangular flow domain which switch sides rapidly at the halfway point of the system's runtime; and
2. A coherent-mixing-coherent (2-3-2) Double Gyre, another flow system with two counter-rotating gyres, this time of equal size, which undergoes two regime changes: one from a steady, coherent regime where the two rectangular gyre chambers are invariant to a regime where fluid is mixing between these gyre chambers, and then another back to the coherent regime again. Each flow regime lasts for an equal amount of time.

Included within this repository is a Julia function file used for the numerical construction of the inflated Dynamic Laplacian, separate script files used to implement the method on each version of the Double Gyre system, and a function used to implement the SEBA algorithm on eigenvectors of the inflated Dynamic Laplacian to isolate individual finite-time coherent sets encoded within these eigenvectors.

# Downloading the Repository and Running the Code

