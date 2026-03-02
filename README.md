# Inflated-Dynamic-Laplacian
This repository contains Julia code used for the numerical implementation of the trajectory-based version of the inflated Dynamic Laplacian, a mathematical operator used to identify the birth and death of coherent sets within a velocity system over a finite time interval. This trajectory-based method is described in the paper: 

Jason Atnip, Gary Froyland, and Péter Koltai. "An inflated dynamic Laplacian to track the emergence and disappearance of semi-material coherent sets", 2024. https://arxiv.org/abs/2403.10360

The code provided here is for two systems, with definitions and descriptions of both of these systems provided in Sections 6.1 and 6.3 (respectively) of the above paper. The two systems in question are:

1. The Switching Double Gyre, a flow system involving two counter-rotating gyres of unequal size within a rectangular flow domain which switch sides rapidly at the halfway point of the system's runtime; and
2. A coherent-mixing-coherent (2-3-2) Double Gyre, another flow system with two counter-rotating gyres, this time of equal size, which undergoes two regime changes: one from a steady, coherent regime where the two rectangular gyre chambers are invariant to a regime where fluid is mixing between these gyre chambers, and then another back to the coherent regime again. Each flow regime lasts for an equal amount of time.

Included within this repository is a Julia function file used for the numerical construction of the inflated Dynamic Laplacian, separate script and plotting function files used to implement the method on each version of the Double Gyre system, and a function used to implement the SEBA algorithm on eigenvectors of the inflated Dynamic Laplacian to isolate individual finite-time coherent sets encoded within these eigenvectors.

Also included are "Project.toml" and "Manifest.toml" files which detail the packages and version/subdependency information for these packages pertinent to the Julia environment generated for the successful execution of this code.

# Downloading the Repository and Running the Code

Here are the instructions to follow for the successful execution of this code on your system: 

**Using a Stand Alone Julia REPL Window**

1. Download and save the "Inflated-Dynamic-Laplacian" repository to your system.
2. Open a new Julia REPL window and move to the "Inflated-Dynamic-Laplacian" directory.
3. Type "]", followed by the commands "activate ." and "instantiate" to set up the Julia environment for this repository.
4. If you wish to run the script for the switching Double Gyre system, move to the "DG_Switching" subfolder and run the script with the command: include("iDL_script_switching.jl"). Otherwise, to run the script for the 2-3-2 Double Gyre system, move to the "DG_232" subfolder and execute the command: include("iDL_script_232.jl").
5. Results data and images pertinent to the inflated Dynamic Laplacian method will be stored within either the "DG_Switching" or "DG_232" subdirectories, depending on which system you have chosen.

**Using VS Code**

1. Download and save the "Inflated-Dynamic-Laplacian" repository to your system.
2. Open VS Code, and open the "Inflated-Dynamic-Laplacian" folder in your workspace.
3. Start a new Julia REPL in VS Code. Click on "Julia env" at the bottom of your VS Code window, select "(pick a folder)" from the drop down menu appearing at the top of the window, and find and select the "Inflated-Dynamic-Laplacian" folder in the system dialog.
4. Type "]" followed by the commands "activate ." and "instantiate" to complete set up of the Julia environment for this repository.
5. In the VS Code explorer sidebar, left-click either on the "DG_Switching" subfolder and then left-click the "iDL_script_switching.jl" file to execute this method on the switching Double Gyre system, or left-click on the "DG_232" subfolder and then left-click the "iDL_script_232.jl" file to execute this method on the 2-3-2 Double Gyre system. Whichever code you select should open at the right. Click on the right-pointing triangle icon near the top right to run.
6. Results data and images pertinent to the inflated Dynamic Laplacian method will be stored within either the "DG_Switching" or "DG_232" subdirectories, depending on which system you have chosen.