# OMM_GMX
OpenMM TI using Gromacs input

## How does ParmEd read OPLS-AA force field .top file?
Comb_rule = 3 for OPLS-AA force field .top. ParmEd generates two force objects, one "NonBondedForce" object with all coulomb parameters and the other "CustomNonbondedForce" object with all L-J parameters. In the "NonBondedForce" object, "epsilon" and "sigma" for all the atoms have the same values, and are not correct. Instead, the "CustomNonbondedForce" object keeps all the correct parameters.

## How does ParmEd "translate" L-J parameters?
The L-J parameters in the "CustomNonbondedForce" object are not the same as the values in the OPLS-AA forcefield, which means ParmEd "translates" the parameters and stores them in its own comb_rule. Not very sure about the exact translation rule though. See below how ParmEd define the "CustomNonbondedForce" object.

NB_cus.getEnergyFunction()  

epsilon_1 * epsilon_2 * (sigr6^2 - sigr6);  
sigr6 = sigr2 * sigr2 * sigr2;  
sigr2 = (sigc/r)^2;  
sigc= sigma_1 * sigma_2  

## soft-core potential in TI
Based on the description of soft-core potential in Gromacs manual https://manual.gromacs.org/2019.1/reference-manual/functions/free-energy-interactions.html, the soft-core potential has been modified using the function form of the "CustomNonbondedForce" object that ParmEd generates.


