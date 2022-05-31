# OMM_GMX
OpenMM TI using Gromacs input

## How to use this example?
ghost.top: the topology file used in Gromacs thermodynamic integration.  
solv_ions.gro: the coonfiguration in Gromacs thermodynamic integration, with the very last two particle Na+ and Cl- as the alchemical particles.  
ti.py: the OpenMM codes that do thermodynamic integration. ParmEd is used to read .top and .gro files. You need to edit the path of your gromacs top library.  
ti.txt: Output. The derivatives dG/dlambda in each iteration and the delta G.  

## How does ParmEd read OPLS-AA force field .top file?
Comb_rule = 3 for OPLS-AA force field .top. ParmEd generates two force objects, one "NonBondedForce" object with all coulomb parameters and the other "CustomNonbondedForce" object with all L-J parameters. In the "NonBondedForce" object, "epsilon" and "sigma" for all the atoms have the same values, and are not correct. Instead, the "CustomNonbondedForce" object keeps all the correct parameters.  
If you don't use OPLS-AA force field, then you need to delete the "CustomNonbondedForce" object CustomLJbyParmed and read the parameters from the "NonBondedForce" object directly.  

## How does ParmEd "translate" L-J parameters?
The L-J parameters in the "CustomNonbondedForce" object are not the same as the values in the OPLS-AA forcefield, which means ParmEd "translates" the parameters and stores them in its own comb_rule. Not very sure about the exact translation rule though. See below how ParmEd define the "CustomNonbondedForce" object.

NB_cus.getEnergyFunction()  

$$\begin{eqnarray} 
V_{LJ}(r) \ &=& \epsilon * (sigr6^2 - sigr6);      \nonumber \\
sigr6 \ &=& sigr2^3; \nonumber \\
sigr2 \ &=& (sigc/r)^2; \nonumber \\
sigc \ &=& \sigma_1 * \sigma_2 \nonumber \\
\epsilon \ &=& \epsilon_1 * \epsilon_2 \nonumber \\
\end{eqnarray}$$

## soft-core potential in TI
Based on the description of soft-core potential in Gromacs manual https://manual.gromacs.org/2019.1/reference-manual/functions/free-energy-interactions.html, the soft-core potential has been modified using the function form of the "CustomNonbondedForce" object that ParmEd generates.

## Amber force field example
If Amber force field is used in Gromacs, the form of custom nonbonded force we create must be modified. The ParmEd will store the epsilon and sigma parameters that we have in amber forcefield.itp, and won't translate them as it does for OPLS-AA force field. So that the van der Waals interaction will be simply

$$\begin{eqnarray}
V_{LJ}(r) \ &=& 4 * \epsilon * (sigr6^2 - sigr6);      \nonumber \\
sigr6 \ &=& sigr2^3; \nonumber \\
sigr2 \ &=& (sigc/r)^2; \nonumber \\
sigc \ &=& 0.5 * (\sigma_1 + \sigma_2) \nonumber \\
\epsilon \ &=& \sqrt{\epsilon_1 * \epsilon_2} \nonumber \\
\end{eqnarray}$$

And the soft-core potential should also be modified accordingly.
