# RWSP-VM
Random-walk shielding-potential viscosity model (RWSP-VM) is used to calculate the viscosities for elemental metals in warm dense state. 

% --------------Copyright(c)--------------

Version: 3.3

Created by: Yuqing Cheng, ORCID: 0000-0002-3371-9164

Created date: March 10, 2020

Last modified date: January 10, 2022

% -----------------------------------------

These are the codes for Random-Walk Shielding-Potential Vicosity Model, which can be run through Matlab.

This program contains three files: RWSPVM.m, Z_TF.m and Inti2ii.m.

------- RWSPVM.m -------

This is the main program.

There are serveral examples that can be used (Al, Fe, U, Be).

Temperature range can be changed when necessary. 	


------- Z_TF.m -------

This is the function that calculates the average ionization as a function of temperature.


------- Inti2ii.m -------

This is the function that calculates the integral term (Equation 7) as a function of temperature.

------- Implementation -------

To obtain the viscosity data of a certain metal, follow the steps:

1. Set the right quantities of the metal in RWSPVM.m, including Z (nuclear charge number), Am (atomic weight) and Rho (mass density).
2. Change the temperature range if necessary
3. run RWSPVM.m

% -----------------------------------------

