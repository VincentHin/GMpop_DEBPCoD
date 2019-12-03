GMpop_DEBPCoD
==========
A individual-based population model for _Globicephala melas_ (pilot whales) based on Dynamic Energy Budgets
----------

Code from: 
---------
Hin, V, Harwood, J. & De Roos, A.M.
_Density-dependence can obscure nonlethal effects of disturbance on life history of medium-sized cetaceans_


Execution of the code requires the installation of EBTtool software package See: (https://staff.fnwi.uva.nl/a.m.deroos/EBT/Software/index.html) for code & installation

Code was executed on 2017 MacBook Pro running macOS High Sierra (10.13.6) with R version 3.6.0. Note that successful execution of the code might depend on your specific installation.

Code is released under GNU General Public License v3.0. A copy of the license is attached.

Please cite the above paper when reusing any part of the model or the code in a publication

Last modified: VH - 3 December, 2019

How to use
----------

* Obtain and install the EBTtool software package from https://staff.fnwi.uva.nl/a.m.deroos/EBT/Software/index.html. Make sure to define the appropriate environmental variables (EBTPATH and PATH) to run the program from the command line (terminal or Terminal tab from within Rstudio).
* In the terminal set the working directory to the folder that contains the model files 
* Invoke `sh run_GMpopDEBPCoD.sh` in terminal to fully clean and rebuild the model and execute the jobs listed in arglist.txt (requires 'GNU parallel')
* This runs the model with the listed parameter files (.cvf):
  
* Run the R-file `GMpop_DEBPCoD.R` to import and plot the model output
* run `make allclean` in terminal to clean all model output files and c object files and executables
