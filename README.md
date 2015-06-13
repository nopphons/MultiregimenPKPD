# In silico stochastic process efficacy predictions on gentamicin treatment of neonatal E. Coli infections

forExploratoryCalculations.R is independent of all the other scripts and none of them depend on it. It is included mainly for its very detailed commenting (on model, parameters, variables, functions, etc). It includes copyright and credit information. It is sometimes still used for preliminary estimates and for making figures.

The three main scripts are the following.
   1) outputPKPD.R - Contains three main functions that we used in our calulation: concentration, conccurve, & countfunc. 
   2) outputMultiregimen.R - This code is dependent on outputPKPD.R, which it calls.  It takes a list of many dosing regimens, not just one regimen, as input and produces a data frame for script 3) to use.
   3) correlation.R - This is the main code. It depends on 1) and 2) above, and on two input files: a dosing scheme and an initialization table. It calls everything it needs. It outputs everything we need, as decribed in comments to the script.

The following .csv files are included
   1) miley_input_50doses.csv - 50 input dosing schemes used for the prototypical GA= 40 wk neonate.
   2) initializingTable.csv - This is used as an intialization. It contains the results using 50 dosing schemes from the input above applied to miley (GA = 40 with standard weight). This can be regenerated using outputMultiregimen.R
   3) miley_output_50doses.csv - an example of a table output by correlation.R 