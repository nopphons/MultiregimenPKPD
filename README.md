The codes used for a submitted paper “In Silico Stochastic Process Efficacy Predictions On Gentamicin Treatment Of Neonatal E. Coli Infections.” This submitted paper extends the paper “Radivoyevitch T, Siranart N, Hlatky L, Sachs R (2015) “Stochastic Process Pharmacodynamics: Dose Timing in Neonatal Gentamicin Therapy as an Example. AAPS Journal 17 (2):447-456. doi:10.1208/s12248-014-9715-3.”

forExploratoryCalculations.R is independent of all the other scripts and none of them depend on it. It is included mainly for its very detailed commenting (on model, parameters, variables, functions, etc). It includes copyright and credit information. It is sometimes still used for preliminary estimates and for making figures.

The three main scripts are the following.
   1) outputPKPD.R - Contains three main functions that we used in our calculations: concentration, conccurve, & countfunc. 
   2) outputMultiregimen.R - This code is dependent on outputPKPD.R, which it calls.  It takes a list of many dosing regimens, not just one regimen, as input and produces a data frame for script 3) to use.
   3) correlation.R - This is the main code. It depends on 1) and 2) above, and on two input files: a dosing scheme and an initialization table. It calls everything it needs. It outputs almost everything we need, as decribed in comments to the script.

The following .csv files are included
   1) miley_input_56_doses.csv - 56 input dosing schemes used for the prototypical GA= 40 wk neonate.
   2) initializingTable.csv - This is used as an intialization. It contains the results using 56 dosing schemes from the input above applied to miley (GA = 40 with standard weight). This can be regenerated using outputMultiregimen.R
   3) miley_output_56_doses.csv - an example of a table output by correlation.R 

Auxiliary subdirectory contains scripts that are used to compare analytic results with Monte Carlo results by considering a virtual neonate which has only about 50 bacteria at t=0 and only give 1, very small, dose.
These values give an E_F ~40-60% and facilitate checking this script with a Monte Carlo script.

Note: Miley refers to a baby with gestational age of 40 weeks.
