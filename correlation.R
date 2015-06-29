#This script performs the main calculations for the submitted paper Siranart et al., qv. 
#The script is used first to parameterize a Hill (sometimes called log-logistic or Fisk) model for
#approximating E. Coli extinction probabilities. We use simulated annealing and information about
#the prototypical neonate with gestation age GA=40 wk to find optimal Hill model parameters. 
#We output a table ("miley_output_56_doses.csv" and a figure; both compare
#the Hill model with calculated extinction probabilities for the prototypical neonate. The table also contains
#gentamicin concentration information relevant to toxicity and to more conventional measures of efficacy.
#Then we compare the parametrized Hill model with a more complicated model when both are applied
#to 100 neonates of the same GA but different birthweights WT, drawn from a standard neonate WT distribution.
#The output is a histogram showing the distribution of 100 correlation coefficients.
#For more details, copyright & availability information, etc.,see the introductory comments for forExploratoryCalculations.R
#The script calls everything it needs: outputMultiregimen.R and .csv inputfiles

library(GenSA)
library(visreg)
library(VGAM)

#Preprocess data
rm(list=ls())
set.seed(22222)# comment out this line to get different results for each run
source("outputMultiregimen.R")#loads many parameters; refer to forExploratoryCalculations.R for details on parameters.

###choose GA and number of neonates ######
GA = 40
miley = 3.5; mu = 3486; sigm=434#weight of prototype neonate (kg) and weight distribution parameters (g) for given GA
#"miley" is a hypothetical nicknames for the (hypothetical) prototypical GA 40 wk female neonate.
number =2# number=a few neonates when testing the script; but usually=100 in actual calculations.
weights=(rnorm(number,mu,sigm))*0.001#defines a virtual population of number neonates, differing by birth weight and
#therefore by many parameters that depend on birthweight. 

####################Preprocessing#########################
dose_input = 'miley_input_56_doses.csv'#loads the chosen dose cycles
#table = outputMultiregimen(dose_input, 'initializingTable.csv', miley)
table = read.csv('initializingTable.csv')# used as initialization; will be overwritten by script below
l = nrow(table)
E40 = as.double(table[1:l, "E40"])#E_F for the prototypical GA 40 wk female neonate and all 50 dose cycles
rank = order(E40)# order the 50 rows by E_F
dose_s = as.character(unlist(table[1:l, "Doses"])) #get a vector of doses.
dose_c = paste('c(',dose_s,')',sep='') #make all the doses to be vectors in string
dose_vec = sapply(dose_c,function(x) parse(text=x)) #change from string to expression
doses = sapply(dose_vec, eval) #evaluate the expression
dose1 = doses[1,]; dose2 = doses[2,]; dose3 = doses[3,]; dose4 = doses[4,] #get the doses from the dosing cycle
E40 = as.double(table[1:l, "E40"]) #get E40

#Find B coefficients. These are obtained using E_F values for the prototypical neonate as (virtual) data; then used for other neonates.

#This function is used to find the correlation of the predicted value from the model
#given the coefficients and the E_F values for the prototypical baby
rrho <- function(coeff) {
  a= coeff[1];b=coeff[2];c=coeff[3];d=coeff[4]
  x = a*dose1+b*dose2+c*dose3
  val = pfisk(x, 1, d)
  1-cor(val,E40)
}
#We used simulated annealing to optimize parameters of the Hill function
global.min = 0; tol = 0.0001; lower = rep(0,4); upper = c(1,1,1,20)
## during simulated annealing, do that and simplify the following. However the following works
out <- GenSA(par = c(0.5,0.5,0.5,10), lower = lower, upper = upper, fn = rrho,
             control=list(threshold.stop=global.min+tol,verbose=TRUE))
result = out[c("value","par","counts")] #result from simulated annealing
boncoeff = result$par #optimal coefficients, but renormalization is needed due to parameter redundancy 
#renormalized result for GA=40 on 5/25/2015: boncoeff[1]/boncoeff[3]=0.2606252, boncoeff[2]/boncoeff[3]=0.005040724, boncoeff[4]=9.871
i_sim = boncoeff[1]*dose1+ boncoeff[2]*dose2+boncoeff[3]*dose3 #B for the Hill (i.e. Fisk) function 1/[1+(1/B)^b], where the shape parameter b=boncoeff[4]
E40_p = pfisk(i_sim,1,boncoeff[4])  #the predicted probability E40=E_F

#This is the function that is used to predict the probability given the value B for the hill.
E40_func <- function(B){
  pfisk(B,1,boncoeff[4])
}

#for the protoypical neonate, plot the virtual data points and the optimized Hill function
rank = order(i_sim) #get the order of i_sim to generate the plot
plot(i_sim[rank],E40[rank], type = 'p', xlab = 'B=b1*d1+b2*d2+b3*d3', ylab = 'E40')# to properly renormalize, x-axis labels must be divided by 80.5
legend("bottomright", legend = c('E40','E40_predicted'), col = c('black','red'), lty=1)
curve(E40_func, min(i_sim),max(i_sim),add=TRUE, col = 'red') #plot the curve of the predicted prob E_F
correlation = cor(E40_p, E40)#Pearson correlation applied to non-linear curve.
table[,'B=b1d1+b2d2+b3d3'] = i_sim #add B to the table. Comment out following line to retain 'miley_output_50_doses.csv'
write.csv(table, 'miley_output_56_doses.csv')

#########################SIMULATION###########################
#We consider neonates with the same GA but different weights. For each such neonate we compare the Hill function,
#with the same parameters as above to the E_F values calculated with the complex, biologically motivated model.
#Using 50 comparisons we comput Pearson's r. Repeating for the other 99 neonates gives a distribution, which
#we plot as a histogram.
rho=c(); #Empty lists to store the correlation of each neonate
for(j in 1: length(weights)){
  print(j)
  sim_table = outputMultiregimen(dose_input, 'output', weights[j]) #a table for each neonate
  sim_table = sim_table[rank,] 
  E40_sim = sim_table['E40'] #get E40 of each cycle
  E40_sim = as.double(E40_sim[,1]) 
  rho = c(rho, cor(i_sim ,E40_sim))
}
distribution = array(rho) 
hist(rho, breaks = 40, main="Histogram of Correlation") #plot the histogram of correlations

### Kolmogorov-Smirnov Test (may be biased; see supplement to Siranart et. al. for more details) ###
set.seed(1234567)
ks.pval.exact = ks.test(jitter(E40), jitter(E40_p))$p.value
### The "jitter" function was used to add a little bit of random noise to E40 and E40_p.
### Since there are duplicate values in E40 and E40_p, the use of the "jitter" function is
### used to remove duplicate values, which would otherwise prevent the computation of exact
### p-value for the Kolmogorov-Smirnov Test.

ks.pval.approx = ks.test(E40, E40_p)$p.value
### Alternatively, an approximate p-value can be found for the Kolmogorov-Smirnov Test even if
### E40 and E40_p contain duplicate values.