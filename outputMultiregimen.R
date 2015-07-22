#This script depends on outputPKPD.R, whose comments apply also to this script.
#Input is a list of many dosing regimens, not just one regimen as for outputPKPD.R.
#Output is a data frame, for all regimens input, passed to correlation.R. 
rm(list=ls())
source("outputPKPD.R")
GA = 40
outputMultiregimen = function(filename, output, weight){
  input = read.csv(filename) #load all the data
  dose_s = as.character(unlist(input)) #get a vector of doses
  dose_c = paste('c(',dose_s,')',sep='') #make all the dose to be vectors in string
  dose_vec = sapply(dose_c,function(x) parse(text=x)) #evaluate the vectors in string
  WT <<- weight
  #made WT global
  ## combined two functions: conccurve and countfunc -- this function will make it easier to create a data table##
  cure = function(dose,interval,GA){
    GAV1= -0.0114; GAmed= 28.9 #constant parameters; these and other parameters are explained in detail in Radivoyevitch
    #et al. and in forExploratoryCalculations.R
    out = conccurve(dose,interval,GA) #getting concentration from all doses
    conc = out[1:nrow(out),4][-1] #nrow=1+sum(interval)/step; 4=concentration in central PK compartment vc
    count = countfunc(conc) #function to get the count output
    nonzero<-count[ ,"m"]/(1+count[,"corr"]) #stochastic estimate for probability 1 lineage has not gone extinct.
    factor<-.406*WT*4.83e8*(1+GAV1*(GA-GAmed)) #intial average number of bacteria per neonate
    cure<-exp(-factor*nonzero) #probability every lineage has gone extinct; zero Poisson class.
    return(list(cure, conc))
  }
  
  table = data.frame(Doses =character(0), E40=numeric(0), c0.1=numeric(0),c11.99=numeric(0),#column names
                     c12.01=numeric(0), c23.99=numeric(0),c24.1=numeric(0),c35.99=numeric(0),c36.1=numeric(0),
                     c47.99=numeric(0),c48.1=numeric(0),c59.99=numeric(0),c60.1=numeric(0),
                     c71.99=numeric(0), c72.1=numeric(0),stringsAsFactors=F)
  
  ################### generating data from each set of doses ######################
  for(i in 1:length(dose_vec)){#loop thru the input dosing regimens
    print(i) #just print to check R is running when we create a table
    dose <<- array(eval(dose_vec[i])) #convert list of dose to array
    interval=c(0,rep(12,6),4) #this is used to determine what intervals correspond to the dose
    if(length(dose)==4) interval= c(0,rep(24,3),4) 
    if(length(dose)==3) interval= c(0,rep(36,2),4)
    times <<- cumsum(interval[-length(interval)])
    time_step <<- seq(0,sum(interval),by=step)
    cure_indices =c(dose_s[i]) #create a vector containg infor of each row. The first column is doses
    
    cure_vec = array(unlist(cure(dose,interval,GA)[1])) #probability of cure=E_F
    conc = array(unlist(cure(dose,interval,GA)[2])) #concentration
    cure_prob = cure_vec[7500]
    if(length(cure_vec) < 7500) cure_prob = cure_vec[length(cure_vec)]
    cure_indices = c(cure_indices,round(cure_prob,3)
                     ,round(conc[c(8,1199,1208, 2399,2408,3599,3608,4799,4808,5999,6008,7199,7208)],3))
    #adding other info we want to the row
    table[nrow(table)+1,] = cure_indices #add the row to the table
  }
  write.csv(table, output)
  return(table)
}

