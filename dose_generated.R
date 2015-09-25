rm(list=ls())
set.seed(5888)# comment out this line to get different results for each run
source("outputMultiregimen.R") #loads many parameters; refer to forExploratoryCalculations.R for details on parameters.
dose_input = 'csv_files/miley_input_55_doses.csv' #loads the chosen dose cycles; other input files of the same structure can be used.
input = read.csv(dose_input) #load all the data
dose_s = as.character(unlist(input)) #get a vector of doses
dose_c = paste('c(',dose_s,')',sep='') #make all the dose to be vectors in string
dose_vec = sapply(dose_c,function(x) parse(text=x)) #evaluate the vectors in string
doses = sapply(dose_vec, eval) #evaluate the expression
dose1 = as.vector(doses[1,]) #get the doses from the dosing cycle
table = data.frame(doses =character(0),stringsAsFactors=F) #initialize an empty table to add info for each dose
sort_table = data.frame(doses =character(0),stringsAsFactors=F) #initialize an empty table to add info for each front-loaded dose
for(i in 1:length(dose1)){ #for each dose in dose_input
  for(j in 1:15){  #we simulate 15 times (can be changed)
    d1 = dose1[i] #set dose1 equal to the first dose of dose_input
    d2 = 0; d3 = 0; #d2 and d3 are chosen from distribution explained in the paper
    if(runif(1)>0.2) d2 = round(runif(1,3,5.5),1) 
    if(runif(1)>0.2) d3 = round(runif(1,3,5.5),1)
    dose_indices = paste(toString(d1),toString(d2),toString(d3),sep=",") #concatenate strings
    table[nrow(table)+1,] = dose_indices #add to the table
    
    d = sort(c(d1,d2,d3), decreasing = TRUE) #sort the dose we got in the previous simulation
    d1 = d[1];d2=d[2]; d3 = d[3];  
    dose_indices = paste(toString(d1),toString(d2),toString(d3),sep=",") #concatenate strings
    sort_table[nrow(sort_table)+1,] = dose_indices #add to the sort_table
  } 
}
#write tables to csv files
write.csv(table, 'csv_files/input_825_doses.csv', row.names=FALSE)
write.csv(sort_table, 'csv_files/front_input_825_doses.csv', row.names=FALSE)




