rm(list=ls())
set.seed(22222)# comment out this line to get different results for each run
source("outputMultiregimen.R")
dose_input = 'miley_input_56_doses.csv'
input = read.csv(dose_input)
dose_s = as.character(unlist(input))
dose_c = paste('c(',dose_s,')',sep='')
dose_vec = sapply(dose_c,function(x) parse(text=x))
doses = sapply(dose_vec, eval)
dose1 = as.vector(doses[1,]);

table = data.frame(doses =character(0),stringsAsFactors=F)

for(i in 1:length(dose1)){
  for(j in 1:150){
    d1 = dose1[i]
    d2 = 0; d3 = 0; d4 = 0
    if(runif(1)>0.2) d2 = round(runif(1,3,5.5),1)
    if(runif(1)>0.2) d3 = round(runif(1,3,5.5),1)
    if(runif(1)>0.2) d4 = round(runif(1,3,5.5),1)
    dose_indices = paste(toString(d1),toString(d2),toString(d3),toString(d4),sep=",")
    table[nrow(table)+1,] = dose_indices
  } 
}
write.csv(table, 'input_dose.csv', row.names=FALSE)
new_input = 'input_dose.csv'
GA = 40; miley = 3.5
table = outputMultiregimen(new_input, 'output_dose.csv', miley)
l = nrow(table)
E40 = as.double(table[1:l, "E40"])
dose_s = as.character(unlist(table[1:l, "Doses"])) #get a vector of doses.
dose_c = paste('c(',dose_s,')',sep='') #make all the doses to be vectors in string
dose_vec = sapply(dose_c,function(x) parse(text=x)) #change from string to expression
doses = sapply(dose_vec, eval) #evaluate the expression
dose1 = doses[1,]; dose2 = doses[2,]; dose3 = doses[3,]; dose4 = doses[4,] 

rank = order(E40/dose1)
sort_table = table[rank,]
write.csv(sort_table, 'output_dose.csv')






