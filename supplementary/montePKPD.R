#script compares analytic results with Monte Carlo results, as a check.Need to run mildlyInfectedBaby.R first
#Following lines are used to get drug rate causing cells to die at each time. 
conc= out[1:nrow(out),4] #concentration vector
ARon = count[,'ARon'] #ARon vector
Emax = with(as.list(parameters),Emax0*(1-ARon/(ARon+AR50))) #Emax vector
conc = conc[1:length(ARon)] #concentration vector with the same length as ARon
DRUG <- with(as.list(parameters),(Emax*conc^gamma)/(EC500^gamma+conc^gamma)) #drug rate vector
DRUG = DRUG[2:7601] #the first component is equal to zero

Kgrowth= 2
Kdeath=0.179

### starting and ending times of stochastic Monte Carlo estimate of the number of cells
cells = factor*dtrmn #number of cells vector
t0 = 1 #the first point where the average number of cells is between 9 and 10.5
t1 = max(which(cells <50)) #the point where the number of cells goes back up to 50
 
Ntrvl<-1000 #number of time intervals per hour; 100 might work
intrvl=1.0/Ntrvl #time in each interval
M=3# number of sample paths (Monte-Carlo runs). Larger M gives higher accuracy
time = t1 #time at the ending point
y=array(0,dim=c(M,time))#used to store the value for each path at each time

#In each monte carlo path
for (i in 1:M){if(abs(100*i/M -round(100*i/M)) < 0.5/M) print(i)
  j<-50 #the number of cells at the starting point
  for (t in t0:t1){ #for each time t
    if (j==0) { #if all cells died
      y[M,t:time] = 0 #the rest will be 0 too
      break()#once eradicated no ressurection occurs
    }
    for(k in 1:(Ntrvl/100)){ #divide into intervals here 
    #Monte Carlo, see if each cell is born or died
      x<-runif(1)
      if(x >intrvl*(Kgrowth+Kdeath+DRUG[t])*j) j<-j else{
        if(x >intrvl*Kgrowth*j) j<-j-1 else j<-j+1}
      if (j==0) { #if all cells died
        y[M,t:time] = 0 #the rest will be 0 too
        break()#once eradicated no ressurection occurs
      }
    }
    y[i,t] = j #store the value in the path
      }
  
}

#plot the graph of the number of cells by calculating the mean of all paths
  
plot(colMeans(y),type="l",ylim = c(0,50))
#the graph of the number of cells (deterministically)
lines(0:(t1-t0),factor*dtrmn[t0:t1], col='red')

####percent gone extinct between t0 and t1#######
1-length(which(y[,t1]>0))/M

plot(factor*dtrmn, type="l",xlim=c(0,800), ylim=c(0,50), col="blue")
#in above, adjust xlim and ylim to show the most interesting region, where 
#the average is at its minimum
lines(y[1,], cex=0.4,col="red")
lines(y[2,], cex=.4)#keep going till find two different sample paths,
#with one going extinct and the other above average. 

#Zoom in (specify the range below)
range = 200:320
plot(range, y[1,range], type="s", ylim = c(min(y[1,range]),max(y[1,range])), col = 'red')


