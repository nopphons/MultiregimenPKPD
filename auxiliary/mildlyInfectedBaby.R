#This script is almost the same as forExploratoryCalculations.R, whose comments should be consulted.
#However, we consider a virtual neonate which has only about 50 bacteria at t=0 and only give 1, very small, dose.
#These values give an E_F ~40-60% and facilitate checking this script with a Monte Carlo script.
###################################################################

library (deSolve)
GA =40 #gestational age (months).
dose<-c(0.310,0,0,0); interval<-c(0,rep(24,3),4) 
time = cumsum(interval[-length(interval)])
step<-.01; SPH<-1/step
times<-seq(0,sum(interval),by=step)
#**************end of user adjusted parameters************
weight=c('25'=0.778,'29'= 1.200,'34'=2.710,'40'=3.500)
concentration<-function(t,state,parameter){ 
  with(as.list(c(state,parameter)),{
    WT = unname(weight[toString(GA)])
    CL=thetaCL*(WT^0.75)*(1+thetaGACL*(GA-GAmed))*(1+((t/24)^thetaPNACL)) 
    VC = thetaVC*WT*(1+GAVc*(GA-GAmed)) 
    dose = parameter['doses']$doses; time = parameter['times']$times
    concs=with(as.list(parameters),dose/(thetaVC*(1+GAVc*(GA-GAmed))))
    val = concs[1]*(t<delta)+concs[2]*(time[2]< t & t < time[2]+delta) + concs[3]*(time[3]< t & t < time[3]+delta)
    if(length(dose)==4) val = val + concs[4]*(time[4]< t & t < time[4]+delta)
    val = val/delta
    dp1 <- Q1*c/VC - Q1*p1/V1 
    dp2 <- Q2*c/VC - Q2*p2/V2 
    dc <- - Q1*c/VC + Q1*p1/V1 - Q2*c/VC + Q2*p2/V2 -CL*c/VC + val
    list(c(dp1,dp2,dc)) 
  })
}
parameters <-c(thetaVC = 0.406,V1=0.305, V2=4.55, thetaCL =0.01, Q1 = 0.0930 , Q2 = 0.0155,
               WT = unname(weight[toString(GA)]), GA, GAmed= 28.9,thetaGACL = 0.0871, 
               GAVc= -0.0114, thetaPNACL=0.517, Kgrowth= 2,Kdeath=0.179, Bmax=8.26e8, BP =2.09e6, gamma = 1,
               Emax0= 51,EC500= 9.85, AR50= 0.113,Kon=0.0426,Koff=0.0139, infect = 0.04029,
               doses = list(dose), times = list(time), delta = 0.08)#infect, CFU/ml, is chosen so that
               # the total number of bacteria in the central compartment is initially ~50 instead of ~10^9

conccurve = function(dose, interval,GA){
  #concs=with(as.list(parameters),dose/(thetaVC*(1+GAVc*(GA-GAmed))))
  #eventdat <- data.frame(var = "c", time = cumsum(interval[-length(interval)]), value = concs, method = "add")
  out2<-ode(y=c(p1=0,p2=0,c=0),times=times,func=concentration,parms=parameters)#events = list(data = eventdat))
  return(out2)
}

out = conccurve(dose,interval,GA) 
conc = out[1:nrow(out),4] 
conc = conc[-1] 

times<-seq(0,sum(interval),by=step) 
countfunc <- function(conc){
  concfunc <-splinefun(conc) 
  state <-c(s=0.04029, r=0,m=1,corr=0, ARoff = 1, ARon = 0) #initial values for 6 ode in 6 unknowns
  withdrug<-function(t,state,parameters){ 
    with(as.list(c(state,parameters)),{
      C = concfunc(SPH*t) 
      EC50 =  EC500
      dARoff <- Koff*ARon-Kon*ARoff*C 
      dARon <- Kon*ARoff*C - Koff*ARon 
      Emax <- Emax0*(1-ARon/(ARon+AR50)) 
      DRUG <- (Emax*C^gamma)/(EC50^gamma+C^gamma) 
      ksr<-(Kgrowth-Kdeath)*(s+r)/Bmax
      if(s<BP) ksr = 0 #ksr with cutoff; for our purposes gives same results as no cutoff
      ds<-Kgrowth*s-(Kdeath +DRUG)*s-ksr*s 
      dr<-ksr*s-Kdeath*r 
      dm<-(Kgrowth-Kdeath-DRUG)*m
      dcorr<-Kgrowth+(Kgrowth-Kdeath-DRUG)*corr #corr is a stochastic correction used in stochastic models.
      list(c(ds,dr,dm,dcorr,dARoff,dARon)) 
    })
  }
  count<-ode(y=state, times=times,func=withdrug,parms=parameters) 
  return(count)
}
count = countfunc(conc)

dtrmn<-count[ ,'s']/0.04029
on<-count[,"ARon"]

######### plot results 
par(fin=c(4,4),ann=FALSE, mai=rep(.1,4))# parameters specifying layout
thetaCL = 0.01
WT = unname(weight[toString(GA)]) 
GAmed= 28.9 
thetaVC = 0.406
thetaGACL = 0.0871 
thetaPNACL=0.517 
t= seq(0,sum(interval),by=step) 
a = round(t/24)*24 
CL=thetaCL*(WT^0.75)*(1+thetaGACL*(GA-GAmed))*(1+((t/24)^thetaPNACL))

#now make plot similar to Nielsen Fig. 3b
par(fin=c(4,5),ann=FALSE, mai=rep(.1,4),xaxp=c(0,sum(interval[-length(interval)]),1))# parameters specifying layout
plot(out[ ,"time"],out[ ,"c"],ylim=c(0,17),xlim=c(-1,sum(interval)),yaxp=c(0,16,4),type='l') # concentration as a function of time
abline(h=c(0,2,8,10,12,14,16),lty=2)
ex<-cumsum(interval)+1
points(ex,conc[SPH*ex],pch=8) # points one hour after each dose fraction shown as asterisks
legend(20,9,legend=c(dose,GA), cex=.5)
title(paste0("Doses=",paste(dose,collapse=","),"; GA=",GA))

GAmed= 28.9;GAVc= -0.0114
factor<- thetaVC*WT*0.04029*1000*(1+GAVc*(GA-GAmed))
print(min(factor*dtrmn))
#plot(times,factor*dtrmn,log='y',ylim=c(.1,1e13),xlim=c(step, sum(interval)), type='l')
plot(times,factor*dtrmn,ylim=c(.1,50),xlim=c(step, 8), type='l')
abline(h=c(1,100,1e8), lty=2) 
legend(x=5,y=1e11, c(dose,GA),cex=.65)
lines(times,factor*count[ ,"m"], lty=3) #deterministic ksr=0 approximation

nonzero<-count[ ,"m"]/(1+count[,"corr"]) 
finalnz<-nonzero[(sum(interval)-1)*SPH]
one<-count[ ,"m"]/(1+count[,"corr"])^2 
cure<-exp(-factor*nonzero) #probability every lineage has gone extinct
print(c("cured",(cure[(sum(interval)-1)*SPH])))
plot(cure, type="l", log="y", ylim=c(1e-3,1)) #extinction probability E(t)
legend(.7*sum(interval)*SPH,.05,c(dose,interval[2], GA,round(cure[(sum(interval)-1)*SPH],5)),cex=.75)
plot(cure, type="l", xlim =c(0,800)) ##non-logarithmic plot of key region
#plot(on, ylim=c(0,1),type="l")
#lines(CL,type='l', lty=2) 
#legend(20*SPH,.6,c(dose,interval[2], GA,(cure[(sum(interval)-1)*SPH])),cex=.75)
