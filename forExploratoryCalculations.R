#Written & commented by N. Siranart 04/14; tested, modified and further commented by R.Sachs, Y.Chen & T.Radivoyevitch.
#Distributed under the GNU GPLv3; Contact Sachs. sachs@math.berkeley.edu., for .txt or.R versions.
#This script uses methods similar to (Nielsen et al 2009 Clin Pharmacokinet 48:253) for neonatal gentamicin dosing and uses
#corresponding PK/PD calculations in Mohamed et al. 2012, Antimicrob Agents Chemother 56:179. It adds stochastic calculations.
#In the stochastic calculation of eradication probability E(t) the simplified model with ksr=0=r was always used.
#The script was heavily used in Radivoyevitch et al. (2015) AAPS Journal 17 (2):447-456. doi:10.1208/s12248-014-9715-3
#It is now just a self-contained seperate unit used for preliminary calculations, to output figures, and as a tutorial.
rm(list=ls())
library (deSolve) #We will solve systems of non-linear first order ODE (ordinary differential equations)

#The next few commands are adjustable by the user; they set the dosing cycle and the neonate gestational age GA.
#For details, including meaning of parameters, see Nielsen et al. For self consistency and consistency with other literature,
#the notation of Nielsen et al. was modified as follows: (1) subscripts written on the main line. (2)(greek lower case theta)
#goes to theta. An example which uses both these changes is (greek lower case theta)_GACL->thetaGACL;
#(3) CL->thetaCL; (4) V1->thetaVC (with C standing for the central compartment); (5) TVV1->VC, Q2->Q1, V2->V1, Q3->Q2, V3->V2.
GA =40 #gestational age (weeks); possible values are 25, 29, 34, 40; determines neonate parameters such as weight

#dose<-rep(2.2,4); interval<-c(0,rep(12,3),4) #
#drug doseS (mg/kgm); intervalS (hours) between doses. For example dose=c(5,4,4) and interval=c(0,24,24,4)
#would mean a bolus of 5 at t=0, then 2 fractions of 4 at 24h intervals, and then a 4 h dose-free recovery period.
dose<-c(4,4); interval<-c(0,24,52) 
#dose<-c(4,4,0,0); interval<-c(0,rep(24,3),4)#should give the same answer
#dose<-rep(4,6); interval<-c(0,rep(36,5),60)#Fig. From Nielsen et al.
time = cumsum(interval[-length(interval)])
step<-.01; SPH<-1/step #ODE output times: step=time (h) for one output step; SPH=steps per hour.
#Generally in this script, time can refer to hours, days, weeks, or number of steps; one has to keep track.
times<-seq(0,sum(interval),by=step) #Time ranges for the ODE output, for example 76 hours in 7601 steps
#************** end of user adjusted parameters *************

weight=c('25'=0.778,'29'= 1.200,'34'=2.710,'40'=3.500) #WT=weight (kg) of neonate at post natal age (PNA) 0 days; 

#Next define a PK response to one dose, assuming 3 compartments: central, 1, and 2
concentration<-function(t,state,parameter){ #t=time; state is drug concentration (c,p1,p2) as function of t
  with(as.list(c(state,parameter)),{#now define some auxiliary fctns, then number (3) of ODE & of unknown functions.
    WT = unname(weight[toString(GA)])#next uses parameters from Nielsen and Mohamed
    CL=thetaCL*(WT^0.75)*(1+thetaGACL*(GA-GAmed))*(1+((t/24)^thetaPNACL)) #allometric clearance; theta indicates a covariate.
    VC = thetaVC*WT*(1+GAVc*(GA-GAmed)) #volume VC of central compartment depends on GA; see Nielsen et al
    dose = parameter['doses']$doses; time = parameter['times']$times
    concs=with(as.list(parameters),dose/(thetaVC*(1+GAVc*(GA-GAmed))))
     val = 0
     for (ii in 1:length(dose)){    
       val = val + concs[ii]*(time[ii]< t & t < time[ii]+delta)
     }
    val = val/delta
    dp1 <- Q1*c/VC - Q1*p1/V1 #dp1, dp2, dc= time derivatives; Q1,Q2=wash rates; V1,2=volume of peripheral cmpt. 1,2. 
    dp2 <- Q2*c/VC - Q2*p2/V2 #time rate of change of concentration in peripheral compartment 2.
    dc <- - Q1*c/VC + Q1*p1/V1 - Q2*c/VC + Q2*p2/V2 -CL*c/VC + val #rate of change of concentration in central cpt.
    list(c(dp1,dp2,dc)) #output of the function concentration; defines the 3-compartment model
  })
}
# Following numbers are from Table III "combined dataset" in Nielsen et al. and from Mohamed et. al. Central values 
# are used, not "population" PK/PD. Mohamed Koff=.0139, may be unrealistically small. 
parameters <-c(thetaVC = 0.406,V1=0.305, V2=4.55, thetaCL =0.01, Q1 = 0.0930 , Q2 = 0.0155,
               GAmed= 28.9,thetaGACL = 0.0871, 
               GAVc= -0.0114, thetaPNACL=0.517, Kgrowth= 2,Kdeath=0.179, Bmax=8.26e8, BP =2.09e6, gamma = 1,
              Emax0= 51,EC500= 9.85, AR50= 0.113,Kon=0.0426,Koff=0.0139, doses = list(dose), times = list(time), delta = 0.08)
              #Emax0= 51,EC500= 9.85, AR50= 0.113,Kon=0.0426,Koff=0.0278, doses = list(dose), times = list(time), delta = 0.08)

#PK model for getting concentration from all doses. The above concentration function is only for one dose. 
#This new function will help superimpose concentrations from different doses by looping through intervals.
conccurve = function(dose, interval,GA){
  out2<-ode(y=c(p1=0,p2=0,c=0),times=times,func=concentration,parms=parameters)#events = list(data = eventdat))
  return(out2)
}

out = conccurve(dose,interval,GA) #output of above
conc = out[1:nrow(out),4] #nrow=1+sum(interval)/step; 1=time (h), 2,3, & 4 = conc. in v1, v2 & vc, Here 4 used.
conc = conc[-1] #remove the first component (which is 0)

#Now predict susceptible bacteria number as a function of time.
#Some of the calculations are similar to those used in Mohamed et al. The main method will be
# solving a system of non-linear first order ODE, usable for both deterministic and stochastic results.
#conc is now a vector with, for example, 7601 components

countfunc <- function(conc){
  conbacterianc <-splinefun(conc) # temporarily make concentration into a function
  state <-c(s=4.83e5, r=0,m=1,corr=0, ARoff = 1, ARon = 0) #initial (t=0) values for 6 ODE in 6 unknowns
  #s=bacteria per mL; r refers to quiescent bacteria. Initially r=0 (preceding line). For us, r usually remains 0 throughout.
  #corr is an auxiliary function required for stochastic corrections to deterministic estimates (see Radivoyevitch et al.)
  #AR is adaptive resistance; ARon, ARoff are fractions with and without adaptive resistance
  withdrug<-function(t,state,parameters){#model for number of bacteria taking into account drug effects
    with(as.list(c(state,parameters)),{#following lines give the derivatives in the ODE and some auxiliary quantitites 
      C = conbacterianc(SPH*t) #get concentration from spline
      EC50 =  EC500 #not EC500*(1 -exp(-10*Cr*t)); the optional extra factor of the references is never used in our calculations
      dARoff <- Koff*ARon-Kon*ARoff*C #the ODE for no adaptive resistance; d=time derivative
      dARon <- Kon*ARoff*C - Koff*ARon #time derivative of the fraction of bacteria that have acquired adaptive resistance
      Emax <- Emax0*(1-ARon/(ARon+AR50)) #Emax  model (uses Hill i.e.log-logistic function)
      DRUG <- (Emax*C^gamma)/(EC50^gamma+C^gamma) #The drug effect on bacterial killing, another Hill function. For us gamma=1 
      ksr<-(Kgrowth-Kdeath)*(s+r)/Bmax #ksr is transfer rate constant from cycling to quiescent bacteria; in our
      #calculation the number of quiescent bacteria remains 0 until bacteria numbers become so large stochastics no longer matter
      if(s<BP) ksr = 0 #ksr with cutoff; for calculating eradication probability E(t) gives same results as no cutoff
      ds<-Kgrowth*s-(Kdeath +DRUG)*s-ksr*s #deterministic growth and death equation for cycling bacteria
      dr<-ksr*s-Kdeath*r # r is # quiescent bacteria, not susceptible to the drug even for AR off. Set=0 when calculating E(t) 
      dm<-(Kgrowth-Kdeath-DRUG)*m #m*4.83e5 approximates s if ksr is small;
      dcorr<-Kgrowth+(Kgrowth-Kdeath-DRUG)*corr #corr is a stochastic correction function.
      list(c(ds,dr,dm,dcorr,dARoff,dARon)) #output of the function withdrug
    })
  }
  count<-ode(y=state, times=times,func=withdrug,parms=parameters) #output from the withdrug function showing 
  #the number of bacteria and of other relevant quantities at time t
  return(count)
}

count = countfunc(conc)
dtrmn<-count[ ,'s']/4.83e5 # Deterministic s/s(0), where s is number in central compartment
on<-count[,"ARon"]# the fraction of bacteria that have adaptive resistance at time t

####### Collect results needed to plot figures ######################################

thetaCL = 0.01
WT = unname(weight[toString(GA)]) #roundabout way to get weight as a numeric
GAmed= 28.9 #GA in weeks. Median for nmany observations, see Nielsen et al.
thetaGACL = 0.0871 #helps calculate correction of clearance for GA different from median.
thetaPNACL=0.517 #another covariate
t= seq(0,sum(interval),by=step) #for example t= 0, .01, .02, ..., 75.99, 76 in some cases.
a = round(t/24)*24 # using a instead of t in the following line gives step functions as in Nielsen et al. smooth is preferred here
CL=thetaCL*(WT^0.75)*(1+thetaGACL*(GA-GAmed))*(1+((t/24)^thetaPNACL))# gives smooth function; more plausible mechanistically
GAVc= -0.0114
factor<-.406*WT*4.83e8*(1+GAVc*(GA-GAmed))#intial average number of bacteria per neonate NOT per mL;vc = thvc*WT*(1+GAVc*(GA-GAmed))
print(min(factor*dtrmn))#deterministic approximation to the the minimum number of neonates 
nonzero<-count[ ,"m"]/(1+count[,"corr"])#stochastic estimate that a lineage has not gone extinct (ksr=0 approximation)
finalnz<-nonzero[(sum(interval)-1)*SPH]
cure<-exp(-factor*nonzero) #probability every lineage has gone extinct, E(t) in text.
print(cure[(sum(interval)-1)*SPH])#approximates long time limit E_F of E(t)

########## To make various figures, set plot=TRUE
# # Check that figures are similar to those in Nielsen et al. and Mohamed et al. Add stochastic extinction plot
plot = TRUE
if(plot){
  par(fin=c(4,5),ann=FALSE, mai=rep(.1,4),xaxp=c(0,sum(interval[-length(interval)]),1))# parameters specifying layout
  plot(out[ ,"time"],out[ ,"c"],ylim=c(0,17),xlim=c(-1,sum(interval)),yaxp=c(0,16,4),type='l')#concentration as function of t
  abline(h=c(0,2,8,10,12,14,16),lty=2)
  ex<-cumsum(interval)+1
  points(ex,conc[SPH*ex],pch=8) #points one hour after each dose fraction shown as asterisks
  legend(20,9,legend=c(dose,GA), cex=.5)
  #title(paste0("Doses=",paste(dose,collapse=","),"; GA=",GA))
  
  plot(times,factor*dtrmn,log='y',ylim=c(.1,1e15),xlim=c(step, sum(interval)), type='l')#plot deterministic number bacteria vs. t
  abline(h=c(1,100,1e8), lty=2) #horizontal lines as in Mohamed Figs 4, 5, 6 bottom row.
  legend(x=5,y=1e11, c(dose,GA),cex=.65)
  lines(times,factor*count[ ,"m"], lty=3) #deterministic ksr=0 approximation
  
  plot(on, ylim=c(0,1),type="l") #fraction of bacteria that have acquired adaptive resistance
  lines(CL,type='l', lty=2) #Clearance for central compartment increases as neonate ages
  legend(20*SPH,.6,c(dose,interval[2], GA),cex=.75)
  
  plot(cure, type="l",ylim=c(0,1)) #extinction probability E(t)
  abline(h=cure[length(cure)])
  #round(conc[c(10,1198,1210,2398,2410,3598,3610,4798,4810,5998,6010,7198,7210)],2)#used for Table III of Radivoyevitch paper.
}