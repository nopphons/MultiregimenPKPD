#This file contains three main functions that we used in our calulation: concentration, conccurve, and countfunc.
#concentration is used to define a PK response to one dose, assuming 3 compartments: central, 1, and 2. 
#conccurve is used to help superimpose concentration from different doses in one regimen (usually just one cycle)
#by looping through intervals. Finally, countfunc is used to predict bacteria number as a function of time. 
#For more details, copyright & availability information, etc.,see the introductory comments in forExploratoryCalculations.R

library (deSolve) 
step<-.01; SPH<-1/step#See comments in forExploratoryCalculations.R
concentration<-function(t,state,parameter){#t=time; "state" is drug concentration (c,p1,p2) as function of t.
  with(as.list(c(state,parameter)),{#now define some auxiliary fctns, then number (3) of ODE & of unknown functions.
    #For more details on the following, see comments in forExploratoryCalculations.R
    CL=thetaCL*(WT^0.75)*(1+thetaGACL*(GA-GAmed))*(1+((t/24)^thetaPNACL)) #allometric clearance; theta indicates a covariate.
    VC = thetaVC*WT*(1+GAVc*(GA-GAmed)) #volume VC of central compartment depends on GA
    concs=with(as.list(parameters),dose/(thetaVC*(1+GAVc*(GA-GAmed))))
    val = 0
    for (ii in 1:length(dose)){    
      val = val + concs[ii]*(times[ii]< t & t < times[ii]+delta)
    }
    val = val/delta
    dp1 <- Q1*c/VC - Q1*p1/V1 #dp1,dp2,dc= time derivatives; Q1,Q2=wash rates; V_1,2=volume of peripheral cmptmt. 1,2. 
    dp2 <- Q2*c/VC - Q2*p2/V2 #time rate of change of concentration in peripheral compartment 2.
    dc <- - Q1*c/VC + Q1*p1/V1 - Q2*c/VC + Q2*p2/V2 -CL*c/VC + val #rate of change of concentration in central cpt.
    list(c(dp1,dp2,dc))#output of the function concentration; defines the 3-compartment model
  })
}
#for more details on the following, see comments in forExploratoryCalculations.R
parameters <-c(thetaVC = 0.406,V1=0.305, V2=4.55, thetaCL =0.01, Q1 = 0.0930 , Q2 = 0.0155,
               GAmed= 28.9,thetaGACL = 0.0871, 
               GAVc= -0.0114, thetaPNACL=0.517, Kgrowth= 2,Kdeath=0.179, Bmax=8.26e8, BP =2.09e6, gamma = 1,
               Emax0= 51,EC500= 9.85, AR50= 0.113,Kon=0.0426,Koff=0.0139, doses = list(), times = list(), delta = 0.08)

conccurve = function(dose, interval,GA){
  out2<-ode(y=c(p1=0,p2=0,c=0),times=time_step,func=concentration,parms=parameters)#events = list(data = eventdat))
  return(out2)
}

countfunc <- function(conc){
  concfunc <-splinefun(conc) # temporarily make concentration into a function
  #for more details on the following, see comments in forExploratoryCalculations.R
  state <-c(s=4.83e5, r=0,m=1,corr=0, ARoff = 1, ARon = 0) #initial (t=0) values for 6 ODE in 6 unknowns
  withdrug<-function(t,state,parameters){ #model for number of CFU taking into account drug effects
    with(as.list(c(state,parameters)),{#following lines give the derivatives in the ODE and some auxiliary quantitites 
      C = concfunc(SPH*t) #get concentration from spline
      EC50 =  EC500
      dARoff <- Koff*ARon-Kon*ARoff*C #the ODE for no adaptive resistance; d=time derivative
      dARon <- Kon*ARoff*C - Koff*ARon 
      Emax <- Emax0*(1-ARon/(ARon+AR50)) #Emax model
      DRUG <- (Emax*C^gamma)/(EC50^gamma+C^gamma) #The drug effect on bacterial killing 
      ksr<-(Kgrowth-Kdeath)*(s+r)/Bmax 
      if(s<BP) ksr = 0 #ksr with cutoff; for calculating eradication probability E(t) gives same results as no cutoff
      ds<-Kgrowth*s-(Kdeath +DRUG)*s-ksr*s #deterministic equation.
      dr<-ksr*s-Kdeath*r # r is # quiescent bacteria, not susceptible to the drug even for AR off. Set=0 when calculating E(t) 
      dm<-(Kgrowth-Kdeath-DRUG)*m #m*4.83e5 approximates s if ksr is small;
      dcorr<-Kgrowth+(Kgrowth-Kdeath-DRUG)*corr #corr is a stochastic correction used in stochastic models.
      list(c(ds,dr,dm,dcorr,dARoff,dARon)) #output of the function withdrug
    })
  }
  count<-ode(y=state, times=time_step,func=withdrug,parms=parameters) #output from the withdrug function showing 
  #the number of bacteria and of other relevant quantities at time t
  return(count)
}
#count is passed to outputMultiregimin.R