setwd("~/OneDrive - Oxford University Clinical Research Unit/Material_courses/MIDSEA summer school/simulation materials/")


library(odin) 
source('code/sim_code_01a-helper.r')


# Models ------------------------------------------------------------------
# Deterministic, for comparison:
sir_ode <- odin({
  ## Rates of change
  deriv(S) <- -beta * S * I 
  deriv(I) <-  beta * S * I - gamma * I 
  deriv(R) <- gamma * I 
  
  ## Total population
  N <- S + I + R
  
  ## Initial conditions
  initial(S) <- N_initial - I_initial - R_initial
  initial(I) <- I_initial
  initial(R) <- R_initial
  
  #N_initial <- user(99400000) # population size of Viet Nam
  #N_initial <- user(457400) # population size of Quy Nhon
  N_initial <- user(50000) # population size of small town
  I_initial <- user(1) # Infectious
  R_initial <- user(0)# Immune
  
  ## user defined parameters (if blank they are specified outside the function)
  R0 <- user()
  gamma <- user() 

  ## derived parameters if any
  beta <- R0*gamma/N
})



gillespie = function(maxt,state,p)
{
  #Q01: What is the difference between Sv and state$S?
  populationsize = state['N']
  # to help memory, precreate the vector at maximum possible length
  
  Sv = Iv = Rv = tv = rep(0, populationsize*2) 
  
  Iv[1] = state['I']
  Rv[1] = state['R']
  Sv[1] = populationsize-state['I']-state['R']
  tv[1] = 0
  
  
  
  beta <- p$R0*p$gamma/populationsize
  
  nevents = 1 #Q02: What does this do?
  continuerunning = TRUE
  printing_waypoint = 0
  while(continuerunning) #Q03: What stops this loop?
  {
    infectionrate = beta*Sv[nevents]*Iv[nevents]
    removalrate = p$gamma*Iv[nevents]
    delta = rexp(1,infectionrate+removalrate) #Q04: What does this do?
    if((tv[nevents]+delta)>maxt)continuerunning=FALSE
    if((tv[nevents]+delta)<=maxt)
    {
      #Q05: what is going on in the next 3 lines?
      u = runif(1)
      eventtype = 'removal'
      if(u<(infectionrate/(infectionrate+removalrate)))eventtype='infection'
      nevents = nevents + 1
      tv[nevents] = tv[nevents-1]+delta
      if(eventtype=='infection')
      {
        Sv[nevents] = Sv[nevents-1]-1
        Iv[nevents] = Iv[nevents-1]+1
        Rv[nevents] = Rv[nevents-1]+0
      }
      if(eventtype=='removal')
      {
        Sv[nevents] = Sv[nevents-1]+0
        Iv[nevents] = Iv[nevents-1]-1
        Rv[nevents] = Rv[nevents-1]+1
      }
      if(Iv[nevents]==0)continuerunning=FALSE
      if(tv[nevents]>printing_waypoint)
      {
        cat('t ',printing_waypoint,' infections ',Iv[nevents],'\n',sep='')
        printing_waypoint=printing_waypoint+1
        }
        
    }
  }
  
  tv = tv[1:nevents]
  Sv = Sv[1:nevents]
  Iv = Iv[1:nevents]
  Rv = Rv[1:nevents]
  
  tv=c(tv,maxt) # fill up to end time
  Sv=c(Sv,Sv[nevents])
  Iv=c(Iv,Iv[nevents])
  Rv=c(Rv,Rv[nevents])
  output = data.frame(t=tv,S=Sv,I=Iv,R=Rv) 
  return(output) #Q06: What is getting output by the function?
}


# Common objects ----------------------------------------------------------
maxt <- 365
times <- 1:maxt
parameters <- list(R0=2,gamma=0.2,N=1000)
state <- c(I = 1,R = 0,N=1000)



# Run deterministic model -------------------------------------------------

sir_model <- sir_ode$new(R0 = parameters$R0, gamma = parameters$gamma,N_initial=parameters$N)
time_horizon <- seq(from = 0,to = maxt,length.out=10001)
modelrun <- sir_model$run(time_horizon)
# Convert to data frame
df_detmodelrun = as.data.frame(modelrun)

plotsirD(df_detmodelrun,'example01a.png')



# Run stochastic model ----------------------------------------------------
df_stochmodelrun = gillespie(maxt,state,parameters)
plotsirS(df_stochmodelrun,'example01b.png')





rm(df_detmodelrun,df_stochmodelrun,modelrun,parameters,sir_model,maxt)
rm(sir_ode,state,time_horizon,times,gillespie,plotsirD,plotsirS)
