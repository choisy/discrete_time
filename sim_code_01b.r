setwd("~/OneDrive - Oxford University Clinical Research Unit/Material_courses/MIDSEA summer school/simulation materials/")

source('code/sim_code_01a-helper.r')



# Models ------------------------------------------------------------------


discretetime = function(maxt,state,p)
{
  Sv = Iv = Rv = tv = rep(0,1+maxt)
  populationsize = state['N']
  Iv[1] = state['I']
  Rv[1] = state['R']
  Sv[1] = populationsize-state['I']-state['R']
  tv[1] = 0
  
  beta <- p$R0*p$gamma/populationsize
  for(timeindex in 1:length(tv))
  {
    infectionrate = beta*Iv[timeindex]
    removalrate = p$gamma
    infectionprob = 1-exp(-infectionrate)
    removalprob = 1-exp(-removalrate)
    n_infections = rbinom(1,Sv[timeindex],infectionprob)
    n_removals = rbinom(1,Iv[timeindex],removalprob)
    Sv[1+timeindex] = Sv[timeindex] - n_infections
    Iv[1+timeindex] = Iv[timeindex] + n_infections - n_removals
    Rv[1+timeindex] = Rv[timeindex] + n_removals
    tv[1+timeindex] = tv[timeindex] + 1
  }
  
  output = data.frame(t=tv,S=Sv,I=Iv,R=Rv) 
  return(output) 
}


# Common objects ----------------------------------------------------------
maxt <- 365
times <- 1:maxt
parameters <- list(R0=2,gamma=0.2)
state <- c(I = 1,R = 0,N=1000)  #457400




# Run stochastic model ----------------------------------------------------
df_stochmodelrun = discretetime(maxt,state,parameters)
plotsirS(df_stochmodelrun,'example01b01.png')
Nsim <- 100
onemonth <- rep(0, Nsim)
Smatrix =matrix(0,Nsim, maxt+2)
Imatrix =matrix(0,Nsim, maxt+2)
Rmatrix =matrix(0,Nsim, maxt+2)
for (simulation in 1:Nsim){
  cat("Simulation", simulation, 'of', Nsim, "\n") 
  #run
  df_stochmodelrun = discretetime(maxt,state,parameters)
  #store
  onemonth[simulation]= df_stochmodelrun$I[30]+df_stochmodelrun$R[30]
  Smatrix[simulation,] = df_stochmodelrun$S
  Imatrix[simulation, ] = df_stochmodelrun$I
  Rmatrix[simulation, ] = df_stochmodelrun$R
}



colMeans(Smatrix)
rm(df_stochmodelrun,parameters,maxt,state,times,discretetime,plotsirD,plotsirS)

set.seed(666)
## Day 2 practice1
t1 <- rexp(50000, rate = 1)
t2 <- rexp(50000, rate = 1)
t3 <- t1 +t2
hist(t1)
hist(t3)
t4 <- rexp(50000, rate = 0.5)
t5 <- rexp(50000, rate = 0.5)
t6 <- rexp(50000, rate = 0.5)

t7 <- t4+t5+t6
hist(t4)
hist(t7)
plot(density(t3))
mean(t3)
s

for (sim in 1:50){
  tsim <- rexp(50000, rate = 1)
  tsum <- tsim +tsim
  tsum
}

plot(density(tsum))

set(123)
#practice 2
NSIM=10000
b=0.3
t= seq(0.01, 8.00, 0.01)
hazard = 0*t
for(i in 1:length(t))
{
  if(t[i] <= 2) hazard[i] = b* t[i]
  if(t[i] > 2) hazard[i] = 2*b-b*( t[i] -2)/3
}
survival = exp(-0.01*cumsum(hazard))
  
  fdensity = hazard*survival
infectiontimes = sample(t, size = NSIM, prob = fdensity, replace = TRUE)
hist(infectiontimes)
plot(density(infectiontimes))

infinitetime =rbinom(NSIM,1, survival[length(survival)])

infectiontimes[infinitetime == 1] = 99999999
#practice 3
latenpop <- c(.15, .4, .3, .1, .05)

p=c()
p[1]= .15
p[2]= 1-.4/p[1]
p[3]= 1-.3/(p[1]*p[2])
p[4]= 1-.3/(p[1]*p[2]*p[3])
p[5]= 1-.3/(p[1]*p[2]*p[3]*p[4])
  
infectiouspop <- c(0, .04,.08
                   ,.14
                   ,.20
                   ,.18
                   ,.14
                   ,.10
                   ,.08
                   ,.04)

q=c()
q[1]=1 - infectiouspop
tinfectious = 1: length(infectiouspop)
for (k in 2:9) q[k] = 1- infectiouspop[k]/prod(q[1:(k-1)])
1- cumprod(q)

cumsum(infectiouspop)
tinfectious = 1-length(infectiouspop)
meaninfectious = sum(tinfectious*infectiouspop)/sum(infectiouspop)
R0=2
n= 50000

beta=R0/(n*meaninfectious)



library(igraph)
install.packages("ggraph")

library(ggraph)
labels <- data.frame(label = 1:6)
pairs <- data.frame(from = c(1,1,1,2,2,3,3,4,5),
                       to = c(2,3,4,3,5,4,6,5,6))

# Packages to work with and plot “graphs”
# Data for six vertices and nine edges
myedges <- c()
for(i in 1:dim(pairs)[1]) myedges <- c(myedges,pairs$from[i], pairs$to[i])
g <- make_graph(n = dim(labels)[1], edges = myedges, directed = FALSE)
V(g)$labels <- labels$label
 degree(g)
as_adjacency_matrix(g)
mynet <- ggraph(g,layout = 'stress') +
geom_edge_link(colour = rgb(0,0,0,0.5)) +
geom_node_point(colour = "red",cex=2)
plot(mynet)

# Creating a graph from those data
# Some stuff you can do with the graph
# To make a plot using ggraph
set.seed(123)
g <- sample_gnp(n = 1000, p = 3/1000, directed = FALSE)
V(g)$group <- sample(1:5, 1000, replace =TRUE)

mynet <- ggraph(g,layout = 'stress') +
  geom_edge_link(colour = rgb(0,0,0,0.5)) +
  geom_node_point(colour = "red2", cex = 2)
plot(mynet)

g1 <- sample_pa(n = 1000, m = 2, directed = FALSE)
mynet <- ggraph(g1,layout = 'stress') +
  geom_edge_link(colour = rgb(0,0,0,0.5)) +
  geom_node_point(colour = "red2", cex = 2)
plot(mynet)

nodes1 <- jefferson_labels%>%
  mutate(Male= ifelse(Male == 1, "Male", "Female"))%>%
  rename(id = SID, group = Male )


visNetwork(nodes1, jefferson_pairs)%>%
  visLegend()%>%
  visGroups(groupname = "Male", color = "blue")%>%
  visGroups(groupname = "Female", color = "red")


myedges <- c()
for(i in 1:dim(jefferson_pairs)[1]) myedges <- c(myedges, jefferson_pairs$from[i], jefferson_pairs$to[i])
g <- make_graph(n= dim(jefferson_labels)[1], edges =myedges, directed = F)
V(g)$sid <- jefferson_labels$SID
V(g)$male <- jefferson_labels$Male
V(g)$colour <- c("pink", "blue")[V(g)$male+1]
mynet_j <- ggraph(g, layout= "stress")+
  geom_edge_link(colour = rgb(0,0,0,0.5)) +
  geom_node_point(colour = V(g)$colour, cex = 2)
plot(mynet_j)  

mean(degree(g)[V(g)$male == 1])
mean(degree(g)[V(g)$male == 0])
m_degree <- degree(g)[V(g)$male==1]
f_degree <- degree(g)[V(g)$male==0]
mean(m_degree>=3)








