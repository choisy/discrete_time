# the function that simulates an SIR model:
sir <- function(beta, gamma, S0 = 999, I0 = 1, R0 = 0, times) {
  require(deSolve) # for the "ode" function
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -beta * I * S
      dI <-  beta * I * S - gamma * I
      dR <-  gamma * I
      return(list(c(dS, dI, dR)))
    })
  }
  
  # the parameters values:
  parameters_values <- c(beta  = beta, gamma = gamma)
  
  # the initial values of variables:
  initial_values <- c(S = S0, I = I0, R = R0)
  
  # solving
  out <- ode(initial_values, times, sir_equations, parameters_values)
  
  # returning the output:
  as.data.frame(out)
}

# running the model:
out <- sir(beta = 0.004, gamma = 0.5, S0 = 999, I0 = 1, R0 = 0, times = seq(0, 10, le = 500))

# plotting the output:
with(out, plot(time, I, type = "l", col = 2, lwd = 3))

#######################################################################################

sir_discrete_time <- function(beta = .004, gamma = .5, S0 = 999, I0 = 1, R0 = 0,
                              duration = 10, step = 1) {
  last <- function(x) tail(x, 1)
  times <- seq(step, duration, step)
  
  S <- S0
  I <- I0
  R <- R0
  p_inf <- NA
  p_rec <- NA
  
  for (i in times) {
    last_S <- last(S)
    last_I <- last(I)
    last_R <- last(R)
    pi <- 1 - exp(- beta * last_I * step)
    pr <- 1 - exp(- gamma * step)
    new_inf <- pi * last_S
    new_rec <- pr * last_I
    S <- c(S, last_S - new_inf)
    I <- c(I, last_I + new_inf - new_rec)
    R <- c(R, last_R + new_rec)
    p_inf <- c(p_inf, pi)
    p_rec <- c(p_rec, pr)
  }
  
  tibble::tibble(time = c(0, times), S = S, I = I, R = R, p_inf = p_inf, p_rec = p_rec)
}

dbinom2 <- function(...)  dbinom(..., log = TRUE)

mLL <- function(par) {
  dis <- sir_discrete_time(beta = par[1], gamma = par[2], step = 1)
  cont <- sir(beta = .004, gamma = .5, times = 0:10)
  - sum(dbinom2(round(tail(dis$S, -1)), 1000, tail(cont$S, -1) / 1000)) -
    sum(dbinom2(round(tail(dis$I, -1)), 1000, tail(cont$I, -1) / 1000)) -
    sum(dbinom2(round(tail(dis$R, -1)), 1000, tail(cont$R, -1) / 1000))
}

mLL(c(.005, .5))

par_dis <- optim(c(.004, .5), mLL)$par
out_dis <- sir_discrete_time(beta = par_dis[1], gamma = par_dis[2])

with(out, plot(time, I, type = "l", col = 2, lwd = 3, ylim = c(0, 1000)))
with(out_dis, lines(time, I, type = "o"))
