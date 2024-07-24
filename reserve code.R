# the function that simulates an SIR model:
sir <- function(beta, gamma, S0, I0, R0, times) {
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

# Note: the parameters beta and gamma and the duration and time step should be
# expressed in the same time unit.

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

sim1 <- sir_discrete_time(step = 1)
sim2 <- sir_discrete_time(step = .1)
sim3 <- sir_discrete_time(step = .01)
sim4 <- sir_discrete_time(step = .001)

with(sim1, plot(time, I, type = "l"))
with(sim2, lines(time, I, col = 2))
with(sim3, lines(time, I, col = 3))
with(sim4, lines(time, I, col = 4))

with(sim1, plot(time, p_inf, type = "l", ylim = 0:1))
with(sim2, lines(time, p_inf, col = 2))
with(sim3, lines(time, p_inf, col = 3))
with(sim4, lines(time, p_inf, col = 4))

with(sim1, plot(time, p_rec, type = "l", ylim = 0:1))
with(sim2, lines(time, p_rec, col = 2))
with(sim3, lines(time, p_rec, col = 3))
with(sim4, lines(time, p_rec, col = 4))

#######################################################################################

# IBM

sir_discrete_time_ibm <- function(beta = .004, gamma = .5, S0 = 999, I0 = 1, R0 = 0,
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
    new_inf <- rbinom(1, last_S, pi)
    new_rec <- rbinom(1, last_I, pr)
    S <- c(S, last_S - new_inf)
    I <- c(I, last_I + new_inf - new_rec)
    R <- c(R, last_R + new_rec)
    p_inf <- c(p_inf, pi)
    p_rec <- c(p_rec, pr)
  }
  
  tibble::tibble(time = c(0, times), S = S, I = I, R = R, p_inf = p_inf, p_rec = p_rec)
}

sim1 <- sir_discrete_time_ibm(step = 1)
sim2 <- sir_discrete_time_ibm(step = .1)
sim3 <- sir_discrete_time_ibm(step = .01)
sim4 <- sir_discrete_time_ibm(step = .001)

with(sim1, plot(time, I, type = "l"))
with(sim2, lines(time, I, col = 2))
with(sim3, lines(time, I, col = 3))
with(sim4, lines(time, I, col = 4))

with(sim1, plot(time, p_inf, type = "l", ylim = 0:1))
with(sim2, lines(time, p_inf, col = 2))
with(sim3, lines(time, p_inf, col = 3))
with(sim4, lines(time, p_inf, col = 4))

with(sim1, plot(time, p_rec, type = "l", ylim = 0:1))
with(sim2, lines(time, p_rec, col = 2))
with(sim3, lines(time, p_rec, col = 3))
with(sim4, lines(time, p_rec, col = 4))

#######################################################################################

f <- function(step, gamma) {
  1 - exp(- gamma * step)  
}

tmax <- 10
gamma <- .5 # /day
steps <- seq(0, tmax, .1)
probs <- purrr::map_dbl(steps, f, gamma = gamma)
plot(steps, probs, type = "l", ylim = 0:1,
     xlab = "time step (days)", ylab = "probability")
abline(v = seq(0, tmax, .5), col = "grey")
abline(v = seq(0, tmax, 1), col = "grey", lwd = 2)
abline(h = seq(0, 1, .1), col = "grey")
abline(h = 0:1, lwd = 2)
lines(steps, probs, lwd = 2, col = 4)
abline(0, gamma, lwd = 2, col = 2)

#######################################################################################

g <- function(step, beta) {
  1 - exp(- beta * max(out$I) * step)  
}

tmax <- 10
beta <- .004 # /day
steps <- seq(0, tmax, .1)
probs <- purrr::map_dbl(steps, g, beta = beta)
plot(steps, probs1, type = "n", ylim = 0:1,
     xlab = "time step (days)", ylab = "probability")
abline(v = seq(0, tmax, .5), col = "grey")
abline(v = seq(0, tmax, 1), col = "grey", lwd = 2)
abline(h = seq(0, 1, .1), col = "grey")
abline(h = 0:1, lwd = 2)
lines(steps, probs, lwd = 2, col = 4)
abline(0, beta * max(out$I), lwd = 2, col = 2)

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
    pi <- 1 - exp(- (exp(beta * last_I * step)) * step)
    pr <- 1 - exp(- (exp(gamma * step)) * step)
    pi <- beta * last_I * step
    pr <- gamma * step
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


sim1 <- sir_discrete_time(step = 1)
sim2 <- sir_discrete_time(step = .1)
sim3 <- sir_discrete_time(step = .01)
sim4 <- sir_discrete_time(step = .001)

with(out, plot(time, I, type = "l", col = 2, lwd = 3))
with(sim1, lines(time, I))
with(sim2, lines(time, I, col = 2))
with(sim3, lines(time, I, col = 3))
with(sim4, lines(time, I, col = 4))

with(sim1, plot(time, p_inf, type = "l", ylim = 0:1))
with(sim2, lines(time, p_inf, col = 2))
with(sim3, lines(time, p_inf, col = 3))
with(sim4, lines(time, p_inf, col = 4))

with(sim1, plot(time, p_rec, type = "l", ylim = 0:1))
with(sim2, lines(time, p_rec, col = 2))
with(sim3, lines(time, p_rec, col = 3))
with(sim4, lines(time, p_rec, col = 4))
