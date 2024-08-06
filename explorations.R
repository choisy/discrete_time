## Packages ###########################################################################

library(deSolve)
library(purrr)
library(parallel)
library(magrittr)

## Functions ##########################################################################

dbinom2 <- function(...)  dbinom(..., log = TRUE)

remove_first <- function(x) tail(x, -1)

last <- function(x) tail(x, 1)


## 1. continuous-time model ###########################################################

# the function that simulates an SIR model:
sir_continuous <- function(beta, gamma, S0, I0, R0, times) {
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      new_cases <- beta * I * S
      recovered <- gamma * I
      dS <- - new_cases
      dI <-   new_cases - recovered
      dR <-   recovered
      return(list(c(dS, dI, dR)))
    })
  }
  
  tibble::as_tibble(
    as.data.frame(
      ode(c(S = S0, I = I0, R = R0), times, sir_equations,
          c(beta  = beta, gamma = gamma))))
}


## 2. discrete-time model #############################################################

sir_discrete <- function(beta, gamma, S0, I0, R0, times) {
  
  step <- mean(diff(times))
  
  S <- S0
  I <- I0
  R <- R0
  p_inf <- NA
  p_rec <- NA
  
  for (i in remove_first(times)) {
    last_S <- last(S)
    last_I <- last(I)
    last_R <- last(R)
    new_cases <- (1 - exp(- beta * last_I * step)) * last_S
    recovered <- (1 - exp(- gamma * step)) * last_I
    S <- c(S, last_S - new_cases)
    I <- c(I, last_I + new_cases - recovered)
    R <- c(R, last_R + recovered)
  }
  
  tibble::tibble(time = times, S = S, I = I, R = R)
}

# sir_discrete(beta = .004, gamma = .5, S0 = 999, I0 = 1, R0 = 0, times = seq(0, 10, .1))

## 3. likelihood ######################################################################

mLL1 <- function(observed, expected, size) {
  - sum(dbinom2(observed, size, expected / size))
}

mLL <- function(par, beta, gamma, S0, I0, R0, size, times) {
  continuous <- sir_continuous(beta, gamma, S0, I0, R0, times)[-1, ]
  discrete <- round(sir_discrete(par[1], par[2], S0, I0, R0, times)[-1, ])
  
#  size <- S0 + I0 + R0
  mLL1(discrete$S, continuous$S, size) +
    mLL1(discrete$I, continuous$I, size) +
    mLL1(discrete$R, continuous$R, size)
}


#######################################################################################

beta <- .004
gamma <- .5
S0 <- 999
I0 <- 1
R0 <- 0
tmin <- 0
tmax <- 10

size <- S0 + I0 + R0

out_continuous <- sir_continuous(beta, gamma, S0, I0, R0, times = seq(tmin, tmax, le = 500))
with(out_continuous, plot(time, I, type = "l", col = 2, lwd = 3, ylim = c(0, size)))

step <- 1
(parameters_discrete <- optim(c(beta, gamma), mLL, beta = beta, gamma = gamma, S0 = S0,
                              I0 = I0, R0 = R0, size = size, times = seq(0, 10, step))$par)

out_discrete <- sir_discrete(parameters_discrete[1], parameters_discrete[2],
                             S0, I0, R0, seq(0, 10, step))

with(out_discrete, lines(time, I, type = "o"))


#######################################################################################

beta <- .004
gamma <- .5
S0 <- 999
I0 <- 1
R0 <- 0
tmin <- 0
tmax <- 10

size <- S0 + I0 + R0

f <- function(x) {
  optim(c(beta, gamma), mLL, beta = beta, gamma = gamma, S0 = S0, I0 = I0, R0 = R0,
        size = size, times = seq(0, 10, x))$par
}

step_size <- seq(.01, 2, .01)
out <- step_size |>
  mclapply(f, mc.cores = detectCores() - 1) |> 
  unlist() |> 
  matrix(ncol = 2, byrow = TRUE) |> 
  as.data.frame() |> 
  setNames(c("beta", "gamma")) %>%
  cbind(step_sizes, .) |> 
  tibble::as_tibble()
  
with(out, plot(step_size, beta, col = 4))
abline(0, 1)
