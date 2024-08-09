continuous_exponential_growth <- function(N0, r, times) {
  tibble::as_tibble(
    as.data.frame(
      ode(c(N = N0),
          times,
          function(time, variables, parameters) {
            with(as.list(c(variables, parameters)), {
              dN <- r * N
              list(dN)
            })
          },
          c(r = r))
    )
  )
}

out <- continuous_exponential_growth(10, .5, seq(0, 10, .1))
with(out, plot(time, N, type = "l", col = 4, lwd = 2))

#######################################################################################

discrete_exponential_growth <- function(N0, r, times) {
  step <- mean(diff(times))
  N <- N0
  for (i in remove_first(times)) {
    N <- c(N, last(N) * (2 - exp(- r * step)))
  }
  tibble::tibble(time = times, N = N)
}

out <- discrete_exponential_growth(10, .5, seq(0, 10, .1))
with(out, plot(time, N, type = "o", col = 4, lwd = 2))

#######################################################################################

compare_dynamics <- function(N0, r, tmin, tmax, length, step) {
  out_continuous <- continuous_exponential_growth(N0, r, seq(tmin, tmax, le = length))
  out_discrete <- discrete_exponential_growth(N0, r, seq(tmin, tmax, step))
  with(out_continuous, plot(time, N, type = "l", col = 4, lwd = 2))
  with(out_discrete, points(time, N, col = 2, lwd = 2))
}

compare_dynamics(N0 = 10, r = .5, tmin = 0, tmax = 10, length = 512, step = .1)

#######################################################################################

mLL <- function(x, r, N0, times) {
  - sum(dpois(round(unlist(discrete_exponential_growth(N0, x, times)[-1, 2])),
              unlist(continuous_exponential_growth(N0, r, times)[-1, 2]), TRUE))
}

mLL(x = .5, r = .5, N0 = 10, times = seq(0, 10, .1))

f <- function(x) {
  optimize(mLL, c(.0001, 100), r = .5, N0 = 10, times = seq(0, 10, x))$minimum
}

f(.1)

step_size <- seq(.01, .5, .01)
out <- step_size |>
  mclapply(f, mc.cores = detectCores() - 1) |> 
  unlist() |> 
  matrix(ncol = 1, byrow = TRUE) |> 
  as.data.frame() |> 
  setNames("r") %>%
  cbind(step_size, .) |> 
  tibble::as_tibble()

with(out, plot(step_size, r, col = 4, ylim = 0:1))
abline(h = .5)

#######################################################################################

compare_dynamics2 <- function(N0, rc, rd, tmin, tmax, length, step) {
  out_continuous <- continuous_exponential_growth(N0, rc, seq(tmin, tmax, le = length))
  out_discrete <- discrete_exponential_growth(N0, rd, seq(tmin, tmax, step))
  with(out_continuous, plot(time, N, type = "l", col = 4, lwd = 2))
  with(out_discrete, points(time, N, col = 2, lwd = 2))
}

walk2(out$r, out$step_size, ~ compare_dynamics2(N0 = 10, rc = .5, rd = .x, tmin = 0,
                                                tmax = 10, length = 512, step = .y))


