# Data generating processes for the toy cases.
# Each returns list(data = data.frame(Y, T, <covariates>), tau = true individual effect,
#                   family = "gaussian" | "binomial")
# Convention: Y is a "good" outcome (higher = better), so tau > 0 means treatment helps.

# Confounding, constant effect, linear world.
# Customer service calls customers with low satisfaction (X1) more often -> T depends on X1, which also drives Y.
dgp_confounded <- function(n, seed) {
  set.seed(seed)
  X1 <- rnorm(n)                      # satisfaction (confounder)
  X2 <- rnorm(n)                      # pure noise
  T  <- rbinom(n, 1, plogis(-1.2 * X1))
  tau <- rep(2, n)
  Y  <- 10 + 3 * X1 + tau * T + rnorm(n)
  list(data = data.frame(Y, T, X1, X2), tau = tau, family = "gaussian")
}

# Randomized, effect varies LINEARLY in X1, linear baseline: a GLM with T x X is correctly specified.
dgp_linear_het <- function(n, seed) {
  set.seed(seed)
  X <- matrix(runif(3 * n, -1, 1), n, 3, dimnames = list(NULL, paste0("X", 1:3)))
  T <- rbinom(n, 1, 0.5)
  tau <- 1 + 1.5 * X[, "X1"]
  Y <- 2 * X[, "X2"] + 2 * X[, "X3"] + tau * T + rnorm(n, sd = 1.5)
  list(data = data.frame(Y, T, X), tau = tau, family = "gaussian")
}

# Nonlinear (U-shaped) confounding, constant effect. Very young and very old customers (|X1| large)
# are contacted more often AND have a lower outcome. Linear adjustment for X1 cannot remove this.
dgp_nonlinear_confounding <- function(n, seed) {
  set.seed(seed)
  X1 <- rnorm(n)                      # age (standardized)
  X2 <- rnorm(n)
  T  <- rbinom(n, 1, plogis(-1 + 1.2 * X1^2))
  tau <- rep(1, n)
  Y  <- 10 - 1.5 * X1^2 + X2 + tau * T + rnorm(n)
  list(data = data.frame(Y, T, X1, X2), tau = tau, family = "gaussian")
}

# Observational data AND heterogeneous effect. Same U-shaped confounding as case 3 (X1), but the effect
# varies linearly in a different variable (X2): confounder and effect modifier are separate.
dgp_obs_het <- function(n, seed) {
  set.seed(seed)
  X1 <- rnorm(n)                      # age (standardized), confounder
  X2 <- rnorm(n)                      # effect modifier
  T  <- rbinom(n, 1, plogis(-1 + 1.2 * X1^2))
  tau <- 1 + 1.5 * X2
  Y  <- 10 - 1.5 * X1^2 + X2 + tau * T + rnorm(n)
  list(data = data.frame(Y, T, X1, X2), tau = tau, family = "gaussian")
}

# Randomized, binary outcome (subscription kept). Churn risk (X1) is NOT the same as uplift:
# the effect sits at medium risk ("persuadables"); high-risk customers are lost causes,
# low-risk ones are sure things. Effect acts on the logit scale, tau is on the probability scale.
dgp_risk_vs_uplift <- function(n, seed) {
  set.seed(seed)
  X1 <- rnorm(n)                      # churn risk score (high = likely to cancel)
  X2 <- rnorm(n)                      # noise
  T  <- rbinom(n, 1, 0.5)
  lp0 <- 0.5 - 1.5 * X1
  lp1 <- lp0 + 1.5 * exp(-X1^2 / (2 * 0.4^2))
  p0 <- plogis(lp0); p1 <- plogis(lp1)
  Y <- rbinom(n, 1, ifelse(T == 1, p1, p0))
  list(data = data.frame(Y, T, X1, X2), tau = p1 - p0, family = "binomial")
}

# Randomized, constant effect, linear, but tiny sample.
dgp_small_n <- function(n, seed) {
  set.seed(seed)
  X <- matrix(rnorm(3 * n), n, 3, dimnames = list(NULL, paste0("X", 1:3)))
  T <- rbinom(n, 1, 0.5)
  tau <- rep(1, n)
  Y <- 1 + X[, "X1"] + 0.5 * X[, "X2"] + tau * T + rnorm(n, sd = 1.5)
  list(data = data.frame(Y, T, X), tau = tau, family = "gaussian")
}

# Bad control: mediator. Randomized campaign T raises the number of app visits M, which raises Y.
# Total effect of T on Y = direct 0.5 + via M 1.0 * 1.0 = 1.5. Adjusting for M leaves only the direct 0.5.
dgp_mediator <- function(n, seed) {
  set.seed(seed)
  X1 <- rnorm(n)
  T  <- rbinom(n, 1, 0.5)
  M  <- 0.5 * X1 + 1 * T + rnorm(n)
  tau <- rep(1.5, n)
  Y  <- X1 + 0.5 * T + 1 * M + rnorm(n)
  list(data = data.frame(Y, T, X1, M), tau = tau, family = "gaussian")
}

# Randomized, many covariates (20), effect only where X1 > 0 AND X2 > 0 (interaction of two
# variables, 18 noise variables), nonlinear baseline. Linear T x X cannot represent this shape.
dgp_many_vars <- function(n, seed, p = 20) {
  set.seed(seed)
  X <- matrix(rnorm(p * n), n, p, dimnames = list(NULL, paste0("X", seq_len(p))))
  T <- rbinom(n, 1, 0.5)
  tau <- 2 * (X[, "X1"] > 0 & X[, "X2"] > 0)
  mu0 <- 2 * sin(X[, "X3"]) + X[, "X4"] * X[, "X5"] + X[, "X6"]
  Y <- mu0 + tau * T + rnorm(n)
  list(data = data.frame(Y, T, X), tau = tau, family = "gaussian")
}

# Unobserved confounding. U (e.g. "motivation", never recorded) raises both the chance of treatment
# and Y. Z is a randomized encouragement (instrument): moves T, has no direct effect on Y.
# Only X1 and Z are observed. True effect = 1.
dgp_unobserved <- function(n, seed) {
  set.seed(seed)
  U  <- rnorm(n)
  X1 <- rnorm(n)
  Z  <- rbinom(n, 1, 0.5)
  T  <- rbinom(n, 1, plogis(-1.5 + 3 * Z + 1.5 * U))
  tau <- rep(1, n)
  Y  <- X1 + 2 * U + tau * T + rnorm(n)
  # placebo outcome (for the placebo test): same confounding by U, but no effect of T. Not part of `data`.
  Yp <- X1 + 2 * U + rnorm(n)
  list(data = data.frame(Y, T, X1, Z), tau = tau, family = "gaussian", placebo = Yp)
}
