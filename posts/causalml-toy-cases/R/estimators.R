# Estimators. Common interface:
#   f(train, newX, family, seed) -> list(tau = CATE estimate per newX row (NA if not an effect model),
#                                        score = targeting score, higher = treat first,
#                                        mu0, mu1 = predicted outcome under T=0 / T=1 (optional, for prediction metric),
#                                        ate_se = standard error of the ATE estimate (optional; only where a standard
#                                                 formula exists, used for the confidence-interval coverage))
# `train` has columns Y, T and covariates; `newX` has covariates only.

covs <- function(df) setdiff(names(df), c("Y", "T"))

est_naive <- function(train, newX, family, seed) {
  m1 <- mean(train$Y[train$T == 1]); m0 <- mean(train$Y[train$T == 0])
  n <- nrow(newX)
  se <- sqrt(var(train$Y[train$T == 1]) / sum(train$T == 1) + var(train$Y[train$T == 0]) / sum(train$T == 0))
  list(tau = rep(m1 - m0, n), score = rep(m1 - m0, n), mu0 = rep(m0, n), mu1 = rep(m1, n), ate_se = se)
}

# Heteroskedasticity-robust (HC1, sandwich) standard error of coefficient j: X design, e residuals
hc1_se <- function(X, e, j) {
  bread <- solve(crossprod(X))
  V <- bread %*% crossprod(X * e) %*% bread * nrow(X) / (nrow(X) - ncol(X))
  sqrt(V[j, j])
}

# LM / GLM with main effect of T only (constant effect, adjusted for covariates)
est_glm_main <- function(train, newX, family, seed) {
  m <- glm(reformulate(c("T", covs(train)), "Y"), data = train, family = family)
  out <- glm_tau(m, newX)
  # robust SE of the T coefficient = SE of the ATE in the linear model; for a logit the coefficient is not the ATE
  if (family == "gaussian") out$ate_se <- hc1_se(model.matrix(m), residuals(m, "response"), "T")
  out
}

# LM / GLM with T x covariates interactions (linear effect heterogeneity)
est_glm_int <- function(train, newX, family, seed) {
  m <- glm(reformulate(paste0("T*(", paste(covs(train), collapse = "+"), ")"), "Y"),
           data = train, family = family)
  glm_tau(m, newX)
}

glm_tau <- function(m, newX) {
  mu1 <- unname(predict(m, cbind(newX, T = 1), type = "response"))
  mu0 <- unname(predict(m, cbind(newX, T = 0), type = "response"))
  list(tau = mu1 - mu0, score = mu1 - mu0, mu0 = mu0, mu1 = mu1)
}

# xgboost with 5-fold CV early stopping (fair-ish tuning, not default params)
fit_xgb <- function(X, y, family, seed, w = NULL) {
  set.seed(seed)
  d <- xgboost::xgb.DMatrix(as.matrix(X), label = y, weight = w)
  p <- list(objective = if (family == "binomial") "binary:logistic" else "reg:squarederror",
            max_depth = 3, eta = 0.1, subsample = 0.8, nthread = 1)
  cv <- xgboost::xgb.cv(p, d, nrounds = 300, nfold = 5, early_stopping_rounds = 10, verbose = FALSE)
  metric <- cv$evaluation_log[[grep("^test_.*_mean$", names(cv$evaluation_log))]]
  xgboost::xgb.train(p, d, nrounds = which.min(metric))
}
pred_xgb <- function(m, X) predict(m, as.matrix(X))

# Generic GBM learners, parameterised by the fit/predict pair (xgboost or LightGBM)
gbm_s <- function(fit, pred) function(train, newX, family, seed) {
  cols <- c(covs(train), "T")
  m <- fit(train[cols], train$Y, family, seed)
  mu1 <- pred(m, cbind(newX, T = 1)[cols]); mu0 <- pred(m, cbind(newX, T = 0)[cols])
  list(tau = mu1 - mu0, score = mu1 - mu0, mu0 = mu0, mu1 = mu1)
}

gbm_t <- function(fit, pred) function(train, newX, family, seed) {
  cv <- covs(train)
  m1 <- fit(train[train$T == 1, cv, drop = FALSE], train$Y[train$T == 1], family, seed)
  m0 <- fit(train[train$T == 0, cv, drop = FALSE], train$Y[train$T == 0], family, seed)
  mu1 <- pred(m1, newX[cv]); mu0 <- pred(m0, newX[cv])
  list(tau = mu1 - mu0, score = mu1 - mu0, mu0 = mu0, mu1 = mu1)
}

# "Pure" GBM as used in practice: model the outcome from X only, ignore T, target the highest risk.
# Not an effect estimate -> tau = NA; only targeting score and outcome prediction are evaluated.
gbm_pure <- function(fit, pred) function(train, newX, family, seed) {
  cv <- covs(train)
  m <- fit(train[cv], train$Y, family, seed)
  mu <- pred(m, newX[cv])
  list(tau = rep(NA_real_, nrow(newX)), score = -mu, mu0 = mu, mu1 = mu)
}

# LightGBM with 5-fold CV early stopping; same depth/learning rate as the xgboost setup,
# otherwise library defaults (incl. min_data_in_leaf = 20, which matters for n = 200)
fit_lgb <- function(X, y, family, seed) {
  p <- list(objective = if (family == "binomial") "binary" else "regression",
            learning_rate = 0.1, num_leaves = 8, max_depth = 3, bagging_fraction = 0.8,
            bagging_freq = 1, num_threads = 1, verbose = -1, seed = seed)
  d <- lightgbm::lgb.Dataset(as.matrix(X), label = y)
  cv <- lightgbm::lgb.cv(p, d, nrounds = 300, nfold = 5, early_stopping_rounds = 10, verbose = -1)
  lightgbm::lgb.train(p, d, nrounds = max(cv$best_iter, 1L), verbose = -1)
}
pred_lgb <- function(m, X) predict(m, as.matrix(X))

est_gbm_s <- gbm_s(fit_xgb, pred_xgb)
est_gbm_t <- gbm_t(fit_xgb, pred_xgb)
est_gbm_pure <- gbm_pure(fit_xgb, pred_xgb)
est_lgb_s <- gbm_s(fit_lgb, pred_lgb)
est_lgb_t <- gbm_t(fit_lgb, pred_lgb)
est_lgb_pure <- gbm_pure(fit_lgb, pred_lgb)

# Causal forest (grf), default settings (honest trees, propensity + outcome nuisance via forests)
est_causal_forest <- function(train, newX, family, seed) {
  cv <- covs(train)
  cf <- grf::causal_forest(as.matrix(train[cv]), train$Y, train$T, num.trees = 1000, seed = seed)
  tau <- predict(cf, as.matrix(newX[cv]))$predictions
  list(tau = tau, score = tau)
}

# Restrict an estimator to a subset of covariates (used to show the effect of domain knowledge)
restrict <- function(f, cols) {
  function(train, newX, family, seed) {
    keep <- c("Y", "T", cols)
    f(train[keep], newX[cols], family, seed)
  }
}

# Instrumental forest (grf): needs an instrument column Z; identifies the effect despite unobserved confounding.
est_iv_forest <- function(train, newX, family, seed) {
  cv <- setdiff(covs(train), "Z")
  f <- grf::instrumental_forest(as.matrix(train[cv]), train$Y, train$T, train$Z,
                                num.trees = 1000, seed = seed)
  tau <- predict(f, as.matrix(newX[cv]))$predictions
  list(tau = tau, score = tau)
}

# GLM with T x (all main effects + all pairwise products): the "brute force" parametric answer
# to unknown interactions. 20 covariates -> ~420 parameters.
est_glm_int2 <- function(train, newX, family, seed) {
  m <- glm(reformulate(paste0("T*(", paste(covs(train), collapse = "+"), ")^2"), "Y"),
           data = train, family = family)
  glm_tau(m, newX)
}

# Two-stage least squares with instrument Z (the classical, linear IV estimator). Constant effect.
est_2sls <- function(train, newX, family, seed) {
  cv <- setdiff(covs(train), "Z")
  train$T_hat <- fitted(lm(reformulate(c("Z", cv), "T"), data = train))
  s2 <- lm(reformulate(c("T_hat", cv), "Y"), data = train)
  tau <- unname(coef(s2)["T_hat"])
  # IV standard error: residuals use the actual T, not the first-stage prediction
  Xh <- model.matrix(s2); X <- Xh; X[, "T_hat"] <- train$T
  se <- hc1_se(Xh, train$Y - drop(X %*% coef(s2)), "T_hat")
  list(tau = rep(tau, nrow(newX)), score = rep(tau, nrow(newX)), ate_se = se)
}

# Double ML, partially linear model (Chernozhukov et al. 2018): cross-fitted xgboost nuisances
# l(X) = E[Y | X] and m(X) = E[T | X]; effect = slope of the Y residuals on the T residuals. Constant effect.
dml_residuals <- function(train, family, seed, k = 5) {
  cv <- covs(train)
  set.seed(seed)
  fold <- sample(rep_len(seq_len(k), nrow(train)))
  rY <- rT <- numeric(nrow(train))
  for (f in seq_len(k)) {
    tr <- fold != f; te <- fold == f
    mY <- fit_xgb(train[tr, cv, drop = FALSE], train$Y[tr], family, seed + f)
    mT <- fit_xgb(train[tr, cv, drop = FALSE], train$T[tr], "binomial", seed + f)
    rY[te] <- train$Y[te] - pred_xgb(mY, train[te, cv, drop = FALSE])
    rT[te] <- train$T[te] - pred_xgb(mT, train[te, cv, drop = FALSE])
  }
  list(rY = rY, rT = rT)
}

est_dml <- function(train, newX, family, seed, k = 5) {
  r <- dml_residuals(train, family, seed, k)
  tau <- sum(r$rT * r$rY) / sum(r$rT^2)
  se <- sqrt(sum(r$rT^2 * (r$rY - tau * r$rT)^2)) / sum(r$rT^2)   # Chernozhukov et al. 2018, robust form
  list(tau = rep(tau, nrow(newX)), score = rep(tau, nrow(newX)), ate_se = se)
}

# R-learner (Nie & Wager 2021): same cross-fitted residuals as Double ML, but tau(x) is a model.
# argmin_tau sum (rY - tau(X) rT)^2  <=>  regress rY / rT on X with weights rT^2 (xgboost, same settings).
est_rlearner <- function(train, newX, family, seed, k = 5) {
  cv <- covs(train)
  r <- dml_residuals(train, family, seed, k)
  m <- fit_xgb(train[cv], r$rY / r$rT, "gaussian", seed, w = r$rT^2)
  tau <- pred_xgb(m, newX[cv])
  list(tau = tau, score = tau)
}

# AIPW / doubly robust ATE from a causal forest (grf::average_treatment_effect). Constant effect.
est_aipw_cf <- function(train, newX, family, seed) {
  cv <- covs(train)
  cf <- grf::causal_forest(as.matrix(train[cv]), train$Y, train$T, num.trees = 1000, seed = seed)
  ate <- grf::average_treatment_effect(cf, target.sample = "all")
  tau <- unname(ate["estimate"])
  list(tau = rep(tau, nrow(newX)), score = rep(tau, nrow(newX)), ate_se = unname(ate["std.err"]))
}
