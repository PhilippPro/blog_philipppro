# Evaluation tools that work WITHOUT knowing the true effect (what one can do in real data).
# They use only observed Y, T, X on a hold-out set. In the simulation we compare their verdict
# with the truth-based metrics (CATE-RMSE) to see whether they would have picked the right method.

# Cross-fitting is replaced by a split: nuisances l(X) = E[Y|X], m(X) = E[T|X] are fitted on the
# training data and evaluated on the hold-out set (needs fit_xgb / pred_xgb from estimators.R).
holdout_residuals <- function(train, test, family, seed) {
  cv <- setdiff(covs(train), "Z")   # the instrument is not a covariate
  mY <- fit_xgb(train[cv], train$Y, family, seed + 501)
  mT <- fit_xgb(train[cv], train$T, "binomial", seed + 502)
  m_hat <- pred_xgb(mT, test[cv])
  list(rY = test$Y - pred_xgb(mY, test[cv]), rT = test$T - m_hat, m_hat = m_hat)
}

# R-loss (Nie & Wager): mean (rY - tau(x) rT)^2. Lower = better. NA for scores that are not effects.
r_loss <- function(tau, res) if (anyNA(tau)) NA_real_ else mean((res$rY - tau * res$rT)^2)

# Uplift (Qini) curve from a randomised hold-out: units sorted by score, observed difference in means
# among the top k. Returns the Qini coefficient (area between the curve and the random line, per unit)
# and the observed uplift in the top `share`.
qini_stats <- function(score, Y, T, share = 0.3) {
  n <- length(score)
  o <- order(score, decreasing = TRUE)   # ties keep data order (data are exchangeable)
  y <- Y[o]; t <- T[o]
  nT <- cumsum(t); nC <- cumsum(1 - t)
  up <- ifelse(nT > 0 & nC > 0, cumsum(y * t) / pmax(nT, 1) - cumsum(y * (1 - t)) / pmax(nC, 1), 0)
  gain <- up * seq_len(n) / n
  k30 <- ceiling(share * n)
  list(qini = mean(gain - seq_len(n) / n * gain[n]), obs_top = up[k30])
}

# Calibration by deciles on a randomised hold-out: slope of observed decile uplift on mean predicted effect.
# 1 = calibrated; NA when the prediction has no spread (constant effect estimators).
calib_decile <- function(tau, Y, T, groups = 10) {
  if (anyNA(tau) || length(unique(tau)) < groups) return(NULL)
  g <- ceiling(groups * rank(tau, ties.method = "first") / length(tau))
  d <- do.call(rbind, lapply(split(seq_along(tau), g), function(i)
    data.frame(pred = mean(tau[i]), obs = mean(Y[i][T[i] == 1]) - mean(Y[i][T[i] == 0]), n = length(i))))
  d$group <- seq_len(nrow(d))
  d
}
calib_slope <- function(tau, Y, T) {
  d <- calib_decile(tau, Y, T)
  if (is.null(d)) NA_real_ else unname(coef(lm(obs ~ pred, d, weights = n))["pred"])
}
