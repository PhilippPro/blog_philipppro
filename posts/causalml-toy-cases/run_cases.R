# Toy cases: when do causal ML models beat LM / GLM / GBM?
# Usage: Rscript run_cases.R [n_reps=30]   (simulate + summarise)
#        Rscript run_cases.R summarise     (only rebuild summary/plots from results/raw_results.csv)
# Needs: grf, xgboost, ggplot2

args <- commandArgs(trailingOnly = TRUE)
SUMMARISE_ONLY <- length(args) >= 1 && args[1] == "summarise"
R_REPS <- if (length(args) >= 1 && !SUMMARISE_ONLY) as.integer(args[1]) else 30L
BUDGET <- 0.3        # share of customers we can treat (for the targeting metric)
N_TEST <- 5000

source("R/dgp.R")
source("R/estimators.R")
source("R/evaluation.R")

source("R/cases.R")

top_share <- function(score, tau, share) {
  sel <- rank(-score, ties.method = "first") <= share * length(score)
  mean(tau[sel])
}

run_rep <- function(case, rep_id, case_idx) {
  seed <- 1000 * case_idx + rep_id
  tr <- case$dgp(case$n, seed)
  te <- case$dgp(N_TEST, seed + 1e6)
  newX <- te$data[setdiff(names(te$data), c("Y", "T"))]
  constant <- sd(te$tau) < 1e-12
  # truth-free evaluation inputs (computed once per rep, shared by all methods)
  hr <- holdout_residuals(tr$data, te$data, tr$family, seed)
  do.call(rbind, lapply(names(case$methods), function(m) {
    out <- case$methods[[m]](tr$data, newX, tr$family, seed)
    # Qini / calibration are only valid on randomised hold-out data
    qs <- if (case$randomized) qini_stats(out$score, te$data$Y, te$data$T) else list(qini = NA_real_, obs_top = NA_real_)
    se <- if (is.null(out$ate_se)) NA_real_ else out$ate_se
    data.frame(
      case = case$id, rep = rep_id, method = m,
      ate_err = mean(out$tau) - mean(te$tau),
      # 95% confidence interval of the ATE: standard error and whether the interval contains the true ATE
      ate_se = se, ate_cover = abs(mean(out$tau) - mean(te$tau)) <= qnorm(0.975) * se,
      cate_rmse = sqrt(mean((out$tau - te$tau)^2)),
      # RMSE of the predicted outcome at the observed T on the test set (what a predictive model is judged on)
      pred_rmse = if (is.null(out$mu0)) NA_real_
                  else sqrt(mean((te$data$Y - ifelse(te$data$T == 1, out$mu1, out$mu0))^2)),
      # mean true effect among the top-BUDGET share we would treat; NA if effect is constant
      policy = if (constant) NA_real_ else top_share(out$score, te$tau, BUDGET),
      oracle = if (constant) NA_real_ else top_share(te$tau, te$tau, BUDGET),
      # truth-free evaluation metrics (see R/evaluation.R)
      r_loss = r_loss(out$tau, hr),
      qini = qs$qini, obs_top30 = qs$obs_top,
      calib_slope = if (case$randomized) calib_slope(out$tau, te$data$Y, te$data$T) else NA_real_
    )
  }))
}

if (SUMMARISE_ONLY) {
  res <- read.csv("results/raw_results.csv")
} else {
  t0 <- Sys.time()
  res <- do.call(rbind, lapply(seq_along(cases), function(i) {
    message("Case ", cases[[i]]$id, " (", R_REPS, " reps)")
    do.call(rbind, lapply(seq_len(R_REPS), run_rep, case = cases[[i]], case_idx = i))
  }))
  message("Done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  dir.create("results", showWarnings = FALSE)
  write.csv(res, "results/raw_results.csv", row.names = FALSE)
}
res$case <- factor(res$case, levels = vapply(cases, `[[`, "", "id"))
res$method <- factor(res$method, levels = unique(res$method))

# Summary: bias and sd of ATE estimate, CATE RMSE, targeting value
# Means over successful runs; an estimator that returns NA in some runs (e.g. AIPW under extreme
# propensities) is counted in n_failed instead of turning the whole mean into NA.
mean_ok <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
summ <- aggregate(cbind(ate_err, cate_rmse, pred_rmse, policy, oracle) ~ case + method, data = res,
                  FUN = mean_ok, na.action = na.pass)
fails <- aggregate(ate_err ~ case + method, data = res, na.action = na.pass,
                   FUN = function(x) if (all(is.na(x))) 0L else sum(is.na(x)))
summ$n_failed <- fails$ate_err[match(paste(summ$case, summ$method), paste(fails$case, fails$method))]
sds <- aggregate(ate_err ~ case + method, data = res, FUN = sd)
summ$ate_sd <- sds$ate_err[match(paste(summ$case, summ$method), paste(sds$case, sds$method))]
summ <- summ[order(summ$case, summ$method),
             c("case", "method", "ate_err", "ate_sd", "cate_rmse", "pred_rmse", "policy", "oracle", "n_failed")]
names(summ)[3:8] <- c("ATE_bias", "ATE_sd", "CATE_RMSE", "Y_pred_RMSE", "top30_uplift", "oracle_top30")
# Monte Carlo standard errors of the means (sd over successful runs / sqrt(number of runs)); appended so the
# columns above keep their positions
mc_se <- function(x) { x <- x[!is.na(x)]; if (length(x) < 2) NA_real_ else sd(x) / sqrt(length(x)) }
se_cols <- c(ATE_bias = "ate_err", CATE_RMSE = "cate_rmse", top30_uplift = "policy")
key <- paste(summ$case, summ$method)
# Coverage of the 95% confidence interval (share of runs in which it contains the true ATE) and mean SE;
# only for estimators that return a standard error (older raw results without these columns: NA)
for (cl in c("ate_cover", "ate_se")) if (!cl %in% names(res)) res[[cl]] <- NA_real_
cov <- aggregate(res[c("ate_cover", "ate_se")], by = res[c("case", "method")],
                 FUN = function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
summ$ATE_ci_coverage <- cov$ate_cover[match(key, paste(cov$case, cov$method))]
summ$ATE_se_mean <- cov$ate_se[match(key, paste(cov$case, cov$method))]
for (nm in names(se_cols)) {
  se <- aggregate(res[se_cols[[nm]]], by = res[c("case", "method")], FUN = mc_se)
  summ[[paste0(nm, "_se")]] <- se[[3]][match(key, paste(se$case, se$method))]
}
write.csv(summ, "results/summary.csv", row.names = FALSE)
options(width = 200)
for (cs in levels(summ$case)) {
  cat("\n=== Case", cs, "===\n")
  s <- summ[summ$case == cs, -1]
  if (all(is.na(s$top30_uplift))) s <- s[, !names(s) %in% c("top30_uplift", "oracle_top30")]
  print(format(s, digits = 2, nsmall = 2), row.names = FALSE)
}

# Plots
library(ggplot2)
theme_set(theme_minimal(base_size = 11))
p1 <- ggplot(res, aes(method, ate_err)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_boxplot(fill = "#8ecae6") + coord_flip() + facet_wrap(~case, scales = "free", ncol = 2) +
  labs(x = NULL, y = "ATE error (estimate - truth)", title = "Average effect: bias and variance over simulations")
ggsave("results/ate_error.png", p1, width = 11, height = 9, dpi = 120)

p2 <- ggplot(res[!is.na(res$cate_rmse), ], aes(method, cate_rmse)) +
  geom_boxplot(fill = "#ffb703") + coord_flip() + facet_wrap(~case, scales = "free", ncol = 2) +
  labs(x = NULL, y = "CATE RMSE (lower = better)", title = "Individual effects: error vs. truth")
ggsave("results/cate_rmse.png", p2, width = 11, height = 9, dpi = 120)

pol <- res[!is.na(res$policy), ]
p3 <- ggplot(pol, aes(method, policy)) +
  geom_boxplot(fill = "#90be6d") + coord_flip() +
  geom_hline(aes(yintercept = oracle), data = aggregate(oracle ~ case, pol, mean),
             linetype = "dashed", colour = "red") +
  facet_wrap(~case, scales = "free", ncol = 1) +
  labs(x = NULL, y = sprintf("True mean effect among top %d%% targeted (red = oracle)", BUDGET * 100),
       title = "Targeting value")
ggsave("results/targeting.png", p3, width = 9, height = 7, dpi = 120)
