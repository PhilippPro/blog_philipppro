# Case registry shared by run_cases.R and report.Rmd.
# Requires dgp.R and estimators.R to be sourced.

base_methods <- list(
  "Naive (diff. in means)"  = est_naive,
  "LM/GLM (T + X)"          = est_glm_main,
  "LM/GLM (T x X)"          = est_glm_int,
  "GBM S-learner"           = est_gbm_s,
  "GBM T-learner"           = est_gbm_t,
  "LightGBM S-learner"      = est_lgb_s,
  "LightGBM T-learner"      = est_lgb_t,
  "Causal forest"           = est_causal_forest,
  "Double ML (xgboost)"     = est_dml,
  "R-learner (xgboost)"     = est_rlearner,
  "AIPW (causal forest)"    = est_aipw_cf
)

# Plain GBMs on X only (no T): outcome prediction / risk targeting, no effect estimate
pure_methods <- list(
  "GBM pure (ignores T)"      = est_gbm_pure,
  "LightGBM pure (ignores T)" = est_lgb_pure
)
all_base <- c(base_methods, pure_methods)

# Order: A simple world -> B which variables (identification) -> C how flexible (estimation)
cases <- list(
  list(id = "1 Small n", block = "A", dgp = dgp_small_n, n = 200, p = 3, methods = all_base,
       randomized = TRUE, setting = "Randomized, small sample",
       truth = "tau = 1 (constant)"),
  list(id = "2 Confounding linear", block = "B", dgp = dgp_confounded, n = 1000, p = 2, methods = all_base,
       randomized = FALSE, setting = "Observational, confounder measured, acts linearly",
       truth = "tau = 2 (constant)"),
  list(id = "3 Confounding nonlinear", block = "B", dgp = dgp_nonlinear_confounding, n = 1000, p = 2,
       methods = all_base,
       randomized = FALSE, setting = "Observational, confounder measured, acts U-shaped",
       truth = "tau = 1 (constant)"),
  list(id = "4 Mediator", block = "B", dgp = dgp_mediator, n = 1000, p = 2,
       methods = c(all_base, list(
         "LM/GLM (T + X), w/o mediator" = restrict(est_glm_main, "X1"),
         "Causal forest, w/o mediator"  = restrict(est_causal_forest, "X1"))),
       randomized = TRUE, setting = "Randomized, one variable is a mediator",
       truth = "tau = 1.5 (total effect, constant)"),
  list(id = "5 Unobs. confounding", block = "B", dgp = dgp_unobserved, n = 2000, p = 1,
       methods = c(setNames(lapply(all_base, restrict, cols = "X1"), names(all_base)),
                   list("2SLS (uses Z)" = est_2sls,
                        "Instrumental forest (uses Z)" = est_iv_forest)),
       randomized = FALSE, setting = "Observational, confounder U not measured, randomized instrument Z",
       truth = "tau = 1 (constant)"),
  list(id = "6 Linear heterogeneity", block = "C", dgp = dgp_linear_het, n = 2000, p = 3, methods = all_base,
       randomized = TRUE, setting = "Randomized",
       truth = "tau = 1 + 1.5 X1 (linear)"),
  list(id = "7 Risk ≠ uplift", block = "C", dgp = dgp_risk_vs_uplift, n = 3000, p = 2,
       methods = all_base,
       randomized = TRUE, setting = "Randomized, binary outcome (subscription kept / cancelled)",
       truth = "tau(X1) = bell shape around medium risk (max. approx. 0.27)"),
  list(id = "8 Many variables", block = "C", dgp = dgp_many_vars, n = 2000, p = 20,
       methods = c(all_base, list("LM/GLM (T x all 2-way X)" = est_glm_int2)),
       randomized = TRUE, setting = "Randomized, 20 variables (18 of them noise)",
       truth = "tau = 2 if X1 > 0 AND X2 > 0, else 0"),
  list(id = "9 Observational + heterogeneity", block = "D", dgp = dgp_obs_het, n = 2000, p = 2, methods = all_base,
       randomized = FALSE, setting = "Observational, confounder measured (U-shaped), effect varies with X2",
       truth = "tau = 1 + 1.5 X2 (linear)")
)
