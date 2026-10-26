# Helpers for report.Rmd: styled HTML tables, DAG drawing, per-case result tables and plots.

# Bootstrap-styled HTML table, auto width, cells top-aligned (see CSS in report.Rmd)
tbl <- function(df, align = NULL, caption = NULL) {
  knitr::kable(df, format = "html", escape = FALSE, align = align, caption = caption, row.names = FALSE,
               table.attr = 'class="table table-condensed table-striped"')
}

# DAG: nodes = data.frame(name, x, y, role); edges = data.frame(from, to)
role_cols <- c(Treatment = "#4361ee", Outcome = "#e63946", Confounder = "#f4a261", Collider = "#9d4edd", Mediator = "#9d4edd",
               Instrument = "#2a9d8f", Covariate = "#adb5bd", `Effect modifier` = "#ffd166",
               Unobserved = "#eef0f2")
draw_dag <- function(nodes, edges) {
  # Approximate label boxes in data units, so arrows end at the box border instead of under the label
  fw <- knitr::opts_current$get("fig.width") %||% 7; fh <- knitr::opts_current$get("fig.height") %||% 2.6
  xspan <- diff(range(nodes$x)) + 1; yspan <- diff(range(nodes$y)) + 0.8
  # labels may contain "\n": width from the longest line, height from the number of lines
  lines <- strsplit(nodes$name, "\n", fixed = TRUE)
  nodes$hw <- (vapply(lines, function(l) max(nchar(l)), 0) * 0.042 + 0.1) * xspan / (fw * 0.78)   # 0.78: share of width left of legend
  nodes$hh <- (0.17 + 0.15 * (lengths(lines) - 1)) * yspan / fh
  box <- setNames(nodes[c("name", "x", "y", "hw", "hh")], c("name", "x", "y", "hw", "hh"))
  e <- merge(merge(edges, setNames(box, c("from", "x0", "y0", "hw0", "hh0"))),
             setNames(box, c("to", "x1", "y1", "hw1", "hh1")))
  dx <- e$x1 - e$x0; dy <- e$y1 - e$y0
  s0 <- pmin(e$hw0 / abs(dx), e$hh0 / abs(dy)); s1 <- pmin(e$hw1 / abs(dx), e$hh1 / abs(dy))
  e$xs <- e$x0 + dx * s0; e$ys <- e$y0 + dy * s0
  e$xe <- e$x1 - dx * (s1 + 0.02); e$ye <- e$y1 - dy * (s1 + 0.02)
  e$lt <- ifelse(e$from %in% nodes$name[nodes$role == "Unobserved"], "dashed", "solid")
  p <- ggplot() +
    geom_segment(data = e, aes(xs, ys, xend = xe, yend = ye, linetype = lt),
                 arrow = arrow(length = unit(2.5, "mm"), type = "closed"), colour = "grey30")
  # optional edge labels (flow charts): edges$label, drawn at 60 % of the edge (spreads sibling labels)
  if ("label" %in% names(e)) {
    el <- e[!is.na(e$label), ]
    p <- p + geom_label(data = el, aes(x = xs + 0.6 * (xe - xs), y = ys + 0.6 * (ye - ys), label = label), size = 3.6,
                        fill = "white", label.size = 0, label.padding = unit(0.15, "lines"))
  }
  p +
    geom_label(data = nodes, aes(x, y, label = name, fill = role), size = 4.2, label.padding = unit(0.35, "lines"),
               label.r = unit(0.4, "lines"), key_glyph = draw_key_rect) +
    scale_fill_manual(values = role_cols, name = NULL) + scale_linetype_identity() +
    coord_cartesian(xlim = range(nodes$x) + c(-0.6, 0.9), ylim = range(nodes$y) + c(-0.4, 0.4)) +
    theme_void(base_size = 12) + theme(legend.position = "right")
}

metric_labels <- c(ATE_bias = "ATE bias", ATE_sd = "ATE SD", CATE_RMSE = "CATE RMSE",
                   Y_pred_RMSE = "Y prediction RMSE", top30_uplift = "Top-30% uplift")

# Results table of one case; best value per metric in bold, NA as en dash
result_table <- function(summ, case_id) {
  s <- summ[summ$case == case_id, ]
  oracle <- s$oracle_top30[1]
  cols <- c("ATE_bias", "ATE_sd", "CATE_RMSE", "Y_pred_RMSE", "top30_uplift")
  cols <- cols[vapply(cols, function(cl) any(!is.na(s[[cl]])), TRUE)]
  digits <- if (max(abs(s$CATE_RMSE), na.rm = TRUE) < 0.2) 3 else 2
  best <- list(ATE_bias = which.min(abs(s$ATE_bias)), ATE_sd = integer(0),  # low SD is no merit if biased
               CATE_RMSE = which.min(s$CATE_RMSE), Y_pred_RMSE = which.min(s$Y_pred_RMSE),
               top30_uplift = which.max(s$top30_uplift))
  lab <- as.character(s$method)
  nf <- if ("n_failed" %in% names(s)) s$n_failed else rep(0, nrow(s))
  lab[nf > 0] <- paste0(lab[nf > 0], " &dagger;")
  out <- data.frame(Method = lab)
  for (cl in cols) {
    v <- ifelse(is.na(s[[cl]]), "&ndash;", formatC(s[[cl]], format = "f", digits = digits))
    v <- sub("^-(0[.]0+)$", "\\1", v)  # no "-0.00"
    if (length(best[[cl]])) v[best[[cl]]] <- paste0("<b>", v[best[[cl]]], "</b>")
    out[[metric_labels[[cl]]]] <- v
  }
  cap <- "Means over all simulations. Bold = best method per column (ATE bias: closest to 0; ATE SD not marked, because a small spread is worth nothing when the bias is large)."
  if (any(nf > 0)) cap <- paste0(cap, " &dagger; No result in ", paste(nf[nf > 0], collapse = "/"),
                                 " runs (extreme propensities); mean over the remaining runs.")
  if (!is.na(oracle)) cap <- paste0(cap, " Top-30% uplift of a perfect selection (oracle): ",
                                    formatC(oracle, format = "f", digits = digits), ".")
  # Monte Carlo uncertainty: largest standard error per shown metric (all methods run on the same data sets,
  # so differences between methods are even more precise than this)
  se <- c(ATE_bias = "ATE bias", CATE_RMSE = "CATE RMSE", top30_uplift = "Top-30% uplift")
  se <- se[names(se) %in% cols & paste0(names(se), "_se") %in% names(s)]
  if (length(se)) {
    mx <- vapply(names(se), function(cl) max(s[[paste0(cl, "_se")]], na.rm = TRUE), 0)
    cap <- paste0(cap, " Monte Carlo standard errors of the means at most: ",
                  paste0(se, " ", formatC(mx, format = "f", digits = 3), collapse = ", "),
                  ". Differences below about twice that count as equally good.")
  }
  tbl(out, align = c("l", rep("r", length(cols))), caption = cap)
}

# Boxplots of one case: ATE error, CATE RMSE, targeting (if effect varies)
result_plots <- function(res, case_id) {
  r <- res[res$case == case_id, ]
  r$method <- factor(r$method, levels = rev(unique(as.character(r$method))))
  long <- rbind(
    data.frame(method = r$method, metric = "ATE error (0 = perfect)", value = r$ate_err),
    data.frame(method = r$method, metric = "CATE RMSE (0 = perfect)", value = r$cate_rmse),
    data.frame(method = r$method, metric = "Top-30% uplift (red = oracle)", value = r$policy))
  long <- long[!is.na(long$value), ]
  ref <- data.frame(metric = c("ATE error (0 = perfect)", "CATE RMSE (0 = perfect)",
                               "Top-30% uplift (red = oracle)"),
                    ref = c(0, 0, mean(r$oracle)))
  ref <- ref[ref$metric %in% long$metric & !is.na(ref$ref), ]
  ggplot(long, aes(value, method)) +
    geom_vline(data = ref, aes(xintercept = ref), colour = "#e63946", linetype = "dashed") +
    geom_boxplot(fill = "#8ecae6", outlier.size = 0.8) +
    facet_wrap(~metric, scales = "free_x", nrow = 1) + labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 13) + theme(strip.text = element_text(face = "bold"))
}
