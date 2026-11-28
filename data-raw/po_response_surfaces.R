# Calibration of the Phillips-Ouliaris statistics Z_t and Z_alpha, from
# po_quantiles.csv.gz (po_simulate.R), for each nvars and case, with k = 0.
# Z_t: response surface by fit_surface() (surface_fit.R) on T >= 25,
#   q(T) = theta_inf(prob) + a(prob) / T + b(prob) / T^2 + c(prob) / T^3
# (fit residuals <= 0.02). Z_alpha: the cubic in 1/T does not fit (residuals
# up to 0.46, about 10 standard errors, with a smooth pattern: the Andrews
# bandwidth grows like T^(1/3), so the bias is not a polynomial in 1/T), so
# the simulated quantiles are tabulated by T (table_calibration(), as the
# Ng-Perron statistics in unitroottools) and interpolated in 1/T; beyond
# T = 1000, towards the limit of the cubic fitted on T >= 100 (residuals
# within about one standard error there).
# Writes po_rs.rds (list by nvars, case, statistic); then run make_sysdata.R.
# Report: fit, monotonicity, and the 1/5/10% values at T = 50, 100, 250, 1000
# and the limit, to compare with Phillips and Ouliaris (1990), Tables Ib-Ic
# and IIb-IIc, and (Z_t) with the Engle-Granger limits.
source("surface_fit.R")
po <- read.csv("po_quantiles.csv.gz")
po_rs <- lapply(split(po, po$nvars), function(dn) {
  lapply(c(n = "n", c = "c", ct = "ct"), function(cs) {
    lapply(c(Zt = "Zt", Za = "Za"), function(s) {
      dd <- dn[dn$case == cs & dn$stat == s, ]
      if (s == "Zt") {
        rs <- fit_surface(dd)
        report_fit(rs, sprintf("n=%d %-2s %s", dn$nvars[1], cs, s), probs = c(0.01, 0.05, 0.10))
      } else {
        rs <- table_calibration(dd, Tlim = 100)
        f <- rs$limit_fit[rs$limit_fit$prob %in% c(0.01, 0.05, 0.10), ]
        cat(sprintf("n=%d %-2s Za  table; limit fit max |residual| (T >= 100): %.4f\n",
                    dn$nvars[1], cs, max(f$maxabs)))
      }
      rs
    })
  })
})
saveRDS(po_rs, "po_rs.rds")

source("../R/calibration.R")
cat("\nNon-monotone fitted quantiles (should be 0):\n")
for (nv in names(po_rs)) for (cs in names(po_rs[[nv]])) for (s in c("Zt", "Za")) {
  rs <- po_rs[[nv]][[cs]][[s]]
  bad <- sum(vapply(c(25, 30, 40, 50, 75, 100, 150, 250, 500, 1000),
                    function(T) is.unsorted(rs_quantiles(rs, T, 0)), logical(1)))
  if (bad) cat(sprintf("  n=%s %s %s: %d\n", nv, cs, s, bad))
}
cat("\nCritical values (1%, 5%, 10%)\n")
for (nv in names(po_rs)) for (cs in names(po_rs[[nv]])) for (s in c("Zt", "Za")) {
  rs <- po_rs[[nv]][[cs]][[s]]
  for (T in c(50, 100, 250, 1000, Inf)) {
    cv <- rs_critical_values(rs, if (is.finite(T)) T else 1e9, 0)
    cat(sprintf("n=%s %-2s %-2s T=%5s %s\n", nv, cs, s, if (is.finite(T)) T else "inf",
                paste(sprintf("%8.3f", cv), collapse = "")))
  }
}
