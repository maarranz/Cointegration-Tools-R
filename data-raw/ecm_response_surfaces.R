# Response surfaces for the ECM test statistic, from the files of
# ecm_simulate.R; fit by fit_surface() (surface_fit.R) on T >= 25, separately
# for each nvars and case ("c", "ct"), with p in the role of k:
#   q_p(T, p) = theta_inf(prob) + a_p(prob) / T + b_p(prob) / T^2 + c_p(prob) / T^3.
# Cells with fewer than 10 residual degrees of freedom are dropped (see below).
# Writes ecm_rs.rds (list by nvars, then case); then run make_sysdata.R.
# Report: fit, monotonicity, and the 1/5/10% values at T = 50, 100, 250,
# 1000 and the limit for each p, to compare with Ericsson and MacKinnon (2002)
# and with Table 2 of the cointegration chapter.
source("surface_fit.R")

files <- sort(list.files(pattern = "^ecm_quantiles_n[0-9]+\\.csv\\.gz$"))
if (!length(files)) stop("No ecm_quantiles_n*.csv.gz files: run ecm_simulate.R first.")
ecm <- do.call(rbind, lapply(files, read.csv))
ecm <- ecm[!is.na(ecm$q), ]
# Cells with very few residual degrees of freedom (large nvars and p at the
# smallest T) have t ratios with extremely heavy tails (quantiles of -70 with
# 1 df) that the polynomial in 1/T cannot fit; they are not used. With
# P = deterministics + y_{t-1} + 2(n-1) + p n columns and T - p - 1 rows, keep
# cells with at least 10 residual degrees of freedom; the design Tmin of each
# p then rises accordingly, and ecm_test() notes samples below it.
ecm$df <- with(ecm, T - k - 1 - ((case == "ct") + 2 + 2 * (nvars - 1) + k * nvars))
ecm <- ecm[ecm$df >= 10, ]
nvars_all <- sort(unique(ecm$nvars))

ecm_rs <- lapply(stats::setNames(nvars_all, nvars_all), function(nv) {
  lapply(c(c = "c", ct = "ct"), function(cs) {
    rs <- fit_surface(ecm[ecm$nvars == nv & ecm$case == cs, ])
    report_fit(rs, sprintf("n=%d %s", nv, cs), probs = c(0.01, 0.05, 0.10))
    rs
  })
})
saveRDS(ecm_rs, "ecm_rs.rds")

# The p-values interpolate between the fitted quantiles, so these must
# increase with p in every (T, k) cell of the design: count the exceptions.
source("../R/calibration.R")
cat("\nNon-monotone fitted quantiles (should be 0 everywhere):\n")
for (nv in names(ecm_rs)) for (cs in names(ecm_rs[[nv]])) {
  rs <- ecm_rs[[nv]][[cs]]
  bad <- 0L; n <- 0L
  for (T in c(25, 30, 40, 50, 75, 100, 150, 250, 500, 1000)) for (k in rs$ks) {
    if (T < rs$design$Tmin[match(k, rs$design$k)]) next
    n <- n + 1L
    if (is.unsorted(rs_quantiles(rs, T, k))) bad <- bad + 1L
  }
  cat(sprintf("  n=%s %-2s %d of %d cells\n", nv, cs, bad, n))
}


cat("\nCritical values (1%, 5%, 10%) by p\n")
for (nv in names(ecm_rs)) for (cs in names(ecm_rs[[nv]])) {
  rs <- ecm_rs[[nv]][[cs]]
  for (p in rs$ks) for (T in c(50, 100, 250, 1000, Inf)) {
    if (is.finite(T) && T < rs$design$Tmin[match(p, rs$design$k)]) next
    cv <- rs_critical_values(rs, if (is.finite(T)) T else 1e9, p)
    cat(sprintf("n=%s %-2s p=%d T=%5s  %s\n", nv, cs, p, if (is.finite(T)) T else "inf",
                paste(sprintf("%7.3f", cv), collapse = "")))
  }
}
