# Response surfaces for the Engle-Granger statistic, from the files of
# eg_simulate.R; fit by fit_surface() (surface_fit.R, copied from
# unitroottools) on T >= 25, separately for each nvars and Step 1 case:
#   q_p(T, k) = theta_inf(p) + a_k(p) / T + b_k(p) / T^2 + c_k(p) / T^3,
# T = length of the series. Cells with no residual degrees of freedom (NA)
# are dropped. Writes eg_rs.rds: a list indexed by nvars ("2", "3", ...) and
# case ("n", "c", "ct"); then run make_sysdata.R.
#
# Check against MacKinnon (1996): his numerical distribution functions for
# niv = nvars integrated variables, as shipped in urca (internal .urcval();
# urca only needs to be installed here, it is not a package dependency).
source("surface_fit.R")

files <- sort(list.files(pattern = "^eg_quantiles_n[0-9]+\\.csv\\.gz$"))
if (!length(files)) stop("No eg_quantiles_n*.csv.gz files: run eg_simulate.R first.")
eg <- do.call(rbind, lapply(files, read.csv))
eg <- eg[!is.na(eg$q), ]
nvars_all <- sort(unique(eg$nvars))

eg_rs <- lapply(stats::setNames(nvars_all, nvars_all), function(nv) {
  lapply(c(n = "n", c = "c", ct = "ct"), function(cs) {
    rs <- fit_surface(eg[eg$nvars == nv & eg$case == cs, ])
    report_fit(rs, sprintf("n=%d %s", nv, cs), probs = c(0.01, 0.05, 0.10))
    rs
  })
})
saveRDS(eg_rs, "eg_rs.rds")

# The p-values interpolate between the fitted quantiles, so these must
# increase with p in every (T, k) cell of the design: count the exceptions.
source("../R/calibration.R")
cat("\nNon-monotone fitted quantiles (should be 0 everywhere):\n")
for (nv in names(eg_rs)) for (cs in names(eg_rs[[nv]])) {
  rs <- eg_rs[[nv]][[cs]]
  bad <- 0L; n <- 0L
  for (T in c(25, 30, 40, 50, 75, 100, 150, 250, 500, 1000)) for (k in rs$ks) {
    if (T < rs$design$Tmin[match(k, rs$design$k)]) next
    n <- n + 1L
    if (is.unsorted(rs_quantiles(rs, T, k))) bad <- bad + 1L
  }
  cat(sprintf("  n=%s %-2s %d of %d cells\n", nv, cs, bad, n))
}

if (requireNamespace("urca", quietly = TRUE)) {
  urcval <- get(".urcval", asNamespace("urca"))
  itv <- c(n = 1L, c = 2L, ct = 3L)
  cat("\nOwn surfaces (k = 0) vs MacKinnon (1996), critical values at 1%, 5%, 10%\n")
  worst <- 0
  for (nv in nvars_all) for (cs in names(itv)) for (T in c(50, 100, 250, 500, 1000)) {
    rs <- eg_rs[[as.character(nv)]][[cs]]
    q <- drop(rs$intercept + rs$slopes["0", , ] %*% c(1 / T, 1 / T^2, 1 / T^3))
    own <- q[match(c(0.01, 0.05, 0.10), round(rs$probs, 6))]
    mk <- vapply(c(0.01, 0.05, 0.10), function(a) urcval(a, nobs = T, niv = nv, itt = 1,
                                                          itv = itv[[cs]], nc = 1), numeric(1))
    worst <- max(worst, abs(own - mk))
    cat(sprintf("n=%d %-2s T=%4d  own %s   MacKinnon %s   max diff %.4f\n", nv, cs, T,
                paste(sprintf("%7.3f", own), collapse = ""),
                paste(sprintf("%7.3f", mk), collapse = ""), max(abs(own - mk))))
  }
  cat(sprintf("Largest difference: %.4f\n", worst))
}
