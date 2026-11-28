# Finite-sample distribution of the Phillips-Ouliaris (1990) residual-based
# statistics Z_t and Z_alpha under the null of no cointegration, with the
# kernel po_kernel.cpp (Bartlett kernel, Andrews bandwidth: the defaults of
# po_test()).
#
# DGP: n = nvars independent driftless Gaussian random walks; cointegrating
# regression of y_1 on [deterministics, y_2..y_n] over t = 1..T, cases "n",
# "c", "ct", as eg_simulate.R. No lag order: the statistics correct for
# serial correlation with the kernel, so the surfaces have k = 0 only.
# T = 20..1000 (the Engle-Granger grid); E experiments of R replications per
# T, all three cases on the same paths. Stored per (nvars, stat, case, T,
# prob): mean and standard error across experiments of the simulated quantile.
#
# Usage, from data-raw/:
#   Rscript po_simulate.R            # nvars = 2..6
# Output: po_quantiles.csv.gz. Takes about ten minutes on two cores.
# PO_E and PO_R (environment variables) override E and R, for a quick test.
Rcpp::sourceCpp("po_kernel.cpp")
RNGkind("L'Ecuyer-CMRG")

args <- commandArgs(trailingOnly = TRUE)
nvars_all <- if (length(args)) as.integer(args) else 2:6
Tgrid <- c(20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000)
probs <- c(seq(0.001, 0.009, 0.001), seq(0.01, 0.99, 0.01), seq(0.991, 0.999, 0.001))
cases <- c("n", "c", "ct")
stats <- c("Zt", "Za")
E <- as.integer(Sys.getenv("PO_E", "40"))
R <- as.integer(Sys.getenv("PO_R", "10000"))
cores <- max(1L, parallel::detectCores())

seeds <- list()
for (nv in nvars_all) {
  set.seed(20261128L + nv)
  s <- .Random.seed
  for (T in Tgrid) { seeds[[paste(nv, T)]] <- s; s <- parallel::nextRNGStream(s) }
}

one <- function(nvars, T) {
  assign(".Random.seed", seeds[[paste(nvars, T)]], envir = .GlobalEnv)
  q <- array(NA_real_, c(E, 2, 3, length(probs)))
  for (e in seq_len(E)) {
    st <- po_null_stats(as.integer(T), R, as.integer(nvars))
    for (a in 1:2) for (s in 1:3)
      q[e, a, s, ] <- quantile(st[, a, s], probs, names = FALSE, type = 8)
  }
  do.call(rbind, lapply(1:2, function(a) do.call(rbind, lapply(1:3, function(s) {
    qq <- matrix(q[, a, s, ], nrow = E)
    data.frame(nvars = nvars, stat = stats[a], case = cases[s], T = T, k = 0L, prob = probs,
               q = colMeans(qq), se = apply(qq, 2, sd) / sqrt(E))
  }))))
}

jobs <- expand.grid(T = Tgrid, nvars = nvars_all)
jobs <- jobs[order(-jobs$T), ]
cat(sprintf("Phillips-Ouliaris simulation: nvars = %s; E = %d, R = %d; %d cores; %d jobs\n",
            paste(nvars_all, collapse = ", "), E, R, cores, nrow(jobs)))
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(jobs)), function(i) one(jobs$nvars[i], jobs$T[i]),
                          mc.cores = cores, mc.preschedule = FALSE)
if (any(vapply(res, inherits, logical(1), "try-error"))) stop("some jobs failed")
out <- do.call(rbind, res)
out <- out[order(out$nvars, out$stat, out$case, out$T, out$prob), ]
write.csv(out, gzfile("po_quantiles.csv.gz"), row.names = FALSE)
cat(sprintf("Done in %s: %d rows -> po_quantiles.csv.gz\n", format(round(Sys.time() - t0, 1)), nrow(out)))
