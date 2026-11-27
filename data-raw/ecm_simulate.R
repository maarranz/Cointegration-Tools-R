# Finite-sample distribution of the ECM test statistic (Kremers, Ericsson and
# Dolado, 1992; Ericsson and MacKinnon, 2002) under the null of no
# cointegration, with the kernel ecm_kernel.cpp.
#
# DGP: n = nvars independent driftless Gaussian random walks y, x_1..x_{n-1}.
# Regression over the maximal sample t = p+2..T:
#   dy_t on [deterministics, y_{t-1}, x_{t-1}, dx_t, dy_{t-1..t-p}, dx_{t-1..t-p}]
# statistic: plain-OLS t ratio of y_{t-1}. Cases "c" (constant) and "ct"
# (constant and trend). p = 0..4 lagged differences of all variables.
#
# Design as eg_simulate.R: T = 20..1000; E = 40 experiments of R = 10,000
# replications per (nvars, T); quantiles at 117 probabilities, stored with
# their standard error across experiments. Column `k` holds p, so that
# fit_surface() can be used unchanged.
#
# Usage, from data-raw/:
#   Rscript ecm_simulate.R            # nvars = 2..6
#   Rscript ecm_simulate.R 4 5        # only these values of nvars
# Output: ecm_quantiles_n<nvars>.csv.gz; partial results in ecm_parts/ (not
# committed), so an interrupted run resumes with the same command.
# ECM_E and ECM_R (environment variables) override E and R, for a quick test.
Rcpp::sourceCpp("ecm_kernel.cpp")
RNGkind("L'Ecuyer-CMRG")

args <- commandArgs(trailingOnly = TRUE)
nvars_all <- if (length(args)) as.integer(args) else 2:6
Tgrid <- c(20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000)
probs <- c(seq(0.001, 0.009, 0.001), seq(0.01, 0.99, 0.01), seq(0.991, 0.999, 0.001))
cases <- c("c", "ct")
E <- as.integer(Sys.getenv("ECM_E", "40"))
R <- as.integer(Sys.getenv("ECM_R", "10000"))
cores <- max(1L, parallel::detectCores())

# One job per (nvars, T), each with its own L'Ecuyer stream (one seed per
# nvars, streams over T). Jobs run largest first across all cores; each
# finished job is saved in ecm_parts/, and a rerun skips the jobs already
# there, so an interrupted run can be resumed with the same command.
dir.create("ecm_parts", showWarnings = FALSE)
part_file <- function(nvars, T) file.path("ecm_parts", sprintf("n%d_T%04d.rds", nvars, T))

lags_for <- function(T) 0:4                     # p, the same for every T

seeds <- list()
for (nv in nvars_all) {
  set.seed(20261127L + nv)
  s <- .Random.seed
  for (T in Tgrid) { seeds[[paste(nv, T)]] <- s; s <- parallel::nextRNGStream(s) }
}

one <- function(nvars, T) {
  assign(".Random.seed", seeds[[paste(nvars, T)]], envir = .GlobalEnv)
  lags <- lags_for(T)
  q <- array(NA_real_, c(E, length(lags), length(cases), length(probs)))
  for (e in seq_len(E)) {
    st <- ecm_null_stats(as.integer(T), as.integer(lags), R, as.integer(nvars))
    for (j in seq_along(lags)) for (s in seq_along(cases)) {
      x <- st[, j, s]
      if (sum(is.finite(x)) > 0.99 * R)            # cells without residual df stay NA
        q[e, j, s, ] <- quantile(x, probs, names = FALSE, type = 8, na.rm = TRUE)
    }
  }
  out <- do.call(rbind, lapply(seq_along(cases), function(s) do.call(rbind, lapply(seq_along(lags), function(j) {
    qq <- matrix(q[, j, s, ], nrow = E)
    data.frame(nvars = nvars, case = cases[s], T = T, k = lags[j], prob = probs,
               q = colMeans(qq), se = apply(qq, 2, sd) / sqrt(E))
  }))))
  saveRDS(out, part_file(nvars, T))
  TRUE
}

jobs <- expand.grid(T = Tgrid, nvars = nvars_all)
jobs$cost <- jobs$T * vapply(jobs$T, function(T) length(lags_for(T)), numeric(1)) * (1 + jobs$nvars / 10)
jobs <- jobs[order(-jobs$cost), ]
todo <- jobs[!file.exists(part_file(jobs$nvars, jobs$T)), ]
cat(sprintf("ECM-test simulation: nvars = %s; E = %d, R = %d; %d cores; %d of %d jobs to run\n",
            paste(nvars_all, collapse = ", "), E, R, cores, nrow(todo), nrow(jobs)))
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(todo)), function(i) one(todo$nvars[i], todo$T[i]),
                          mc.cores = cores, mc.preschedule = FALSE)
failed <- !vapply(res, isTRUE, logical(1))
if (any(failed)) {
  stop(sum(failed), " job(s) failed, e.g. nvars = ", todo$nvars[failed][1], ", T = ",
       todo$T[failed][1], ": ", format(res[failed][[1]]), "\nRerun the same command to retry them.")
}
cat(sprintf("Simulations done in %s\n", format(round(Sys.time() - t0, 1))))

for (nv in nvars_all) {
  out <- do.call(rbind, lapply(Tgrid, function(T) readRDS(part_file(nv, T))))
  out <- out[order(out$case, out$T, out$k, out$prob), ]
  f <- sprintf("ecm_quantiles_n%d.csv.gz", nv)
  write.csv(out, gzfile(f), row.names = FALSE)
  cat(sprintf("nvars = %d: %d rows -> %s\n", nv, nrow(out), f))
}
