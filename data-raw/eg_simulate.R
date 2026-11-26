# Finite-sample distribution of the Engle-Granger residual-based ADF statistic
# under the null of no cointegration, with the kernel tau_kernel.cpp (copied
# from unitroottools/data-raw; keep the two copies identical).
#
# DGP: n = nvars independent driftless Gaussian random walks y_1..y_n. Step 1:
# regression of y_1 on [deterministic terms, y_2..y_n] over t = 1..T, cases
# "n" (none), "c" (constant), "ct" (constant and trend). Step 2: ADF t ratio of
# the residuals u_t (regression of du_t on u_{t-1}, du_{t-1..t-k}, without
# deterministic terms) over the maximal sample t = k+2..T, as
# unitroottools::adf(trend = "n").
#
# Design as tau_simulate.R (unitroottools, D21): T = 20..1000;
# k = 0..min(12, Schwert's kmax), plus k = 13..21 for T >= 150; E experiments
# of R replications per T, all k and all three cases on the same paths.
# Stored per (nvars, case, T, k, prob): mean and standard error across
# experiments of the simulated quantile.
#
# Usage, from data-raw/:
#   Rscript eg_simulate.R            # nvars = 2..6
#   Rscript eg_simulate.R 4 5        # only these values of nvars
# Output: one file per nvars, eg_quantiles_n<nvars>.csv.gz. Partial results
# are kept in eg_parts/ (not committed); see below.
# EG_E and EG_R (environment variables) override E and R, for a quick test.
Rcpp::sourceCpp("tau_kernel.cpp")
RNGkind("L'Ecuyer-CMRG")

args <- commandArgs(trailingOnly = TRUE)
nvars_all <- if (length(args)) as.integer(args) else 2:6
Tgrid <- c(20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000)
probs <- c(seq(0.001, 0.009, 0.001), seq(0.01, 0.99, 0.01), seq(0.991, 0.999, 0.001))
cases <- c("n", "c", "ct")
E <- as.integer(Sys.getenv("EG_E", "40"))
R <- as.integer(Sys.getenv("EG_R", "10000"))
cores <- max(1L, parallel::detectCores())

# One job per (nvars, T), each with its own L'Ecuyer stream (one seed per
# nvars, streams over T). Jobs run largest first across all cores; each
# finished job is saved in eg_parts/, and a rerun skips the jobs already
# there, so an interrupted run can be resumed with the same command.
dir.create("eg_parts", showWarnings = FALSE)
part_file <- function(nvars, T) file.path("eg_parts", sprintf("n%d_T%04d.rds", nvars, T))

lags_for <- function(T) {
  lags <- 0:min(12L, as.integer(floor(12 * (T / 100)^0.25)))
  if (T >= 150) lags <- c(lags, 13:21)
  lags
}

seeds <- list()
for (nv in nvars_all) {
  set.seed(20261126L + nv)
  s <- .Random.seed
  for (T in Tgrid) { seeds[[paste(nv, T)]] <- s; s <- parallel::nextRNGStream(s) }
}

one <- function(nvars, T) {
  assign(".Random.seed", seeds[[paste(nvars, T)]], envir = .GlobalEnv)
  lags <- lags_for(T)
  q <- array(NA_real_, c(E, length(lags), length(cases), length(probs)))
  for (e in seq_len(E)) {
    st <- tau_null_stats(as.integer(T), as.integer(lags), R, as.integer(nvars))
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
cat(sprintf("Engle-Granger simulation: nvars = %s; E = %d, R = %d; %d cores; %d of %d jobs to run\n",
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
  f <- sprintf("eg_quantiles_n%d.csv.gz", nv)
  write.csv(out, gzfile(f), row.names = FALSE)
  cat(sprintf("nvars = %d: %d rows -> %s\n", nv, nrow(out), f))
}
