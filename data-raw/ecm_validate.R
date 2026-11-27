# Size check of the ECM-test response surfaces on fresh simulations, with
# seeds independent of ecm_simulate.R. Run from data-raw/ after
# ecm_response_surfaces.R:  Rscript ecm_validate.R
# For each nvars, case, sample size T (series length) and lag order k, the
# empirical rejection rates at 1/5/10% using (a) the critical values and
# (b) p <= level. The sample sizes lie between the simulated grid points;
# column k is the number p of lagged differences.
Rcpp::sourceCpp("ecm_kernel.cpp")
source("../R/calibration.R")
ecm_rs <- readRDS("ecm_rs.rds")
cases <- c("c", "ct")
R <- as.integer(Sys.getenv("ECM_R", "20000"))

size_table <- function(nvars, T, ks, seed) {
  set.seed(seed)
  st <- ecm_null_stats(as.integer(T), as.integer(ks), R, as.integer(nvars))
  out <- list()
  for (j in seq_along(ks)) for (cs in seq_along(cases)) {
    rs <- ecm_rs[[as.character(nvars)]][[cases[cs]]]
    s <- st[, j, cs]; s <- s[is.finite(s)]
    cv <- unname(rs_critical_values(rs, T, ks[j]))
    p <- vapply(s, function(x) rs_p_value(rs, x, T, ks[j]), numeric(1))
    out[[length(out) + 1]] <- data.frame(
      nvars = nvars, case = cases[cs], T = T, k = ks[j],
      cv01 = mean(s < cv[1]), cv05 = mean(s < cv[2]), cv10 = mean(s < cv[3]),
      p01 = mean(p <= 0.01), p05 = mean(p <= 0.05), p10 = mean(p <= 0.10))
  }
  do.call(rbind, out)
}

cells <- list(list(T = 33L, ks = c(0L, 1L)), list(T = 60L, ks = c(0L, 2L, 4L)),
              list(T = 137L, ks = c(0L, 2L, 4L)), list(T = 290L, ks = c(0L, 4L)),
              list(T = 610L, ks = c(0L, 3L)))
grid <- expand.grid(cell = seq_along(cells), nvars = as.integer(names(ecm_rs)))
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), function(i) {
  cl <- cells[[grid$cell[i]]]
  size_table(grid$nvars[i], cl$T, cl$ks, 6161L + i)
}, mc.cores = parallel::detectCores(), mc.preschedule = FALSE))
cat(sprintf("Empirical size, R = %d replications per cell (MC se at 5%%: %.4f)\n",
            R, sqrt(.05 * .95 / R)))
print(format(res, digits = 3), row.names = FALSE)
cat(sprintf("\nLargest |size - nominal| at 5%%: critical values %.4f, p-values %.4f (%d cells)\n",
            max(abs(res$cv05 - 0.05)), max(abs(res$p05 - 0.05)), nrow(res)))
