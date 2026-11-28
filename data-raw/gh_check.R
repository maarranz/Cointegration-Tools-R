# Check of the Gregory-Hansen critical values typed into R/gh_test.R: a
# simulation at T = 400 (independent random walks, trim 0.15) of ADF* with
# k = 0 and of Z_alpha* computed as the coefficient statistic N(rho - 1),
# its limit under i.i.d. errors, for m = 1 and 2 regressors and the four
# models. The 5% quantiles should be close to the tables (within a few
# hundredths for ADF*; finite-sample and Monte Carlo error).
# Usage, from data-raw/: Rscript gh_check.R   (a few minutes on two cores)
source("../R/gh_test.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(20261128)
T <- 400; R <- 3000; trim <- 0.15
breaks <- seq.int(ceiling(trim * T), floor((1 - trim) * T))
tt <- seq_len(T)
one_model <- function(y, X, model) {
  st <- vapply(breaks, function(tb) {
    D <- as.numeric(tt > tb)
    Z <- cbind(1, D)
    if (model %in% c("C/T", "C/S/T")) Z <- cbind(Z, tt)
    if (model == "C/S/T") Z <- cbind(Z, tt * D)
    Z <- cbind(Z, X)
    if (model %in% c("C/S", "C/S/T")) Z <- cbind(Z, X * D)
    u <- .lm.fit(Z, y)$residuals
    u1 <- u[-T]; du <- diff(u)
    b <- sum(du * u1) / sum(u1^2)
    e <- du - b * u1
    c(b / sqrt(sum(e^2) / (T - 2) / sum(u1^2)), (T - 1) * b)
  }, numeric(2))
  apply(st, 1, min)
}
models <- names(gh_tables)
for (m in 1:2) {
  sims <- parallel::mclapply(seq_len(R), function(r) {
    W <- apply(matrix(rnorm(T * (m + 1)), T), 2, cumsum)
    sapply(models, function(md) one_model(W[, 1], W[, -1, drop = FALSE], md))
  }, mc.cores = parallel::detectCores())
  A <- simplify2array(sims)                       # 2 x models x R
  for (j in seq_along(models)) {
    q <- apply(A[, j, ], 1, quantile, c(0.01, 0.05, 0.10))
    cat(sprintf("m=%d %-5s ADF* sim %s | table %s\n", m, models[j],
                paste(sprintf("%6.2f", q[, 1]), collapse = ""),
                paste(sprintf("%6.2f", gh_tables[[models[j]]]$ADF[m, ]), collapse = "")))
    cat(sprintf("m=%d %-5s Za*  sim %s | table %s\n", m, models[j],
                paste(sprintf("%7.2f", q[, 2]), collapse = ""),
                paste(sprintf("%7.2f", gh_tables[[models[j]]]$Za[m, ]), collapse = "")))
  }
}
