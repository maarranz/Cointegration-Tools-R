# Size of the Phillips-Ouliaris tests with the package's tables, in fresh
# simulations (new seed) through the R code of po_test() (po_stats(),
# rs_p_value()): independent random walks, 5% tests, 4000 replications per
# cell. Usage, from data-raw/ after make_sysdata.R: Rscript po_validate.R
for (f in c("lrv.R", "calibration.R", "po_test.R")) source(file.path("../R", f))
load("../R/sysdata.rda")
RNGkind("L'Ecuyer-CMRG"); set.seed(99)
R <- 4000
cells <- expand.grid(T = c(30, 60, 100, 250, 600), nvars = c(2, 4, 6), case = c("n", "c", "ct"),
                     stringsAsFactors = FALSE)
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(cells)), function(i) {
  T <- cells$T[i]; nv <- cells$nvars[i]; cs <- cells$case[i]
  rej <- replicate(R, {
    Y <- apply(matrix(rnorm(T * nv), T), 2, cumsum)
    X <- cbind(if (cs != "n") 1, if (cs == "ct") seq_len(T), Y[, -1])
    u <- .lm.fit(X, Y[, 1])$residuals
    ps <- po_stats(u)
    c(rs_p_value(po_rs[[as.character(nv)]][[cs]]$Zt, ps$Zt, T, 0) < 0.05,
      rs_p_value(po_rs[[as.character(nv)]][[cs]]$Za, ps$Za, T, 0) < 0.05)
  })
  data.frame(cells[i, ], Zt = mean(rej[1, ]), Za = mean(rej[2, ]))
}, mc.cores = parallel::detectCores()))
print(res, row.names = FALSE)
cat(sprintf("\nMax |size - 0.05|: Zt %.4f, Za %.4f (Monte Carlo s.e. %.4f)\n",
            max(abs(res$Zt - 0.05)), max(abs(res$Za - 0.05)), sqrt(0.05 * 0.95 / R)))
