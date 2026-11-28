// Null distribution of the Phillips-Ouliaris (1990) residual-based
// statistics Z_t and Z_alpha, by sample size T, for nvars = 2..6 independent
// driftless Gaussian random walks and three deterministic cases of the
// cointegrating regression (n, c, ct), as tau_kernel.cpp for Engle-Granger.
//
// Residuals u_t of the regression of y_1 on [deterministics, y_2..y_n] over
// t = 1..T; then the regression u_t = rho u_{t-1} + e_t over t = 2..T
// (N = T - 1 rows, no deterministics), and
//   s2     = sum e^2 / N,            sigma2 = sum e^2 / (N - 1)
//   lambda2 = s2 + 2 sum_j k(j/b) gamma_j,  gamma_j = sum_{t>j} e_t e_{t-j} / N
// with the Bartlett kernel and Andrews' (1991) AR(1) bandwidth, exactly as
// R/lrv.R (long_run_cov(e, "bartlett", "andrews")), and
//   Z_alpha = N (rho - 1) - (lambda2 - s2) / (2 S),   S = sum u_{t-1}^2 / N^2
//   Z_t     = sqrt(s2 / lambda2) t_rho - (lambda2 - s2) / (2 sqrt(lambda2 S))
// where t_rho = (rho - 1) / sqrt(sigma2 / sum u_{t-1}^2). These are the
// formulas of po_stats() in R/po_test.R; keep the two in step.
//
// Interface (R): po_null_stats(T, reps, nvars) returns an array
// reps x 2 (Z_t, Z_alpha) x 3 (cases n, c, ct). The RNG is R's (rnorm).
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Residuals of the regression of y1 on [deterministics (case), y2..yn] over
// t = 1..T (normal equations; small number of regressors).
static std::vector<double> coint_residuals(const std::vector<std::vector<double>>& Y,
                                           int cs) {
  int T = Y[0].size(), nv = Y.size();
  int d = (cs == 0) ? 0 : (cs == 1 ? 1 : 2);
  int p = d + nv - 1;
  std::vector<double> Xc(T * p);
  for (int t = 0; t < T; ++t) {
    int c = 0;
    if (d >= 1) Xc[(c++) * T + t] = 1.0;
    if (d >= 2) Xc[(c++) * T + t] = t + 1;
    for (int v = 1; v < nv; ++v) Xc[(c++) * T + t] = Y[v][t];
  }
  std::vector<double> G(p * p), g(p);
  for (int a = 0; a < p; ++a) {
    for (int b = 0; b <= a; ++b) {
      double s = 0;
      for (int t = 0; t < T; ++t) s += Xc[a * T + t] * Xc[b * T + t];
      G[a * p + b] = s; G[b * p + a] = s;
    }
    double s = 0;
    for (int t = 0; t < T; ++t) s += Xc[a * T + t] * Y[0][t];
    g[a] = s;
  }
  // Cholesky solve
  std::vector<double> L(p * p, 0.0), z(p), beta(p);
  for (int i = 0; i < p; ++i)
    for (int j = 0; j <= i; ++j) {
      double s = G[i * p + j];
      for (int m = 0; m < j; ++m) s -= L[i * p + m] * L[j * p + m];
      if (i == j) L[i * p + i] = std::sqrt(s); else L[i * p + j] = s / L[j * p + j];
    }
  for (int i = 0; i < p; ++i) {
    double s = g[i];
    for (int m = 0; m < i; ++m) s -= L[i * p + m] * z[m];
    z[i] = s / L[i * p + i];
  }
  for (int i = p - 1; i >= 0; --i) {
    double s = z[i];
    for (int m = i + 1; m < p; ++m) s -= L[m * p + i] * beta[m];
    beta[i] = s / L[i * p + i];
  }
  std::vector<double> u(T);
  for (int t = 0; t < T; ++t) {
    double f = 0;
    for (int a = 0; a < p; ++a) f += Xc[a * T + t] * beta[a];
    u[t] = Y[0][t] - f;
  }
  return u;
}


static void po_stats(const std::vector<double>& u, double* zt, double* za) {
  int T = u.size(), N = T - 1;
  double sxy = 0, sxx = 0;
  for (int t = 1; t < T; ++t) { sxy += u[t] * u[t - 1]; sxx += u[t - 1] * u[t - 1]; }
  double rho = sxy / sxx;
  std::vector<double> e(N);
  double see = 0;
  for (int t = 1; t < T; ++t) { e[t - 1] = u[t] - rho * u[t - 1]; see += e[t - 1] * e[t - 1]; }
  double s2 = see / N, sigma2 = see / (N - 1);
  // Andrews (1991) bandwidth, AR(1) approximation (R/lrv.R andrews_bandwidth)
  double a01 = 0, a11 = 0;
  for (int t = 1; t < N; ++t) { a01 += e[t] * e[t - 1]; a11 += e[t - 1] * e[t - 1]; }
  double r = a01 / a11, ss = 0;
  for (int t = 1; t < N; ++t) { double d = e[t] - r * e[t - 1]; ss += d * d; }
  ss /= N;
  double num = 4 * r * r * ss * ss / (std::pow(1 - r, 6) * std::pow(1 + r, 2));
  double den = ss * ss / std::pow(1 - r, 4);
  double b = 1.1447 * std::pow(num / den * N, 1.0 / 3.0);
  if (b > N - 1) b = N - 1;
  double lambda2 = s2;
  for (int j = 1; j < N; ++j) {
    double w = 1.0 - j / b;
    if (w <= 0) break;
    double g = 0;
    for (int t = j; t < N; ++t) g += e[t] * e[t - j];
    lambda2 += 2 * w * g / N;
  }
  double S = sxx / ((double)N * N);
  double trho = (rho - 1) / std::sqrt(sigma2 / sxx);
  *za = N * (rho - 1) - (lambda2 - s2) / (2 * S);
  *zt = std::sqrt(s2 / lambda2) * trho - (lambda2 - s2) / (2 * std::sqrt(lambda2 * S));
}

// [[Rcpp::export]]
NumericVector po_stats_cpp(NumericVector u) {
  std::vector<double> uu(u.begin(), u.end());
  double zt, za;
  po_stats(uu, &zt, &za);
  return NumericVector::create(zt, za);
}

// [[Rcpp::export]]
NumericVector po_null_stats(int T, int reps, int nvars) {
  NumericVector out(reps * 2 * 3);
  out.attr("dim") = IntegerVector::create(reps, 2, 3);
  RNGScope scope;
  std::vector<std::vector<double>> Y(nvars, std::vector<double>(T));
  for (int r = 0; r < reps; ++r) {
    for (int v = 0; v < nvars; ++v) {
      NumericVector e = rnorm(T);
      double acc = 0.0;
      for (int t = 0; t < T; ++t) { acc += e[t]; Y[v][t] = acc; }
    }
    for (int cs = 0; cs < 3; ++cs) {
      std::vector<double> u = coint_residuals(Y, cs);
      double zt, za;
      po_stats(u, &zt, &za);
      out[r + reps * (0 + 2 * cs)] = zt;
      out[r + reps * (1 + 2 * cs)] = za;
    }
  }
  return out;
}
