// Null distribution of the ECM (error-correction) test statistic of
// Kremers, Ericsson and Dolado (1992) / Ericsson and MacKinnon (2002), by
// sample size T and number p of lagged differences, for n = nvars variables.
//
// DGP: n independent driftless Gaussian random walks y, x_1..x_{n-1}
// (starting at zero). Regression, over the maximal sample t = p+2..T:
//   dy_t on [deterministics, y_{t-1}, x_{t-1}, dx_t,
//            dy_{t-1..t-p}, dx_{t-1..t-p}]
// statistic: OLS t ratio of the coefficient on y_{t-1} (plain OLS standard
// error, no HAC). Cases: 0 "c" (constant), 1 "ct" (constant and trend).
//
// Interface (R): ecm_null_stats(T, lags, reps, nvars) returns an array
// reps x length(lags) x 2. All p use the same paths. R's RNG (rnorm).
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// t ratio of column `pos` (0-based within cols) from the Gram system
static double t_ratio(const std::vector<double>& G, const std::vector<double>& g,
                      double yy, int P, int n, const std::vector<int>& cols, int pos) {
  int p = cols.size();
  if (n - p < 1) return NA_REAL;
  std::vector<double> L(p * p, 0.0), z(p), beta(p), w(p, 0.0);
  for (int i = 0; i < p; ++i)
    for (int j = 0; j <= i; ++j) {
      double s = G[cols[i] * P + cols[j]];
      for (int m = 0; m < j; ++m) s -= L[i * p + m] * L[j * p + m];
      if (i == j) {
        if (!(s > 0)) return NA_REAL;
        L[i * p + i] = std::sqrt(s);
      } else {
        L[i * p + j] = s / L[j * p + j];
      }
    }
  for (int i = 0; i < p; ++i) {
    double s = g[cols[i]];
    for (int m = 0; m < i; ++m) s -= L[i * p + m] * z[m];
    z[i] = s / L[i * p + i];
  }
  double fit_ss = 0.0;
  for (int i = 0; i < p; ++i) fit_ss += z[i] * z[i];
  double rss = yy - fit_ss;
  for (int i = p - 1; i >= 0; --i) {
    double s = z[i];
    for (int m = i + 1; m < p; ++m) s -= L[m * p + i] * beta[m];
    beta[i] = s / L[i * p + i];
  }
  for (int i = 0; i < p; ++i) {
    double s = (i == pos) ? 1.0 : 0.0;
    for (int m = 0; m < i; ++m) s -= L[i * p + m] * w[m];
    w[i] = s / L[i * p + i];
  }
  double v = 0.0;
  for (int i = 0; i < p; ++i) v += w[i] * w[i];
  return beta[pos] / std::sqrt(rss / (n - p) * v);
}

// [[Rcpp::export]]
NumericVector ecm_null_stats(int T, IntegerVector lags, int reps, int nvars) {
  int nl = lags.size(), k = nvars - 1;
  NumericVector out(reps * nl * 2);
  out.attr("dim") = IntegerVector::create(reps, nl, 2);
  std::vector<std::vector<double>> W(nvars, std::vector<double>(T));
  for (int r = 0; r < reps; ++r) {
    for (int v = 0; v < nvars; ++v) {
      NumericVector e = rnorm(T);
      double s = 0.0;
      for (int t = 0; t < T; ++t) { s += e[t]; W[v][t] = s; }
    }
    for (int li = 0; li < nl; ++li) {
      int p = lags[li];
      // columns: 0 const, 1 trend, 2 y_{t-1}, 3..2+k x_{t-1}, then k dx_t,
      // then for j = 1..p: dy_{t-j}, k dx_{t-j}
      int P = 3 + 2 * k + p * nvars;
      std::vector<double> G(P * P, 0.0), g(P, 0.0), x(P);
      double yy = 0.0;
      int n = 0;
      for (int t = p + 1; t < T; ++t) {          // 0-based t: needs t - p - 1 >= 0
        int c = 0;
        x[c++] = 1.0;
        x[c++] = t + 1;
        x[c++] = W[0][t - 1];
        for (int v = 1; v < nvars; ++v) x[c++] = W[v][t - 1];
        for (int v = 1; v < nvars; ++v) x[c++] = W[v][t] - W[v][t - 1];
        for (int j = 1; j <= p; ++j)
          for (int v = 0; v < nvars; ++v) x[c++] = W[v][t - j] - W[v][t - j - 1];
        double dy = W[0][t] - W[0][t - 1];
        for (int a = 0; a < P; ++a) {
          g[a] += x[a] * dy;
          for (int b = 0; b <= a; ++b) G[a * P + b] += x[a] * x[b];
        }
        yy += dy * dy;
        ++n;
      }
      for (int a = 0; a < P; ++a)
        for (int b = a + 1; b < P; ++b) G[a * P + b] = G[b * P + a];
      std::vector<int> cc, cct;
      cc.push_back(0); cct.push_back(0); cct.push_back(1);
      for (int a = 2; a < P; ++a) { cc.push_back(a); cct.push_back(a); }
      out[r + reps * (li + nl * 0)] = t_ratio(G, g, yy, P, n, cc, 1);
      out[r + reps * (li + nl * 1)] = t_ratio(G, g, yy, P, n, cct, 2);
    }
  }
  return out;
}
