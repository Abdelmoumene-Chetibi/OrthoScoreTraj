#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fast_knn_transition_cpp(NumericMatrix pca_coords, int k) {
  int n = pca_coords.nrow();
  int d = pca_coords.ncol();
  if (k >= n) k = n - 1;
  
  NumericMatrix P(n, n);
  std::vector<std::pair<double, int>> dists(n);
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) {
        dists[j] = std::make_pair(1e30, j);
        continue;
      }
      double sum_sq = 0.0;
      for (int c = 0; c < d; ++c) {
        double diff = pca_coords(i, c) - pca_coords(j, c);
        sum_sq += diff * diff;
      }
      dists[j] = std::make_pair(std::sqrt(sum_sq), j);
    }
    
    std::nth_element(dists.begin(), dists.begin() + k, dists.end());
    std::sort(dists.begin(), dists.begin() + k);
    
    double sigma = dists[k - 1].first;
    if (sigma < 1e-8) sigma = 1e-8;
    
    double row_sum = 0.0;
    for (int m = 0; m < k; ++m) {
      int nbr = dists[m].second;
      double w = std::exp(-(dists[m].first * dists[m].first) / (sigma * sigma));
      P(i, nbr) = w;
      row_sum += w;
    }
    if (row_sum > 0.0) {
      for (int m = 0; m < k; ++m) {
        P(i, dists[m].second) /= row_sum;
      }
    }
  }
  return P;
}

// [[Rcpp::export]]
NumericVector fast_geodesic_dist_cpp(NumericMatrix P, int root_idx_1based, double alpha = 0.9, int n_steps = 40) {
  int n = P.nrow();
  int root = root_idx_1based - 1;
  if (root < 0 || root >= n) stop("Invalid root_idx_1based");
  
  // Sparse adjacency extraction for O(k*N) Matrix-Vector diffusion instead of O(N^3)
  std::vector<std::vector<std::pair<int, double>>> adj(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double val = P(i, j);
      if (val > 0.0) {
        adj[i].push_back(std::make_pair(j, val));
      }
    }
  }
  
  std::vector<double> prob(n, 0.0);
  std::vector<double> next_prob(n, 0.0);
  prob[root] = 1.0;
  
  for (int step = 0; step < n_steps; ++step) {
    std::fill(next_prob.begin(), next_prob.end(), 0.0);
    next_prob[root] = (1.0 - alpha);
    for (int i = 0; i < n; ++i) {
      double pi = prob[i];
      if (pi <= 1e-15) continue;
      for (size_t idx = 0; idx < adj[i].size(); ++idx) {
        int j = adj[i][idx].first;
        double w = adj[i][idx].second;
        next_prob[j] += alpha * pi * w;
      }
    }
    prob.swap(next_prob);
  }
  
  NumericVector pt(n);
  double min_d = 1e30, max_d = -1e30;
  for (int i = 0; i < n; ++i) {
    double d = -std::log(prob[i] + 1e-12);
    pt[i] = d;
    if (d < min_d) min_d = d;
    if (d > max_d) max_d = d;
  }
  
  double range = max_d - min_d;
  if (range > 1e-12) {
    for (int i = 0; i < n; ++i) {
      pt[i] = (pt[i] - min_d) / range;
    }
  }
  pt[root] = 0.0;
  return pt;
}