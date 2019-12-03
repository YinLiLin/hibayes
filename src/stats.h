#ifndef STATS_H
#define STATS_H

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]

using namespace Rcpp;

double uniform_sample(const double start = 0, const double end = 1);
double norm_sample(const double mean = 0, const double sd = 1);
double gamma_sample(const double shape, const double scale);
double invgamma_sample(const double shape, const double scale);
double chisq_sample(const double shape);
double invchisq_sample(const double shape, const double scale);
double beta_sample(const double a, const double b);
double t_sample(const double shape);
double cauchy_sample(const double location, const double scale);
double exponential_sample(const double scale);
double laplace_sample(const double mean, const double scale);
double rinvgaussian_sample(double mu, double lambda);
NumericVector rdirichlet_sample(const double n, const NumericVector x);
IntegerVector which_c(const NumericVector x, const double value, const int c = 1);
#endif
