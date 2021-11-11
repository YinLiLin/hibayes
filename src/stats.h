#ifndef STATS_H
#define STATS_H

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]

using namespace Rcpp;

double uniform_sample(double start = 0, double end = 1);
double norm_sample(double mean = 0, double sd = 1);
double gamma_sample(double shape, double scale);
double invgamma_sample(double shape, double scale);
double chisq_sample(double shape);
double invchisq_sample(double shape, double scale);
double beta_sample(double a, double b);
double t_sample(double shape);
double cauchy_sample(double location, double scale);
double exponential_sample(double scale);
double laplace_sample(double mean, double scale);
double rinvgaussian_sample(double mu, double lambda);
NumericVector rdirichlet_sample(double n, NumericVector x);
IntegerVector which_c(NumericVector x, double value, int c = 1);
#endif
