#ifndef SOLVER_H
#define SOLVER_H

#if !defined(ARMA_64BIT_WORD)
#define ARMA_64BIT_WORD 1
#endif
#define ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS 1

#include <RcppArmadillo.h>
#include "stats.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

template<typename T> vec CG(const T A, const arma::vec b, const Nullable<NumericVector> x0 = R_NilValue, const Nullable<NumericVector> lambda = R_NilValue,  const double esp = 1e-6, const int outfreq = 100, const bool verbose = true);
template<typename T> void Gibbs(T &A, vec &x, vec &b, double ve);
void solve(arma::mat &A, double lambda = 0.0);
#endif
