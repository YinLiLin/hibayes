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
void Gibbs(arma::mat &A, arma::vec &x, arma::vec &b, double ve);
void Gibbs(arma::sp_mat &A, arma::vec &x, arma::vec &b, double ve);
void Gibbs(arma::sp_mat &A, arma::vec &x, arma::vec &b, double ve, int iter);
void solve(arma::mat &A, double lambda = 0.0);
void eigen_sym_dc(arma::mat &A, arma::vec &eigval);
#endif
