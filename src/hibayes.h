#ifndef HIBAYES_
#define HIBAYES_

#if !defined(ARMA_64BIT_WORD)
#define ARMA_64BIT_WORD 1
#endif
#define ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS 1
#include <RcppArmadillo.h>
#include <iostream>
#include <R_ext/Print.h>
#include <R.h>
#include <Rmath.h>
#define R_NO_REMAP
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include "MyTimer.h"
#include "stats.h"
#include "omp_set.h"

#ifdef ARMA_USE_LAPACK
#if !defined(ARMA_BLAS_CAPITALS)
#define arma_daxpy daxpy
#else
#define arma_daxpy DAXPY
#endif
extern "C"
void arma_fortran(arma_daxpy)(int* n, double* da, double* dx, int* incx, double* dy, int* incy);
#endif

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

#endif
