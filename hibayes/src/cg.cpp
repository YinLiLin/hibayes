#if !defined(ARMA_64BIT_WORD)
#define ARMA_64BIT_WORD 1
#endif

#include <RcppArmadillo.h>
#include <iostream>
#include <R_ext/Print.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

arma::vec solver_ma(
	const arma::mat A, 
	const arma::vec b, 
	const Nullable<NumericVector> x0 = R_NilValue, 
	const Nullable<NumericVector> lambda = R_NilValue,  
	const double esp = 1e-6, 
	const int outfreq = 100, 
	const bool verbose = true
){
	int m = b.n_elem;
	NumericVector x0_(m);
	arma::vec ap(m);
	if(x0.isNotNull()){
		x0_ = Rcpp::as<NumericVector>(x0);
	}else{
		x0_.fill(0);
	}
	arma::vec x = Rcpp::as<arma::vec>(x0_);
	arma::vec r = b - A * x;
	bool adjust = false;
	arma::vec lambda_;
	if(lambda.isNotNull()){
		lambda_ = Rcpp::as<arma::vec>(lambda);
		adjust = true;
	}
	if(adjust){
		r -= (x % lambda_);
	}
	arma::vec p = r;
	double r2 = sum(square(r));
	double r2update = 0.0;

	double alpha, err, beta;
	for(int i = 0; i < m; i++){
		ap = A * p;
		if(adjust){
			ap += (p % lambda_);
		}
		alpha = r2 / sum(p % ap);
		x += (alpha * p);
		r -= (alpha * ap);
		r2update = sum(square(r));
		err = sqrt(r2update);
		if(verbose && ((i + 1) % outfreq == 0)){
			Rcpp::Rcout.precision(6);
			Rcpp::Rcout << "Iter No." << i << ", err = " << std::fixed << err << std::endl;
		}
		if(err < esp){
			break;
		}
		beta = r2update / r2;
		p = r + beta * p;
		r2 = r2update;
	}
	if(err < esp){
		if(verbose){Rcout << "Convergence: YES" << endl;}
	}else{
		if(verbose){Rcout << "Convergence: NO[try to adjust lambda]" << endl;}
	}
	return x;
}

arma::vec solver_sp(
	const arma::sp_mat A, 
	const arma::vec b, 
	const Nullable<NumericVector> x0 = R_NilValue, 
	const Nullable<NumericVector> lambda = R_NilValue, 
	const double esp = 1e-6, 
	const int outfreq = 100, 
	const bool verbose = true
){
	int m = b.n_elem;
	NumericVector x0_(m);
	arma::vec ap(m);
	if(x0.isNotNull()){
		x0_ = Rcpp::as<NumericVector>(x0);
	}else{
		x0_.fill(0);
	}
	arma::vec x = Rcpp::as<arma::vec>(x0_);
	arma::vec r = b - A * x;
	bool adjust = false;
	arma::vec lambda_;
	if(lambda.isNotNull()){
		lambda_ = Rcpp::as<arma::vec>(lambda);
		adjust = true;
	}
	if(adjust){
		r -= (x % lambda_);
	}
	arma::vec p = r;
	double r2 = sum(square(r));
	double r2update = 0.0;

	double alpha, err, beta;
	for(int i = 0; i < m; i++){
		ap = A * p;
		if(adjust){
			ap += (p % lambda_);
		}
		alpha = r2 / sum(p % ap);
		x += (alpha * p);
		r -= (alpha * ap);
		r2update = sum(square(r));
		err = sqrt(r2update);
		if(verbose && ((i + 1) % outfreq == 0)){
			Rcpp::Rcout.precision(6);
			Rcpp::Rcout << "Iter No." << i << ", err = " << std::fixed << err << std::endl;
		}
		if(err < esp){
			break;
		}
		beta = r2update / r2;
		p = r + beta * p;
		r2 = r2update;
	}
	if(err < esp){
		if(verbose){Rcout << "Convergence: YES" << endl;}
	}else{
		if(verbose){Rcout << "Convergence: NO[try to adjust lambda]" << endl;}
	}
	return x;
}

// [[Rcpp::export]]
Rcpp::List conjgt(
	const NumericMatrix sumstat, 
	const SEXP ldm,
	const Nullable<NumericVector> lambda = R_NilValue,
	const double esp = 1e-6,
	const int outfreq = 100,
	const bool verbose = true
){

	bool sparse;
	if(Rf_inherits(ldm, "dgCMatrix")){
		sparse = TRUE;	
	}else if(Rf_isMatrix(ldm)){
		sparse = FALSE;
	}else{
		throw Rcpp::exception("Unknown type of ldm.");
	}

	arma::vec vx(sumstat.nrow());
	for(int i = 0; i < sumstat.nrow(); i++){
		vx[i] = sqrt(2 * sumstat(i, 0) * (1 - sumstat(i, 0)) * sumstat(i, 3));
	}

	int m;
	int n = mean(na_omit(sumstat(_, 3)));
	arma::sp_mat tXX_sp;
	arma::mat tXX_ma;
	sp_mat::const_iterator start, end;
	if(sparse){
		tXX_sp = as<arma::sp_mat>(ldm);
		if(sumstat.nrow() != tXX_sp.n_rows){
			throw Rcpp::exception("Number of SNPs not equals.");
		}
		m = tXX_sp.n_rows;
		for(int i = 0; i < m; i++){
			tXX_sp(i, i) = vx[i] * vx[i];
			start = tXX_sp.begin_col(i);
	    	end = tXX_sp.end_col(i);
			for (; start != end; ++start){
				int j = start.row();
				if(j > i){
					tXX_sp(i, j) = tXX_sp(j, i) = (*start) * vx[i] * vx[j];
				}
	   		}
		}
	}else{
		tXX_ma = as<arma::mat>(ldm);
		sparse = FALSE;
		if(sumstat.nrow() != tXX_ma.n_rows){
			throw Rcpp::exception("Number of SNPs not equals.");
		}
		m = tXX_ma.n_rows;
		for(int i = 0; i < m; i++){
			tXX_ma(i, i) = vx[i] * vx[i];
			for(int j = (i + 1); j < m; j++){
				tXX_ma(i, j) = tXX_ma(j, i) = tXX_ma(j, i) * vx[i] * vx[j];
			}
		}
	}
	arma::vec xpx(m);
	if(sparse){
		for(int i = 0; i < m; i++){
			xpx[i] = tXX_sp(i, i);
		}
	}else{
		xpx = tXX_ma.diag();
	}
	arma::vec xy(m);
	arma::vec yyi(m); yyi.zeros();
	int count_y = 0;
	for(int k = 0; k < m; k++){
		xy[k] = xpx[k] * sumstat(k, 1);
		if(Rcpp::NumericVector::is_na(sumstat(k, 2))){
			// nothing
		}else{
			yyi[k] = xpx[k] * (sumstat(k, 1) * sumstat(k, 1) + (sumstat(k, 3) - 2) * sumstat(k, 2) * sumstat(k, 2));
			count_y++;
		}
	}
	if(count_y == 0){
		throw Rcpp::exception("Lack of SE.");
	}
	double yy = sum(yyi) / (count_y);
	if(verbose){
		Rcpp::Rcout.precision(4);
		Rcpp::Rcout << "Prior parameters:" << std::endl;
		Rcpp::Rcout << "    Model fitted at [Conjugate Gradient]" << std::endl;
		Rcpp::Rcout << "    Maximum iteration number: " << m << std::endl;
		Rcpp::Rcout << "    Phenotypic var " << std::fixed << yy / (n - 1) << std::endl;
	}
	arma::vec g;
	double vara, vare;
	if(sparse){
		g = solver_sp(tXX_sp, xy, R_NilValue, lambda, esp, outfreq, verbose);
		vara = arma::as_scalar(g.t() * tXX_sp * g) / (n - 1);
		vare = (yy / (n - 1)) - vara;
	}else{
		g = solver_ma(tXX_ma, xy, R_NilValue, lambda, esp, outfreq, verbose);
		vara = arma::as_scalar(g.t() * tXX_ma * g) / (n - 1);
		vare = (yy / (n - 1)) - vara;
	}

	if(verbose){
		Rcpp::Rcout.precision(4);
		Rcpp::Rcout << "Prior parameters:" << std::endl;
		Rcpp::Rcout << "    Genetic var " << std::fixed << vara << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare << std::endl;
	}
	return List::create(
			Named("vara") = vara, 
			Named("vare") = vare,
			Named("g") = g
    );
}
