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


// solve the equation Ax=b, x is a variables
arma::vec PCGv(mat A, vec b, size_t maxiter, const double tol){
	vec dA = A.diag();
	// stable checking, using available func to speed up
	for(size_t i = 0; i < dA.n_elem; i++){
		if(dA[i] == 0)
			dA[i] = 1e-4;
	}// end for i
	vec Minv = 1.0/dA;
	// initialize
	vec x = zeros<vec>(b.n_elem);
	vec r = zeros<vec>(b.n_elem);
	vec r1 = zeros<vec>(b.n_elem);
	vec z1 = zeros<vec>(b.n_elem);
	r = b;
	vec z = Minv % r;
	vec p = z;
	size_t iter = 0;
	double sumr2 = norm(r, 2);
	// PCG main loop 
	while( sumr2 > tol && iter < maxiter){
		iter += 1;
		// move direction
		vec Ap = A*p;
		// step size
		double a = dot(r,z)/dot(p,Ap);
		// move
		x = x + a * p;
		r1 = r - a * Ap;
		z1 = Minv % r1;
		double bet = dot(z1, r1)/dot(z, r);
		p = z1 + bet * p;
		z = z1;
		r = r1;
		sumr2 = norm(r, 2);
	}// end while loop
	if (iter >= maxiter){
		Rcerr << "ERROR: Matrix is Singular!" << endl;
	}
	return(x);
}// end function

arma::mat PCGm(mat A, mat B, size_t maxiter, const double tol){
	
	size_t n_iter = B.n_cols;
	mat x = zeros<mat>(A.n_rows, n_iter);
	for (size_t i = 0; i < n_iter; i++){
		x.col(i) = PCGv(A, B.col(i), maxiter, tol);
	}// end for loop
	return(x);
}// end function

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
Rcpp::List conjgt_spa(
	const NumericMatrix sumstat, 
	const arma::sp_mat ldm,
	const Nullable<NumericVector> lambda = R_NilValue,
	const double esp = 1e-6,
	const int outfreq = 100,
	const bool verbose = true
){

	// bool sparse;
	// if(Rf_inherits(ldm, "dgCMatrix")){
	// 	sparse = TRUE;	
	// }else if(Rf_isMatrix(ldm)){
	// 	sparse = FALSE;
	// }else{
	// 	throw Rcpp::exception("Unknown type of ldm.");
	// }

	// arma::vec vx(sumstat.nrow());
	// for(int i = 0; i < sumstat.nrow(); i++){
	// 	vx[i] = sqrt(2 * sumstat(i, 0) * (1 - sumstat(i, 0)) * sumstat(i, 3));
	// }

	int m;
	int n = mean(na_omit(sumstat(_, 3)));
	// arma::sp_mat ldm;
	// arma::mat ldm;
	sp_mat::const_iterator start, end;
	// if(sparse){
		// ldm = as<arma::sp_mat>(ldm);
		if(sumstat.nrow() != ldm.n_rows){
			throw Rcpp::exception("Number of SNPs not equals.");
		}
		m = ldm.n_rows;
		// for(int i = 0; i < m; i++){
		// 	ldm(i, i) = vx[i] * vx[i];
		// 	start = ldm.begin_col(i);
	 //    	end = ldm.end_col(i);
		// 	for (; start != end; ++start){
		// 		int j = start.row();
		// 		if(j > i){
		// 			ldm(i, j) = ldm(j, i) = (*start) * vx[i] * vx[j];
		// 		}
	 //   		}
		// }
	// }else{
	// 	ldm = as<arma::mat>(ldm);
	// 	sparse = FALSE;
	// 	if(sumstat.nrow() != ldm.n_rows){
	// 		throw Rcpp::exception("Number of SNPs not equals.");
	// 	}
	// 	m = ldm.n_rows;
	// 	// for(int i = 0; i < m; i++){
	// 	// 	ldm(i, i) = vx[i] * vx[i];
	// 	// 	for(int j = (i + 1); j < m; j++){
	// 	// 		ldm(i, j) = ldm(j, i) = ldm(j, i) * vx[i] * vx[j];
	// 	// 	}
	// 	// }
	// }
	arma::vec xpx(m);
	arma::vec vx(m);
	// if(sparse){
		for(int i = 0; i < m; i++){
			vx[i] = ldm(i, i);
			xpx[i] = vx[i] * n;
		}
	// }else{
	// 	vx = ldm.diag();
	// 	xpx = ldm.diag() * n;
	// }
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
	// if(sparse){
		g = solver_sp(ldm, xy / n, R_NilValue, lambda, esp, outfreq, verbose);
		vara = n * arma::as_scalar(g.t() * ldm * g) / (n - 1);
		vare = (yy / (n - 1)) - vara;
	// }else{
	// 	g = solver_ma(ldm, xy / n, R_NilValue, lambda, esp, outfreq, verbose);
	// 	vara = n * arma::as_scalar(g.t() * ldm * g) / (n - 1);
	// 	vare = (yy / (n - 1)) - vara;
	// }

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

// [[Rcpp::export]]
Rcpp::List conjgt_den(
	const NumericMatrix sumstat, 
	const arma::mat ldm,
	const Nullable<NumericVector> lambda = R_NilValue,
	const double esp = 1e-6,
	const int outfreq = 100,
	const bool verbose = true
){

	// bool sparse;
	// if(Rf_inherits(ldm, "dgCMatrix")){
	// 	sparse = TRUE;	
	// }else if(Rf_isMatrix(ldm)){
	// 	sparse = FALSE;
	// }else{
	// 	throw Rcpp::exception("Unknown type of ldm.");
	// }

	// arma::vec vx(sumstat.nrow());
	// for(int i = 0; i < sumstat.nrow(); i++){
	// 	vx[i] = sqrt(2 * sumstat(i, 0) * (1 - sumstat(i, 0)) * sumstat(i, 3));
	// }

	int m;
	int n = mean(na_omit(sumstat(_, 3)));
	// arma::sp_mat ldm;
	// arma::mat ldm;
	// sp_mat::const_iterator start, end;
	// if(sparse){
	// 	ldm = as<arma::sp_mat>(ldm);
	// 	if(sumstat.nrow() != ldm.n_rows){
	// 		throw Rcpp::exception("Number of SNPs not equals.");
	// 	}
	// 	m = ldm.n_rows;
	// 	// for(int i = 0; i < m; i++){
	// 	// 	ldm(i, i) = vx[i] * vx[i];
	// 	// 	start = ldm.begin_col(i);
	//  //    	end = ldm.end_col(i);
	// 	// 	for (; start != end; ++start){
	// 	// 		int j = start.row();
	// 	// 		if(j > i){
	// 	// 			ldm(i, j) = ldm(j, i) = (*start) * vx[i] * vx[j];
	// 	// 		}
	//  //   		}
	// 	// }
	// }else{
		// ldm = as<arma::mat>(ldm);
		// sparse = FALSE;
		if(sumstat.nrow() != ldm.n_rows){
			throw Rcpp::exception("Number of SNPs not equals.");
		}
		m = ldm.n_rows;
		// for(int i = 0; i < m; i++){
		// 	ldm(i, i) = vx[i] * vx[i];
		// 	for(int j = (i + 1); j < m; j++){
		// 		ldm(i, j) = ldm(j, i) = ldm(j, i) * vx[i] * vx[j];
		// 	}
		// }
	// }
	arma::vec xpx(m);
	arma::vec vx(m);
	// if(sparse){
	// 	for(int i = 0; i < m; i++){
	// 		vx[i] = ldm(i, i);
	// 		xpx[i] = vx[i] * n;
	// 	}
	// }else{
		vx = ldm.diag();
		xpx = ldm.diag() * n;
	// }
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
	// if(sparse){
	// 	g = solver_sp(ldm, xy / n, R_NilValue, lambda, esp, outfreq, verbose);
	// 	vara = n * arma::as_scalar(g.t() * ldm * g) / (n - 1);
	// 	vare = (yy / (n - 1)) - vara;
	// }else{
		g = solver_ma(ldm, xy / n, R_NilValue, lambda, esp, outfreq, verbose);
		vara = n * arma::as_scalar(g.t() * ldm * g) / (n - 1);
		vare = (yy / (n - 1)) - vara;
	// }

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
