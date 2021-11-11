#include "solver.h"

vec PCGv(mat A, vec b, size_t maxiter, const double tol){
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

mat PCGm(mat A, mat B, size_t maxiter, const double tol){
	
	size_t n_iter = B.n_cols;
	mat x = zeros<mat>(A.n_rows, n_iter);
	for (size_t i = 0; i < n_iter; i++){
		x.col(i) = PCGv(A, B.col(i), maxiter, tol);
	}// end for loop
	return(x);
}// end function

template<typename T>
vec solver(
	const T A, 
	const arma::vec b, 
	const Nullable<NumericVector> x0, 
	const Nullable<NumericVector> lambda,  
	const double esp, 
	const int outfreq, 
	const bool verbose
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
template vec solver(const sp_mat A, const arma::vec b, const Nullable<NumericVector> x0, const Nullable<NumericVector> lambda,  const double esp, const int outfreq, const bool verbose);
template vec solver(const mat A, const arma::vec b, const Nullable<NumericVector> x0, const Nullable<NumericVector> lambda,  const double esp, const int outfreq, const bool verbose);

template<typename T>
void Gibbs(T &A, vec &x, vec &b){
    int n = b.n_elem;
    for(int i = 0; i < n; i++){
        double aii = A(i,i);
        if(aii){
            double invlhs  = 1.0 / aii;
            double μ = invlhs * (b[i] - as_scalar(A.col(i).t() * x)) + x[i];
            x[i] = norm_sample()*sqrt(invlhs) + μ;
        }
    }
}
template void Gibbs(sp_mat &A, vec &x, vec &b);
template void Gibbs(mat &A, vec &x, vec &b);
