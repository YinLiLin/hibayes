#include "solver.h"

vec PCGv(mat A, vec b, int maxiter, const double tol){
	vec dA = A.diag();
	// stable checking, using available func to speed up
	for(int i = 0; i < dA.n_elem; i++){
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
	int iter = 0;
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

mat PCGm(mat A, mat B, int maxiter, const double tol){
	
	int n_iter = B.n_cols;
	mat x = zeros<mat>(A.n_rows, n_iter);
	for (int i = 0; i < n_iter; i++){
		x.col(i) = PCGv(A, B.col(i), maxiter, tol);
	}// end for loop
	return(x);
}// end function

template<typename T>
vec CG(
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
template vec CG(const sp_mat A, const arma::vec b, const Nullable<NumericVector> x0, const Nullable<NumericVector> lambda,  const double esp, const int outfreq, const bool verbose);
template vec CG(const mat A, const arma::vec b, const Nullable<NumericVector> x0, const Nullable<NumericVector> lambda,  const double esp, const int outfreq, const bool verbose);

template<typename T>
void Gibbs(T &A, vec &x, vec &b, double ve){
    int n = b.n_elem;
    for(int i = 0; i < n; i++){
        double aii = A(i, i);
        if(aii){
            double invlhs  = 1.0 / aii;
            double u = invlhs * (b[i] - dot(A.col(i), x)) + x[i];
            x[i] = norm_sample(u, sqrt(invlhs * ve));
        }
    }
}
template void Gibbs(sp_mat &A, vec &x, vec &b, double ve);
template void Gibbs(mat &A, vec &x, vec &b, double ve);

bool solver_chol(
	arma::mat &A, 
	double lambda = 0.0
)
{
	int n = A.n_cols;
    int info = 0, int_n = (int) n;
    char uplo = 'L';
	vec A_diag = diagvec(A);
	if(lambda){
		A.diag() += lambda;
	}
	double * Aiptr = A.memptr();

	dpotrf_(&uplo, &int_n, Aiptr, &int_n, &info);

    if(info){
		for (int j = 0; j < n; j++) {
			for (int i = (j + 1); i < n; i++) {
				A(i, j) = A(j, i);
			}
		}
		A.diag() = A_diag;
		return false;
    } else {
		double deti = 0.0;
		double toler = abs(Aiptr[0]) * sqrt(datum::eps);
        for (int i = 0; i < n; i++) {
            double d_buf = Aiptr[i * n + i];
			deti = d_buf * d_buf;
			if(deti < toler){
				for(int j = 0; j <= i; j++){
					Aiptr[i * n + j] = 0.0;
				}
			}
        }
		dpotri_(&uplo, &int_n, Aiptr, &int_n, &info);
        if(info){
			for (int j = 0; j < n; j++) {
				for (int i = (j + 1); i < n; i++) {
					A(i, j) = A(j, i);
				}
			}
			A.diag() = A_diag;
			return false;
        }else{
			for (int j = 0; j < n; j++) {
                for (int i = (j + 1); i < n; i++) {
					A(j, i) = A(i, j);
				}
            }
        }
    }
    return true;
}

void solver_lu(
	arma::mat &A,
	double lambda = 0.0
)
{
	int n = A.n_cols;
	int *IPIV = new int[n];
    int LWORK;
    double wkopt;
    double *WORK = new double[4 * n];
    int INFO;
	if(lambda){
		A.diag() += lambda;
	}

	double * Aiptr = A.memptr();
	dgetrf_(&n, &n, Aiptr, &n, IPIV, &INFO);

    if(INFO){
		throw Rcpp::exception("matrix is not invertible, try to specify parameter 'lambda' with a small value, eg: 0.001 or bigger");
    }else{
		double rcond;
		double anorm = dlange_("1", &n, &n, Aiptr, &n, WORK);
		int *IPIVn = new int[n];
		dgecon_("1", &n, Aiptr, &n, &anorm, &rcond, WORK, IPIVn, &INFO);
		delete[] IPIVn;
		if(rcond <= datum::eps){
			throw Rcpp::exception("matrix is not invertible, try to specify parameter 'lambda' with a small value, eg: 0.001 or bigger");
		}else{
			LWORK = -1;
			dgetri_(&n, Aiptr, &n, IPIV, &wkopt, &LWORK, &INFO);
			LWORK = (int)wkopt;
			WORK = (double*)malloc(LWORK * sizeof(double));
			dgetri_(&n, Aiptr, &n, IPIV, WORK, &LWORK, &INFO);
			if(INFO){
				throw Rcpp::exception("U matrix of LU decomposition is singular.");
			}
		}
	}
    delete[] IPIV;
    delete[] WORK;
}

void solve(
	arma::mat &A,
	double lambda
)
{
	if(!solver_chol(A, lambda)){
		solver_lu(A, lambda);
	}
}
