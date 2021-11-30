#include "solver.h"

// [[Rcpp::export]]
Rcpp::List conjgt_spa(
	const NumericMatrix sumstat, 
	const arma::sp_mat ldm,
	const Nullable<NumericVector> lambda = R_NilValue,
	const double esp = 1e-6,
	const int outfreq = 100,
	const bool verbose = true
){
	int m;
	int n = mean(na_omit(sumstat(_, 3)));
	sp_mat::const_iterator start, end;
	if(sumstat.nrow() != ldm.n_rows){
		throw Rcpp::exception("Number of SNPs not equals.");
	}
	m = ldm.n_rows;
	arma::vec xpx(m);
	arma::vec vx(m);
	for(int i = 0; i < m; i++){
		vx[i] = ldm(i, i);
		xpx[i] = vx[i] * n;
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
	double vg, ve;
	g = CG(ldm, xy / n, R_NilValue, lambda, esp, outfreq, verbose);
	vg = n * arma::as_scalar(g.t() * ldm * g) / (n - 1);
	ve = (yy / (n - 1)) - vg;
	if(verbose){
		Rcpp::Rcout.precision(4);
		Rcpp::Rcout << "Prior parameters:" << std::endl;
		Rcpp::Rcout << "    Genetic var " << std::fixed << vg << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << ve << std::endl;
	}
	return List::create(
			Named("vg") = vg, 
			Named("ve") = ve,
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

	int m;
	int n = mean(na_omit(sumstat(_, 3)));
	if(sumstat.nrow() != ldm.n_rows){
		throw Rcpp::exception("Number of SNPs not equals.");
	}
	m = ldm.n_rows;

	arma::vec xpx(m);
	arma::vec vx(m);
	vx = ldm.diag();
	xpx = ldm.diag() * n;
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
	double vg, ve;

	g = CG(ldm, xy / n, R_NilValue, lambda, esp, outfreq, verbose);
	vg = n * arma::as_scalar(g.t() * ldm * g) / (n - 1);
	ve = (yy / (n - 1)) - vg;

	if(verbose){
		Rcpp::Rcout.precision(4);
		Rcpp::Rcout << "Prior parameters:" << std::endl;
		Rcpp::Rcout << "    Genetic var " << std::fixed << vg << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << ve << std::endl;
	}
	return List::create(
			Named("vg") = vg, 
			Named("ve") = ve,
			Named("g") = g
    );
}
