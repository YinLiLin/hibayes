// #if !defined(ARMA_64BIT_WORD)
// #define ARMA_64BIT_WORD 1
// #endif

#include <RcppArmadillo.h>
#include <iostream>
#include <R_ext/Print.h>
#include "MyTimer.h"
#include "stats.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

double spsubmatmtp(const arma::sp_mat X, const arma::vec g, const arma::uvec windxi){

	arma::vec gmtpX(windxi.n_elem);
	double summed, varw;
    for(int i = 0; i < windxi.n_elem; i++){
    	summed = 0;
    	for(int j = 0; j < windxi.n_elem; j++){
    		if(g[windxi[j]]){summed += g[windxi[j]] * X(windxi[j], windxi[i]);}
    	}
    	gmtpX[i] = summed;
    }
    varw = sum(gmtpX % g(windxi));
	// Return result
	return varw;
}

// [[Rcpp::export]]
Rcpp::List SBayesRR_spa(
	const NumericMatrix sumstat, 
	const arma::sp_mat & ldm,
	const int niter = 50000,
	const int nburn = 20000,
	const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vara = R_NilValue,
    const Nullable<double> dfvara = R_NilValue,
    const Nullable<double> s2vara = R_NilValue,
    const Nullable<double> vare = R_NilValue,
    const Nullable<double> dfvare = R_NilValue,
    const Nullable<double> s2vare = R_NilValue,
	const int outfreq = 100,
	const bool verbose = true
){
	// arma::vec vx(sumstat.nrow());
	// for(int i = 0; i < sumstat.nrow(); i++){
	// 	vx[i] = sqrt(2 * sumstat(i, 0) * (1 - sumstat(i, 0)) * sumstat(i, 3));
	// }

	int m;
	// arma::sp_mat ldm;
	// arma::mat tXX_ma;
	sp_mat::const_iterator start, end;

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
	int n = mean(na_omit(sumstat(_, 3)));
	bool WPPA = false;
	int count = 0;
	double xx, gi, rhs, lhs, uhat, v, hsq;
	double vara_, vargi, dfvara_, s2vara_, vare_, dfvare_, s2vare_, sumvg, sumvargi;
	arma::vec g(m); g.zeros();
	arma::vec g_store(m); g_store.zeros();
	arma::vec xpx(m); xpx.zeros();
	arma::vec vx(m); vx.zeros();
	for(int i = 0; i < m; i++){
		vx[i] = ldm(i, i);
		xpx[i] = vx[i] * n;
	}
	arma::vec xy(m);
	arma::vec yyi(m); yyi.zeros();
	NumericVector sumvg_store(niter - nburn + 1), vara_store(niter - nburn + 1), vare_store(niter - nburn + 1), hsq_store(niter - nburn + 1);
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
	double vary = yy / n;
	if(dfvara.isNotNull()){
		dfvara_ = as<double>(dfvara);
	}else{
		dfvara_ = 4;
	}
	if(dfvara_ <= 2){
        throw Rcpp::exception("dfvara should not be less than 2.");
    }
	if(vara.isNotNull()){
		vara_ = as<double>(vara);
	}else{
		vara_ = ((dfvara_ - 2) / dfvara_) * yy / (n - 1);
	}
	if(vare.isNotNull()){
		vare_ = as<double>(vare);
	}else{
		vare_ = yy / (n - 1);
	}
	if(dfvare.isNotNull()){
		dfvare_ = as<double>(dfvare);
	}else{
		dfvare_ = -2;
	}
	if(s2vara.isNotNull()){
		s2vara_ = as<double>(s2vara);
	}else{
		s2vara_ = (vara_ / (sum(xpx) / (n - 1))) * (dfvara_ - 2) / dfvara_;
	}
	if(s2vare.isNotNull()){
		s2vare_ = as<double>(s2vare);
	}else{
		s2vare_ = 0;
	}
	arma::vec r_hat = xy;
	vec vargx(m);
    // double vargfilter = yy / (n - 1);
    double varei;
	vec varediff(m);

	for(int i = 0; i < m; i++){
		start = ldm.begin_col(i);
		end = ldm.end_col(i);
		double j = 0;
		for (; start != end; ++start){
			j++;
		}
		varediff[i] = (m - j) / m;
	}

	double varg = (vara_ / (sum(xpx) / (n - 1)));

	if(niter < nburn){
		throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
	}

    // for gwas
    R_xlen_t nw;
    int ni;
    double varw;
    arma::vec windindx_(m);
    Rcpp::List windx(m);
    arma::uvec windxi;
    NumericVector wppai(m);
    NumericVector wgvei(m);
    if(windindx.isNotNull()){
        windindx_ = as<arma::vec>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = wppai[seq(0, (nw - 1))];
        wgvei = wgvei[seq(0, (nw - 1))];
        for(R_xlen_t w = 0; w < nw; w++){
            windx[w] = find(windindx_ == (w+1));
        }
    }

	if(verbose){
		Rcpp::Rcout.precision(4);
		Rcpp::Rcout << "Prior parameters:" << std::endl;
		Rcpp::Rcout << "    Model fitted at [SBayes Ridge Regression]" << std::endl;
		Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
		Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
		Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
		Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << s2vara_ << std::endl;
		Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
		Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << s2vare_ << std::endl;
		Rcpp::Rcout << "    Phenotypic var " << std::fixed << yy / (n - 1) << std::endl;
		if(WPPA){
            Rcpp::Rcout << "    Number of windows " << nw << std::endl;
            Rcpp::Rcout << "    GVE threshold for windows " << wppa << std::endl;
        }
		Rcpp::Rcout << "MCMC started: " << std::endl;
        Rcpp::Rcout << " Iter" << "\t";
        Rcpp::Rcout << "SigmaSq" << "\t";
		Rcpp::Rcout << "GenVar" << "\t";
		Rcpp::Rcout << "ResVar" << "\t";
		Rcpp::Rcout << "hsq" << "\t";
		Rcpp::Rcout << "Timeleft" << std::endl;
	}

	MyTimer timer;
    timer.step("start");
    MyTimer timer_loop;
    timer_loop.step("start");
    double tt0; 
    int tt, hor, min, sec;

	// MCMC procedure
	for(int iter = 0; iter < niter; iter++){

		sumvg = 0;
		sumvargi = 0;

		// loop on snps
		for(int i = 0; i < m; i++){
			xx = xpx[i];
			gi = g[i];

			varei = varediff[i] * vara_ + vare_;

			// right hand
			rhs = r_hat[i];
			if(gi){rhs += xx * gi;}

			// left hand
			lhs = xx / varei;
            v = xx + varei / varg;

            // sample snp effect
            gi = norm_sample(rhs / v, sqrt(varei / v));

            vargi = gi * gi;
			// if(sparse && vargi > vargfilter){
			// 	int ii = 0;
			// 	while(vargi > vargfilter){
			// 		gi = norm_sample(rhs / v, sqrt(varei / v));
			// 		vargi = gi * gi;
			// 		ii++;
			// 		if(ii > 100){
			// 			gi = 0;
			// 			vargi = 0;
			// 		}
			// 	}
			// }
			sumvargi += vargi;
			sumvg += vargi * vx[i];

			// X'e = X'(y - Xb)
			if(gi != g[i]){
				start = ldm.begin_col(i);
    			end = ldm.end_col(i);
    			for (; start != end; ++start){
      				r_hat[start.row()] -= (gi - g[i]) * n * (*start);
   				}

				// update snp effect
				g[i] = gi;
			}
		}

        // if(sparse){
        // 	vargfilter = sum(vargx) / m + stddev(vargx) * 5;
        // }
        varg = (sumvargi + s2vara_ * dfvara_) / (dfvara_ + m);
		varg = invchisq_sample(m + dfvara_, varg);

		// genetic variance (betaX'Xbeta)
		vara_ = sum(g % (xy - r_hat)) / (n - 1);

		// sample residual variance from inv-chisq distribution
		vare_ = (yy - sum(g % (xy + r_hat)) + s2vare_ * dfvare_) / (n + dfvare_);
		vare_ = invchisq_sample(n + dfvare_, vare_);

		if(iter >= nburn){
			count++;
			sumvg_store[iter - nburn] = sumvg;
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            hsq_store[iter - nburn] = vara_ / (vara_ + vare_);
            g_store += g;

            // record pve for each window
            if(WPPA){
                for(R_xlen_t w = 0; w < nw; w++){
                	windxi = as<arma::uvec>(windx[w]); 
                    varw = n * spsubmatmtp(ldm, g, windxi);
                    varw /= (n - 1);
                    wgvei[w] += (varw / vara_);
                    if((varw / vara_) >= wppa){
                        wppai[w] += 1;
                    }
                }
            }

		}
		
		// print iteration details
		if(verbose){
			if((iter + 1) % outfreq == 0){
				timer_loop.step("end");
				NumericVector tdiff(timer_loop);
			    tt0 = (tdiff[1] - tdiff[0]) / (iter + 1) / 1e9;
			    tt = floor(tt0 * (niter - iter));
				hor = floor(tt / 3600);
				min = floor(tt % 3600 / 60);
				sec = floor(tt % 3600 % 60);
				Rcpp::Rcout << " " << iter + 1 << " ";
				Rcpp::Rcout << std::fixed << sumvg << " ";
				Rcpp::Rcout << std::fixed << vara_ << " ";
				Rcpp::Rcout << std::fixed << vare_ << " ";
				Rcpp::Rcout << std::fixed << (vara_ / (vara_ + vare_)) << " ";
				Rprintf("%02dh%02dm%02ds \n", hor, min, sec);
			}
		}
	}

	// get average of all parameters
    g = g_store / count;
    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    sumvg = mean(sumvg_store);
    double sumvgsd = sd(sumvg_store);
    vara_ = mean(vara_store);
    double varasd = sd(vara_store);
    vare_ = mean(vare_store);
    double varesd = sd(vare_store);
    hsq = mean(hsq_store);
    double hsqsd = sd(hsq_store);

    if(verbose){
        Rcpp::Rcout << "Posterior parameters:" << std::endl;
        Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << " ± " << std::fixed << sumvgsd << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << " ± " << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << " ± " << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated hsq " << std::fixed << hsq << " ± " << std::fixed << hsqsd << std::endl;
    }

    timer.step("end");
	NumericVector res(timer);
    tt = floor((res[1] - res[0]) / 1e9);
	hor = floor(tt / 3600);
	min = floor(tt % 3600 / 60);
	sec = floor(tt % 3600 % 60);
	if(verbose){Rprintf("Finished within total run time: %02dh%02dm%02ds \n", hor, min, sec);}
	  if(WPPA){
        return List::create(
            Named("vara") = vara_, 
            Named("vare") = vare_,
            Named("g") = g,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("vara") = vara_, 
            Named("vare") = vare_,
            Named("g") = g
        );
    }
}

// [[Rcpp::export]]
Rcpp::List SBayesA_spa(
	const NumericMatrix sumstat, 
	const arma::sp_mat & ldm,
	const int niter = 50000,
	const int nburn = 20000,
	const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vara = R_NilValue,
    const Nullable<double> dfvara = R_NilValue,
    const Nullable<double> s2vara = R_NilValue,
    const Nullable<double> vare = R_NilValue,
    const Nullable<double> dfvare = R_NilValue,
    const Nullable<double> s2vare = R_NilValue,
	const int outfreq = 100,
	const bool verbose = true
){

	int m;
	// arma::sp_mat ldm;
	// arma::mat tXX_ma;
	sp_mat::const_iterator start, end;
	if(sumstat.nrow() != ldm.n_rows){
		throw Rcpp::exception("Number of SNPs not equals.");
	}
	m = ldm.n_rows;
	int n = mean(na_omit(sumstat(_, 3)));
	bool WPPA = false;
	int count = 0;
	double xx, gi, rhs, lhs, uhat, v, hsq;
	double vara_, dfvara_, s2vara_, vare_, dfvare_, s2vare_, sumvg;
	arma::vec g(m); g.zeros();
	arma::vec g_store(m); g_store.zeros();
	arma::vec xpx(m); xpx.zeros();
	arma::vec vx(m); vx.zeros();
	for(int i = 0; i < m; i++){
		vx[i] = ldm(i, i);
		xpx[i] = vx[i] * n;
	}

	arma::vec xy(m);
	arma::vec yyi(m); yyi.zeros();
	NumericVector sumvg_store(niter - nburn + 1), vara_store(niter - nburn + 1), vare_store(niter - nburn + 1), hsq_store(niter - nburn + 1);
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
	if(dfvara.isNotNull()){
		dfvara_ = as<double>(dfvara);
	}else{
		dfvara_ = 4;
	}
	if(dfvara_ <= 2){
        throw Rcpp::exception("dfvara should not be less than 2.");
    }
	if(vara.isNotNull()){
		vara_ = as<double>(vara);
	}else{
		vara_ = ((dfvara_ - 2) / dfvara_) * yy / (n - 1);
	}
	if(vare.isNotNull()){
		vare_ = as<double>(vare);
	}else{
		vare_ = yy / (n - 1);
	}
	if(dfvare.isNotNull()){
		dfvare_ = as<double>(dfvare);
	}else{
		dfvare_ = -2;
	}
	if(s2vara.isNotNull()){
		s2vara_ = as<double>(s2vara);
	}else{
		s2vara_ = (vara_ / (sum(xpx) / (n - 1))) * (dfvara_ - 2) / dfvara_;
	}
	if(s2vare.isNotNull()){
		s2vare_ = as<double>(s2vare);
	}else{
		s2vare_ = 0;
	}
	arma::vec r_hat = xy;

	double varei;
	vec varediff(m);
	for(int i = 0; i < m; i++){
		start = ldm.begin_col(i);
		end = ldm.end_col(i);
		double j = 0;
		for (; start != end; ++start){
			j++;
		}
		varediff[i] = (m - j) / m;
	}

	double varg = (vara_ / (sum(xpx) / (n - 1)));

	if(niter < nburn){
		throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
	}

    // for gwas
    R_xlen_t nw;
    int ni;
    double varw;
    arma::vec windindx_(m);
    Rcpp::List windx(m);
    arma::uvec windxi;
    NumericVector wppai(m);
    NumericVector wgvei(m);
    if(windindx.isNotNull()){
        windindx_ = as<arma::vec>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = wppai[seq(0, (nw - 1))];
        wgvei = wgvei[seq(0, (nw - 1))];
        for(R_xlen_t w = 0; w < nw; w++){
            windx[w] = find(windindx_ == (w+1));
        }
    }

	if(verbose){
		Rcpp::Rcout.precision(4);
		Rcpp::Rcout << "Prior parameters:" << std::endl;
		Rcpp::Rcout << "    Model fitted at [SBayesA]" << std::endl;
		Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
		Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
		Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
		Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << s2vara_ << std::endl;
		Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
		Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << s2vare_ << std::endl;
		Rcpp::Rcout << "    Phenotypic var " << std::fixed << yy / (n - 1) << std::endl;
		if(WPPA){
            Rcpp::Rcout << "    Number of windows " << nw << std::endl;
            Rcpp::Rcout << "    GVE threshold for windows " << wppa << std::endl;
        }
		Rcpp::Rcout << "MCMC started: " << std::endl;
        Rcpp::Rcout << " Iter" << "\t";
        Rcpp::Rcout << "SigmaSq" << "\t";
		Rcpp::Rcout << "GenVar" << "\t";
		Rcpp::Rcout << "ResVar" << "\t";
		Rcpp::Rcout << "hsq" << "\t";
		Rcpp::Rcout << "Timeleft" << std::endl;
	}

	MyTimer timer;
    timer.step("start");
    MyTimer timer_loop;
    timer_loop.step("start");
    double tt0; 
    int tt, hor, min, sec;

	// MCMC procedure
	for(int iter = 0; iter < niter; iter++){

		sumvg = 0;

		// loop on snps
		for(int i = 0; i < m; i++){
			xx = xpx[i];
			gi = g[i];

			varei = varediff[i] * vara_ + vare_;

            // sample variance
            varg = (gi * gi + s2vara_ * dfvara_) / (dfvara_ + 1);
            varg = invchisq_sample(1 + dfvara_, varg);

			// right hand
			rhs = r_hat[i];
			if(gi){rhs += xx * gi;}

			// left hand
			lhs = xx / varei;
            v = xx + varei / varg;

            // sample snp effect
            gi = norm_sample(rhs / v, sqrt(varei / v));

            sumvg += (gi * gi) * vx[i];

			// X'e = X'(y - Xb)
			if(gi != g[i]){

				start = ldm.begin_col(i);
    			end = ldm.end_col(i);
    			for (; start != end; ++start){
      				r_hat[start.row()] -= (gi - g[i]) * n * (*start);
   				}

   				// update snp effect
				g[i] = gi;
			}
		}

		// genetic variance (betaX'Xbeta)
		vara_ = sum(g % (xy - r_hat)) / (n - 1);

		// sample residual variance from inv-chisq distribution
		vare_ = (yy - sum(g % (xy + r_hat)) + s2vare_ * dfvare_) / (n + dfvare_);
		vare_ = invchisq_sample(n + dfvare_, vare_);

		if(iter >= nburn){
			count++;
			sumvg_store[iter - nburn] = sumvg;
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            hsq_store[iter - nburn] = vara_ / (vara_ + vare_);
            g_store += g;

            // record pve for each window
            if(WPPA){
                for(R_xlen_t w = 0; w < nw; w++){
                	windxi = as<arma::uvec>(windx[w]); 
                    varw = n * spsubmatmtp(ldm, g, windxi);
                    varw /= (n - 1);
                    wgvei[w] += (varw / vara_);
                    if((varw / vara_) >= wppa){
                        wppai[w] += 1;
                    }
                }
            }
		}
		
		// print iteration details
		if(verbose){
			if((iter + 1) % outfreq == 0){
				timer_loop.step("end");
				NumericVector tdiff(timer_loop);
			    tt0 = (tdiff[1] - tdiff[0]) / (iter + 1) / 1e9;
			    tt = floor(tt0 * (niter - iter));
				hor = floor(tt / 3600);
				min = floor(tt % 3600 / 60);
				sec = floor(tt % 3600 % 60);
				Rcpp::Rcout << " " << iter + 1 << " ";
				Rcpp::Rcout << std::fixed << sumvg << " ";
				Rcpp::Rcout << std::fixed << vara_ << " ";
				Rcpp::Rcout << std::fixed << vare_ << " ";
				Rcpp::Rcout << std::fixed << (vara_ / (vara_ + vare_)) << " ";
				Rprintf("%02dh%02dm%02ds \n", hor, min, sec);
			}
		}
	}

	// get average of all parameters
    g = g_store / count;
    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    sumvg = mean(sumvg_store);
    double sumvgsd = sd(sumvg_store);
    vara_ = mean(vara_store);
    double varasd = sd(vara_store);
    vare_ = mean(vare_store);
    double varesd = sd(vare_store);
    hsq = mean(hsq_store);
    double hsqsd = sd(hsq_store);

    if(verbose){
        Rcpp::Rcout << "Posterior parameters:" << std::endl;
        Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << " ± " << std::fixed << sumvgsd << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << " ± " << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << " ± " << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated hsq " << std::fixed << hsq << " ± " << std::fixed << hsqsd << std::endl;
    }

    timer.step("end");
	NumericVector res(timer);
    tt = floor((res[1] - res[0]) / 1e9);
	hor = floor(tt / 3600);
	min = floor(tt % 3600 / 60);
	sec = floor(tt % 3600 % 60);
	if(verbose){Rprintf("Finished within total run time: %02dh%02dm%02ds \n", hor, min, sec);}
	  if(WPPA){
        return List::create(
            Named("vara") = vara_, 
            Named("vare") = vare_,
            Named("g") = g,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("vara") = vara_, 
            Named("vare") = vare_,
            Named("g") = g
        );
    }
}

// [[Rcpp::export]]
Rcpp::List SBayesBpi_spa(
	const NumericMatrix sumstat, 
	const arma::sp_mat & ldm,
	const double pi = 0.95,
	const int niter = 50000,
	const int nburn = 20000,
	const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vara = R_NilValue,
    const Nullable<double> dfvara = R_NilValue,
    const Nullable<double> s2vara = R_NilValue,
    const Nullable<double> vare = R_NilValue,
    const Nullable<double> dfvare = R_NilValue,
    const Nullable<double> s2vare = R_NilValue,
	const int outfreq = 100,
	const bool fixpi = false,
	const bool verbose = true
){

	int m;
	// arma::sp_mat ldm;
	// arma::mat tXX_ma;
	sp_mat::const_iterator start, end;
	if(sumstat.nrow() != ldm.n_rows){
		throw Rcpp::exception("Number of SNPs not equals.");
	}
	m = ldm.n_rows;
	int n = mean(na_omit(sumstat(_, 3)));
	bool WPPA = false;
	int count = 0, indistflag, NnzSnp;
	double xx, gi, rhs, lhs, logdetV, uhat, acceptProb, r, v, hsq;
	double vara_, dfvara_, s2vara_, vare_, dfvare_, s2vare_, sumvg;
	NumericVector pi_;
	NumericVector snptracker(m);
	NumericVector nzrate(m);
	arma::vec g(m); g.zeros();
	arma::vec g_store(m); g_store.zeros();
	arma::vec xpx(m); xpx.zeros();
	arma::vec vx(m); vx.zeros();
	for(int i = 0; i < m; i++){
		vx[i] = ldm(i, i);
		xpx[i] = vx[i] * n;
	}
	arma::vec xy(m);
	arma::vec yyi(m); yyi.zeros();
	NumericVector sumvg_store(niter - nburn + 1), vara_store(niter - nburn + 1), vare_store(niter - nburn + 1), hsq_store(niter - nburn + 1);
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
	double vary = yy / n;
    if(pi <= 0 || pi >= 1){
        throw Rcpp::exception("pi should be at 0 < pi < 1.");
    }else{
        pi_ = {pi, 1 - pi};
    }
	if(dfvara.isNotNull()){
		dfvara_ = as<double>(dfvara);
	}else{
		dfvara_ = 4;
	}
	if(dfvara_ <= 2){
        throw Rcpp::exception("dfvara should not be less than 2.");
    }
	if(vara.isNotNull()){
		vara_ = as<double>(vara);
	}else{
		vara_ = ((dfvara_ - 2) / dfvara_) * yy / (n - 1);
	}
	if(vare.isNotNull()){
		vare_ = as<double>(vare);
	}else{
		vare_ = yy / (n - 1);
	}
	if(dfvare.isNotNull()){
		dfvare_ = as<double>(dfvare);
	}else{
		dfvare_ = -2;
	}
	if(s2vara.isNotNull()){
		s2vara_ = as<double>(s2vara);
	}else{
		s2vara_ = (vara_ / ((1 - pi_[0]) * sum(xpx) / (n - 1))) * (dfvara_ - 2) / dfvara_;
	}
	if(s2vare.isNotNull()){
		s2vare_ = as<double>(s2vare);
	}else{
		s2vare_ = 0;
	}
	int n_fold = 2;
	NumericMatrix pi_store(niter - nburn + 1, n_fold);
	NumericVector fold_snp_num(n_fold);
	NumericVector logpi(n_fold), s(n_fold);
	arma::vec r_hat = xy;
	double varei;
	vec varediff(m);
	for(int i = 0; i < m; i++){
		start = ldm.begin_col(i);
		end = ldm.end_col(i);
		double j = 0;
		for (; start != end; ++start){
			j++;
		}
		varediff[i] = (m - j) / m;
	}

	double varg = (vara_ / ((1 - pi_[0]) * sum(xpx) / (n - 1)));

	if(niter < nburn){
		throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
	}

    // for gwas
    R_xlen_t nw;
    int ni;
    double varw;
    arma::vec windindx_(m);
    Rcpp::List windx(m);
    arma::uvec windxi;
    NumericVector wppai(m);
    NumericVector wgvei(m);
    if(windindx.isNotNull()){
        windindx_ = as<arma::vec>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = wppai[seq(0, (nw - 1))];
        wgvei = wgvei[seq(0, (nw - 1))];
        for(R_xlen_t w = 0; w < nw; w++){
            windx[w] = find(windindx_ == (w+1));
        }
    }

	if(verbose){
		Rcpp::Rcout.precision(4);
		Rcpp::Rcout << "Prior parameters:" << std::endl;
		if(fixpi){
			Rcpp::Rcout << "    Model fitted at [SBayesB]" << std::endl;
		}else{
			Rcpp::Rcout << "    Model fitted at [SBayesBπ]" << std::endl;
		}
		Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
		Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
		Rcpp::Rcout << "    Pi at " << pi_ << std::endl;
		Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
		Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << s2vara_ << std::endl;
		Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
		Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << s2vare_ << std::endl;
		Rcpp::Rcout << "    Phenotypic var " << std::fixed << yy / (n - 1) << std::endl;
		if(WPPA){
            Rcpp::Rcout << "    Number of windows " << nw << std::endl;
            Rcpp::Rcout << "    GVE threshold for windows " << wppa << std::endl;
        }
		Rcpp::Rcout << "MCMC started: " << std::endl;
		Rcpp::Rcout << " Iter" << "\t";
		Rcpp::Rcout << "NnzSnp" << "\t";
		for(int kk = 0; kk < n_fold; kk++){
			Rcpp::Rcout << "Pi" << kk + 1 << "\t";
		}
		Rcpp::Rcout << "SigmaSq" << "\t";
		Rcpp::Rcout << "GenVar" << "\t";
		Rcpp::Rcout << "ResVar" << "\t";
		Rcpp::Rcout << "hsq" << "\t";
		Rcpp::Rcout << "Timeleft" << std::endl;
	}

	MyTimer timer;
    timer.step("start");
    MyTimer timer_loop;
    timer_loop.step("start");
    double tt0; 
    int tt, hor, min, sec;

	// MCMC procedure
	for(int iter = 0; iter < niter; iter++){

		// log transform
		logpi = log(pi_);
		s[0] = logpi[0];
		sumvg = 0;

		// loop on snps
		for(int i = 0; i < m; i++){
			xx = xpx[i];
			gi = g[i];

			varei = varediff[i] * vara_ + vare_;

            // sample variance
            varg = (gi * gi + s2vara_ * dfvara_) / (dfvara_ + 1);
            varg = invchisq_sample(1 + dfvara_, varg);

			// right hand
			rhs = r_hat[i];
			if(gi){rhs += xx * gi;}

			// left hand
			lhs = xx / varei;

			logdetV = log(varg * lhs + 1);
			uhat = rhs / (xx + varei / varg);
			s[1] = -0.5 * (logdetV - (rhs * uhat / varei)) + logpi[1];
			acceptProb = 1 / sum(exp(s - s[0]));

			// group assignation
			r = uniform_sample();
			indistflag = r < acceptProb ? 0 : 1;
			snptracker[i] = indistflag;

			if(indistflag == 0){

				// zero effect
				gi = 0;
			}else{
				v = xx + varei / varg;

				// gibbs sample snp effect
				gi = norm_sample(rhs / v, sqrt(varei / v));

				sumvg += gi * gi * vx[i];
			}

			// X'e = X'(y - Xb)
			if(gi != g[i]){
				start = ldm.begin_col(i);
    			end = ldm.end_col(i);
    			for (; start != end; ++start){
      				r_hat[start.row()] -= (gi - g[i]) * n * (*start);
   				}

				// update snp effect
				g[i] = gi;
			}
		}

        fold_snp_num[1] = sum(snptracker);
        fold_snp_num[0] = m - fold_snp_num[1];
        NnzSnp = fold_snp_num[1];

		// genetic variance (betaX'Xbeta)
		vara_ = sum(g % (xy - r_hat)) / (n - 1);
	
		// sample residual variance from inv-chisq distribution
		vare_ = (yy - sum(g % (xy + r_hat)) + s2vare_ * dfvare_) / (n + dfvare_);
		vare_ = invchisq_sample(n + dfvare_, vare_);

		// update pi
		if(!fixpi){
			pi_ = rdirichlet_sample(n_fold, fold_snp_num + 1);
		}
		// Rcout << fold_snp_num << endl;

		if(iter >= nburn){
			count++;
			pi_store(iter - nburn, _) = pi_;
            sumvg_store[iter - nburn] = sumvg;
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            hsq_store[iter - nburn] = vara_ / (vara_ + vare_);
            g_store += g;
            nzrate += snptracker;

            // record pve for each window
            if(WPPA){
                for(R_xlen_t w = 0; w < nw; w++){
                	windxi = as<arma::uvec>(windx[w]); 
                    varw = n * spsubmatmtp(ldm, g, windxi);
                    varw /= (n - 1);
                    wgvei[w] += (varw / vara_);
                    if((varw / vara_) >= wppa){
                        wppai[w] += 1;
                    }
                }
            }

		}
		
		// print iteration details
		if(verbose){
			if((iter + 1) % outfreq == 0){
				timer_loop.step("end");
				NumericVector tdiff(timer_loop);
			    tt0 = (tdiff[1] - tdiff[0]) / (iter + 1) / 1e9;
			    tt = floor(tt0 * (niter - iter));
				hor = floor(tt / 3600);
				min = floor(tt % 3600 / 60);
				sec = floor(tt % 3600 % 60);
				Rcpp::Rcout << " " << iter + 1 << " ";
				Rcpp::Rcout << NnzSnp << " ";
				for(int kk = 0; kk < n_fold; kk++){
					Rcpp::Rcout << std::fixed << pi_[kk] << " ";
				}
				Rcpp::Rcout << std::fixed << sumvg << " ";
				Rcpp::Rcout << std::fixed << vara_ << " ";
				Rcpp::Rcout << std::fixed << vare_ << " ";
				Rcpp::Rcout << std::fixed << (vara_ / (vara_ + vare_)) << " ";
				Rprintf("%02dh%02dm%02ds \n", hor, min, sec);
			}
		}
	}

	// get average of all parameters
    g = g_store / count;
    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    nzrate = nzrate / count;
    for(int i = 0; i < n_fold; i++){
        pi_[i] = mean(pi_store(_, i));
    }
    sumvg = mean(sumvg_store);
    double sumvgsd = sd(sumvg_store);
    vara_ = mean(vara_store);
    double varasd = sd(vara_store);
    vare_ = mean(vare_store);
    double varesd = sd(vare_store);
    hsq = mean(hsq_store);
    double hsqsd = sd(hsq_store);

    if(verbose){
        Rcpp::Rcout << "Posterior parameters:" << std::endl;
        for(int i = 0; i < n_fold; i++){
            Rcpp::Rcout << "    Pi" << i + 1 << " " << std::fixed << pi_[i] << " ± " << std::fixed << sd(pi_store(_, i)) << std::endl;;
        }
        Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << " ± " << std::fixed << sumvgsd << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << " ± " << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << " ± " << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated hsq " << std::fixed << hsq << " ± " << std::fixed << hsqsd << std::endl;
    }

    timer.step("end");
	NumericVector res(timer);
    tt = floor((res[1] - res[0]) / 1e9);
	hor = floor(tt / 3600);
	min = floor(tt % 3600 / 60);
	sec = floor(tt % 3600 % 60);
	if(verbose){Rprintf("Finished within total run time: %02dh%02dm%02ds \n", hor, min, sec);}
	  if(WPPA){
        return List::create(
            Named("pi") = pi_, 
            Named("vara") = vara_, 
            Named("vare") = vare_,
            Named("g") = g,
            Named("nzrate") = nzrate,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("pi") = pi_, 
            Named("vara") = vara_, 
            Named("vare") = vare_,
            Named("g") = g,
            Named("nzrate") = nzrate
        );
    }
}

// [[Rcpp::export]]
Rcpp::List SBayesB_spa(
	const NumericMatrix sumstat, 
	const arma::sp_mat & ldm,
	const double pi = 0.95,
	const int niter = 50000,
	const int nburn = 20000,
	const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vara = R_NilValue,
    const Nullable<double> dfvara = R_NilValue,
    const Nullable<double> s2vara = R_NilValue,
    const Nullable<double> vare = R_NilValue,
    const Nullable<double> dfvare = R_NilValue,
    const Nullable<double> s2vare = R_NilValue,
	const int outfreq = 100,
	const bool verbose = true
){
	Rcpp::List res = SBayesBpi_spa(sumstat, ldm, pi, niter, nburn, windindx, wppa, vara, dfvara, s2vara, vare, dfvare, s2vare, outfreq, true, verbose);
	return res;
}

// [[Rcpp::export]]
Rcpp::List SBayesCpi_spa(
	const NumericMatrix sumstat, 
	const arma::sp_mat & ldm,
	const double pi = 0.95,
	const int niter = 50000,
	const int nburn = 20000,
	const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vara = R_NilValue,
    const Nullable<double> dfvara = R_NilValue,
    const Nullable<double> s2vara = R_NilValue,
    const Nullable<double> vare = R_NilValue,
    const Nullable<double> dfvare = R_NilValue,
    const Nullable<double> s2vare = R_NilValue,
	const int outfreq = 100,
	const bool fixpi = false,
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
	// arma::sp_mat ldm;
	// arma::mat tXX_ma;
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
	// 	tXX_ma = as<arma::mat>(ldm);
	// 	sparse = FALSE;
	// 	if(sumstat.nrow() != tXX_ma.n_rows){
	// 		throw Rcpp::exception("Number of SNPs not equals.");
	// 	}
	// 	m = tXX_ma.n_rows;
	// 	// for(int i = 0; i < m; i++){
	// 	// 	tXX_ma(i, i) = vx[i] * vx[i];
	// 	// 	for(int j = (i + 1); j < m; j++){
	// 	// 		tXX_ma(i, j) = tXX_ma(j, i) = tXX_ma(j, i) * vx[i] * vx[j];
	// 	// 	}
	// 	// }
	// }
	int n = mean(na_omit(sumstat(_, 3)));
	bool WPPA = false;
	int count = 0, indistflag, NnzSnp;
	double xx, gi, rhs, lhs, logdetV, uhat, acceptProb, r, v, hsq;
	double vara_, vargi, dfvara_, s2vara_, vare_, dfvare_, s2vare_, sumvargi, sumvg;
	NumericVector pi_;
	NumericVector snptracker(m);
	NumericVector nzrate(m);
	arma::vec g(m); g.zeros();
	arma::vec g_store(m); g_store.zeros();
	arma::vec xpx(m); xpx.zeros();
	arma::vec vx(m); vx.zeros();
	// if(sparse){
		for(int i = 0; i < m; i++){
			vx[i] = ldm(i, i);
			xpx[i] = vx[i] * n;
		}
	// }else{
	// 	vx = tXX_ma.diag();
	// 	xpx = tXX_ma.diag() * n;
	// }
	arma::vec xy(m);
	arma::vec yyi(m); yyi.zeros();
	NumericVector sumvg_store(niter - nburn + 1), vara_store(niter - nburn + 1), vare_store(niter - nburn + 1), hsq_store(niter - nburn + 1);
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
	double vary = yy / n;
    if(pi <= 0 || pi >= 1){
        throw Rcpp::exception("pi should be at 0 < pi < 1.");
    }else{
        pi_ = {pi, 1 - pi};
    }
	if(dfvara.isNotNull()){
		dfvara_ = as<double>(dfvara);
	}else{
		dfvara_ = 4;
	}
	if(dfvara_ <= 2){
        throw Rcpp::exception("dfvara should not be less than 2.");
    }
	if(vara.isNotNull()){
		vara_ = as<double>(vara);
	}else{
		vara_ = ((dfvara_ - 2) / dfvara_) * yy / (n - 1);
	}
	if(vare.isNotNull()){
		vare_ = as<double>(vare);
	}else{
		vare_ = yy / (n - 1);
	}
	if(dfvare.isNotNull()){
		dfvare_ = as<double>(dfvare);
	}else{
		dfvare_ = -2;
	}
	if(s2vara.isNotNull()){
		s2vara_ = as<double>(s2vara);
	}else{
		s2vara_ = (vara_ / ((1 - pi_[0]) * sum(xpx) / (n - 1))) * (dfvara_ - 2) / dfvara_;
	}
	if(s2vare.isNotNull()){
		s2vare_ = as<double>(s2vare);
	}else{
		s2vare_ = 0;
	}
	int n_fold = 2;
	NumericMatrix pi_store(niter - nburn + 1, n_fold);
	NumericVector fold_snp_num(n_fold);
	NumericVector logpi(n_fold), s(n_fold);
	arma::vec r_hat = xy;
	// arma::vec varei(m); varei.fill(vare_);
	// vec vargx(m);
    // double vargfilter = yy / (n - 1);
    double varei;
	vec varediff(m);
	// if(sparse){
		for(int i = 0; i < m; i++){
			start = ldm.begin_col(i);
			end = ldm.end_col(i);
			double j = 0;
			for (; start != end; ++start){
				j++;
			}
			varediff[i] = (m - j) / m;
		}
	// }
	double varg = (vara_ / ((1 - pi_[0]) * sum(xpx) / (n - 1)));

	if(niter < nburn){
		throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
	}

    // for gwas
    R_xlen_t nw;
    int ni;
    double varw;
    arma::vec windindx_(m);
    Rcpp::List windx(m);
    arma::uvec windxi;
    NumericVector wppai(m);
    NumericVector wgvei(m);
    if(windindx.isNotNull()){
        windindx_ = as<arma::vec>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = wppai[seq(0, (nw - 1))];
        wgvei = wgvei[seq(0, (nw - 1))];
        for(R_xlen_t w = 0; w < nw; w++){
            windx[w] = find(windindx_ == (w+1));
        }
    }

	if(verbose){
		Rcpp::Rcout.precision(4);
		Rcpp::Rcout << "Prior parameters:" << std::endl;
		if(fixpi){
			Rcpp::Rcout << "    Model fitted at [SBayesC]" << std::endl;
		}else{
			Rcpp::Rcout << "    Model fitted at [SBayesCπ]" << std::endl;
		}
		Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
		Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
		Rcpp::Rcout << "    Pi at " << pi_ << std::endl;
		Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
		Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << s2vara_ << std::endl;
		Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
		Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << s2vare_ << std::endl;
		Rcpp::Rcout << "    Phenotypic var " << std::fixed << yy / (n - 1) << std::endl;
		if(WPPA){
            Rcpp::Rcout << "    Number of windows " << nw << std::endl;
            Rcpp::Rcout << "    GVE threshold for windows " << wppa << std::endl;
        }
		Rcpp::Rcout << "MCMC started: " << std::endl;
		Rcpp::Rcout << " Iter" << "\t";
		Rcpp::Rcout << "NnzSnp" << "\t";
		for(int kk = 0; kk < n_fold; kk++){
			Rcpp::Rcout << "Pi" << kk + 1 << "\t";
		}
		Rcpp::Rcout << "SigmaSq" << "\t";
		Rcpp::Rcout << "GenVar" << "\t";
		Rcpp::Rcout << "ResVar" << "\t";
		Rcpp::Rcout << "hsq" << "\t";
		Rcpp::Rcout << "Timeleft" << std::endl;
	}

	MyTimer timer;
    timer.step("start");
    MyTimer timer_loop;
    timer_loop.step("start");
    double tt0; 
    int tt, hor, min, sec;

	// MCMC procedure
	for(int iter = 0; iter < niter; iter++){

		// log transform
		logpi = log(pi_);
		s[0] = logpi[0];
		sumvargi = 0;
		sumvg = 0;
		// vargx.fill(0);

		// loop on snps
		for(int i = 0; i < m; i++){
			xx = xpx[i];
			gi = g[i];
			// if(sparse){
			// 	vare_ = varei[i];
			// }
			// if(sparse){
				varei = varediff[i] * vara_ + vare_;
			// }else{
				// varei = vare_;
			// }

			// right hand
			rhs = r_hat[i];
			if(gi){rhs += xx * gi;}

			// left hand
			lhs = xx / varei;

			logdetV = log(varg * lhs + 1);
			uhat = rhs / (xx + varei / varg);
			s[1] = -0.5 * (logdetV - (rhs * uhat / varei)) + logpi[1];
			acceptProb = 1 / sum(exp(s - s[0]));	

			// group assignation
			r = uniform_sample();
			indistflag = r < acceptProb ? 0 : 1;
			snptracker[i] = indistflag;

			if(indistflag == 0){

				// zero effect
				gi = 0;
			}else{
				v = xx + varei / varg;

				// gibbs sample snp effect
				gi = norm_sample(rhs / v, sqrt(varei / v));
				vargi = gi * gi;
				if((vargi * vx[i]) > vary){
					int ii = 0;
					while((vargi * vx[i]) > vary){
						gi = norm_sample(rhs / v, sqrt(varei / v));
						vargi = gi * gi;
						ii++;
						if(ii > 100){
							gi = 0;
							vargi = 0;
						}
					}
				}
				sumvargi += vargi;
				sumvg += vargi * vx[i];
			}

			// X'e = X'(y - Xb)
			if(gi != g[i]){
				// if(sparse){
					start = ldm.begin_col(i);
	    			end = ldm.end_col(i);
	    			// double sumvare = 0;
	    			for (; start != end; ++start){
	    				// sumvare += (g[start.row()] * (xy[start.row()] + r_hat[start.row()]));
	    				r_hat[start.row()] -= (gi - g[i]) * (*start) * n;
	   				}
					// varei[i] = (yy - sumvare) / (n - 1);
   				// }else{
   					// r_hat -= (gi - g[i]) * tXX_ma.col(i) * n;
   				// }
   				// update snp effect
				g[i] = gi;
			}
		}

        fold_snp_num[1] = sum(snptracker);
        fold_snp_num[0] = m - fold_snp_num[1];
        NnzSnp = fold_snp_num[1];

        // if(sparse){
        // 	vargfilter = sum(vargx) / NnzSnp + stddev(vargx.elem(find(vargx))) * 5;
        // }
        varg = (sumvargi + s2vara_ * dfvara_) / (dfvara_ + NnzSnp);
		varg = invchisq_sample(NnzSnp + dfvara_, varg);

		// genetic variance (betaX'Xbeta)
		vara_ = sum(g % (xy - r_hat)) / (n - 1);

		// sample residual variance from inv-chisq distribution
		vare_ = yy - sum(g % (xy + r_hat));
		if(vare_ < 0)	vare_ = 0;
		vare_ = (vare_ + s2vare_ * dfvare_) / (n + dfvare_);
		vare_ = invchisq_sample(n + dfvare_, vare_);

		// update pi
		if(!fixpi){
			pi_ = rdirichlet_sample(n_fold, fold_snp_num + 1);
		}

		if(iter >= nburn){
			count++;
			pi_store(iter - nburn, _) = pi_;
            sumvg_store[iter - nburn] = sumvg;
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            hsq_store[iter - nburn] = vara_ / (vara_ + vare_);
            g_store += g;
            nzrate += snptracker;

            // record pve for each window
            if(WPPA){
                for(R_xlen_t w = 0; w < nw; w++){
                	windxi = as<arma::uvec>(windx[w]); 
                    // if(sparse){
                    	varw = n * spsubmatmtp(ldm, g, windxi);
                    // }else{
                    // 	varw = n * arma::as_scalar(g(windxi).t() * tXX_ma(windxi, windxi) * g(windxi));
                    // }
                    varw /= (n - 1);
                    wgvei[w] += (varw / vara_);
                    if((varw / vara_) >= wppa){
                        wppai[w] += 1;
                    }
                }
            }

		}
		
		// print iteration details
		if(verbose){
			if((iter + 1) % outfreq == 0){
				
				// Rcout << mean(vargx.elem(find(vargx))) << " ~ " << stddev(vargx.elem(find(vargx))) << " [" << vargx.elem(find(vargx)).min() << ", " << vargx.elem(find(vargx)).max() << "]" << endl;
		
				timer_loop.step("end");
				NumericVector tdiff(timer_loop);
			    tt0 = (tdiff[1] - tdiff[0]) / (iter + 1) / 1e9;
			    tt = floor(tt0 * (niter - iter));
				hor = floor(tt / 3600);
				min = floor(tt % 3600 / 60);
				sec = floor(tt % 3600 % 60);
				Rcpp::Rcout << " " << iter + 1 << " ";
				Rcpp::Rcout << NnzSnp << " ";
				for(int kk = 0; kk < n_fold; kk++){
					Rcpp::Rcout << std::fixed << pi_[kk] << " ";
				}
				Rcpp::Rcout << std::fixed << sumvg << " ";
				Rcpp::Rcout << std::fixed << vara_ << " ";
				Rcpp::Rcout << std::fixed << vare_ << " ";
				Rcpp::Rcout << std::fixed << (vara_ / (vara_ + vare_)) << " ";
				Rprintf("%02dh%02dm%02ds \n", hor, min, sec);
			}
		}
	}

	// get average of all parameters
    g = g_store / count;
    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    nzrate = nzrate / count;
    for(int i = 0; i < n_fold; i++){
        pi_[i] = mean(pi_store(_, i));
    }
    sumvg = mean(sumvg_store);
    double sumvgsd = sd(sumvg_store);
    vara_ = mean(vara_store);
    double varasd = sd(vara_store);
    vare_ = mean(vare_store);
    double varesd = sd(vare_store);
    hsq = mean(hsq_store);
    double hsqsd = sd(hsq_store);

    if(verbose){
        Rcpp::Rcout << "Posterior parameters:" << std::endl;
        for(int i = 0; i < n_fold; i++){
            Rcpp::Rcout << "    Pi" << i + 1 << " " << std::fixed << pi_[i] << " ± " << std::fixed << sd(pi_store(_, i)) << std::endl;;
        }
        Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << " ± " << std::fixed << sumvgsd << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << " ± " << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << " ± " << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated hsq " << std::fixed << hsq << " ± " << std::fixed << hsqsd << std::endl;
    }

    timer.step("end");
	NumericVector res(timer);
    tt = floor((res[1] - res[0]) / 1e9);
	hor = floor(tt / 3600);
	min = floor(tt % 3600 / 60);
	sec = floor(tt % 3600 % 60);
	if(verbose){Rprintf("Finished within total run time: %02dh%02dm%02ds \n", hor, min, sec);}
	  if(WPPA){
        return List::create(
            Named("pi") = pi_, 
            Named("vara") = vara_, 
            Named("vare") = vare_,
            Named("g") = g,
            Named("nzrate") = nzrate,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("pi") = pi_, 
            Named("vara") = vara_, 
            Named("vare") = vare_,
            Named("g") = g,
            Named("nzrate") = nzrate
        );
    }
}

// [[Rcpp::export]]
Rcpp::List SBayesC_spa(
	const NumericMatrix sumstat, 
	const arma::sp_mat & ldm,
	const double pi = 0.95,
	const int niter = 50000,
	const int nburn = 20000,
	const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vara = R_NilValue,
    const Nullable<double> dfvara = R_NilValue,
    const Nullable<double> s2vara = R_NilValue,
    const Nullable<double> vare = R_NilValue,
    const Nullable<double> dfvare = R_NilValue,
    const Nullable<double> s2vare = R_NilValue,
	const int outfreq = 100,
	const bool verbose = true
){
	Rcpp::List res = SBayesCpi_spa(sumstat, ldm, pi, niter, nburn, windindx, wppa, vara, dfvara, s2vara, vare, dfvare, s2vare, outfreq, true, verbose);
	return res;
}

// [[Rcpp::export]]
Rcpp::List SBayesLASSO_spa(
	const NumericMatrix sumstat, 
	const arma::sp_mat & ldm,
	const int niter = 50000,
	const int nburn = 20000,
	const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vara = R_NilValue,
    const Nullable<double> dfvara = R_NilValue,
    const Nullable<double> s2vara = R_NilValue,
    const Nullable<double> vare = R_NilValue,
    const Nullable<double> dfvare = R_NilValue,
    const Nullable<double> s2vare = R_NilValue,
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
	// arma::sp_mat ldm;
	// arma::mat tXX_ma;
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
	// 	tXX_ma = as<arma::mat>(ldm);
	// 	sparse = FALSE;
	// 	if(sumstat.nrow() != tXX_ma.n_rows){
	// 		throw Rcpp::exception("Number of SNPs not equals.");
	// 	}
	// 	m = tXX_ma.n_rows;
	// 	// for(int i = 0; i < m; i++){
	// 	// 	tXX_ma(i, i) = vx[i] * vx[i];
	// 	// 	for(int j = (i + 1); j < m; j++){
	// 	// 		tXX_ma(i, j) = tXX_ma(j, i) = tXX_ma(j, i) * vx[i] * vx[j];
	// 	// 	}
	// 	// }
	// }

	int n = mean(na_omit(sumstat(_, 3)));
	bool WPPA = false;
	int count = 0;
	double xx, gi, rhs, lhs, uhat, v, hsq;
	double vara_, vargi, dfvara_, s2vara_, vare_, dfvare_, s2vare_, sumvg;
	arma::vec g(m); g.zeros();
	arma::vec g_store(m); g_store.zeros();
	arma::vec xpx(m); xpx.zeros();
	arma::vec vx(m); vx.zeros();
	// if(sparse){
		for(int i = 0; i < m; i++){
			vx[i] = ldm(i, i);
			xpx[i] = vx[i] * n;
		}
	// }else{
	// 	vx = tXX_ma.diag();
	// 	xpx = tXX_ma.diag() * n;
	// }
	arma::vec xy(m);
	arma::vec yyi(m); yyi.zeros();
	NumericVector sumvg_store(niter - nburn + 1), vara_store(niter - nburn + 1), vare_store(niter - nburn + 1), hsq_store(niter - nburn + 1);
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
	double vary = yy / n;
	if(dfvara.isNotNull()){
		dfvara_ = as<double>(dfvara);
	}else{
		dfvara_ = 4;
	}
	if(dfvara_ <= 2){
        throw Rcpp::exception("dfvara should not be less than 2.");
    }
	if(vara.isNotNull()){
		vara_ = as<double>(vara);
	}else{
		vara_ = ((dfvara_ - 2) / dfvara_) * yy / (n - 1);
	}
	if(vare.isNotNull()){
		vare_ = as<double>(vare);
	}else{
		vare_ = yy / (n - 1);
	}
	if(dfvare.isNotNull()){
		dfvare_ = as<double>(dfvare);
	}else{
		dfvare_ = -2;
	}
	if(s2vara.isNotNull()){
		s2vara_ = as<double>(s2vara);
	}else{
		s2vara_ = (vara_ / (sum(xpx) / (n - 1))) * (dfvara_ - 2) / dfvara_;
	}
	if(s2vare.isNotNull()){
		s2vare_ = as<double>(s2vare);
	}else{
		s2vare_ = 0;
	}
	arma::vec r_hat = xy;
	// vec vargx(m);
    // double vargfilter = yy / (n - 1);
    double varei;
	vec varediff(m);
	// if(sparse){
		for(int i = 0; i < m; i++){
			start = ldm.begin_col(i);
			end = ldm.end_col(i);
			double j = 0;
			for (; start != end; ++start){
				j++;
			}
			varediff[i] = (m - j) / m;
		}
	// }
	NumericVector varg(m); varg.fill(vara_ / (sum(xpx) / (n - 1)));

	if(niter < nburn){
		throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
	}

    double R2 = (dfvara_ - 2) / dfvara_;
    double lambda2 = 2 * (1 - R2) / (R2) * (sum(xpx) / (n - 1));
    double lambda = sqrt(lambda2); 
    double shape, shape0 = 1.1;
    double rate, rate0 = (shape - 1) / lambda2;

    // for gwas
    R_xlen_t nw;
    int ni;
    double varw;
    arma::vec windindx_(m);
    Rcpp::List windx(m);
    arma::uvec windxi;
    NumericVector wppai(m);
    NumericVector wgvei(m);
    if(windindx.isNotNull()){
        windindx_ = as<arma::vec>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = wppai[seq(0, (nw - 1))];
        wgvei = wgvei[seq(0, (nw - 1))];
        for(R_xlen_t w = 0; w < nw; w++){
            windx[w] = find(windindx_ == (w+1));
        }
    }

	if(verbose){
		Rcpp::Rcout.precision(4);
		Rcpp::Rcout << "Prior parameters:" << std::endl;
		Rcpp::Rcout << "    Model fitted at [SBayesLASSO]" << std::endl;
		Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
		Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
		Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
		Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << s2vara_ << std::endl;
		Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
		Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << s2vare_ << std::endl;
		Rcpp::Rcout << "    Phenotypic var " << std::fixed << yy / (n - 1) << std::endl;
		if(WPPA){
            Rcpp::Rcout << "    Number of windows " << nw << std::endl;
            Rcpp::Rcout << "    GVE threshold for windows " << wppa << std::endl;
        }
		Rcpp::Rcout << "MCMC started: " << std::endl;
        Rcpp::Rcout << " Iter" << "\t";
        Rcpp::Rcout << "Lambda" << "\t";
		Rcpp::Rcout << "SigmaSq" << "\t";
		Rcpp::Rcout << "GenVar" << "\t";
		Rcpp::Rcout << "ResVar" << "\t";
		Rcpp::Rcout << "hsq" << "\t";
		Rcpp::Rcout << "Timeleft" << std::endl;
	}

	MyTimer timer;
    timer.step("start");
    MyTimer timer_loop;
    timer_loop.step("start");
    double tt0; 
    int tt, hor, min, sec;

	// MCMC procedure
	for(int iter = 0; iter < niter; iter++){

		sumvg = 0;

		// loop on snps
		for(int i = 0; i < m; i++){
			xx = xpx[i];
			gi = g[i];

			// if(sparse){
				varei = varediff[i] * vara_ + vare_;
			// }else{
				// varei = vare_;
			// }

			// right hand
			rhs = r_hat[i];
			if(gi){rhs += xx * gi;}

			// left hand
			lhs = xx / varei;
            v = xx + 1 / varg[i];

            // sample snp effect
            gi = norm_sample(rhs / v, sqrt(varei / v));
            if(abs(gi) < 1e-6){gi = 1e-6;}

            vargi = 1 / rinvgaussian_sample(sqrt(varei) * lambda / abs(gi), lambda2);
            if(vargi < 0){
                // nothing
            }else{
                varg[i] = vargi;
            }

            sumvg += (gi * gi) * vx[i];

			// X'e = X'(y - Xb)
			if(gi != g[i]){
				// if(sparse){
					start = ldm.begin_col(i);
	    			end = ldm.end_col(i);
	    			for (; start != end; ++start){
	      				r_hat[start.row()] -= (gi - g[i]) * n * (*start);
	   				}
   				// }else{
   					// r_hat -= (gi - g[i]) * n * tXX_ma.col(i);
   				// }

				// update snp effect
				g[i] = gi;
			}
		}

		// genetic variance (betaX'Xbeta)
		vara_ = sum(g % (xy - r_hat)) / (n - 1);

        shape = shape0 + m;
        rate = rate0 + sum(varg) / 2;
        lambda2 = gamma_sample(shape, 1 / rate);
        lambda = sqrt(lambda2);

		// sample residual variance from inv-chisq distribution
		vare_ = (yy - sum(g % (xy + r_hat)) + s2vare_ * dfvare_) / (n + dfvare_);
		vare_ = invchisq_sample(n + dfvare_, vare_);

		if(iter >= nburn){
			count++;
			sumvg_store[iter - nburn] = sumvg;
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            hsq_store[iter - nburn] = vara_ / (vara_ + vare_);
            g_store += g;

            // record pve for each window
            if(WPPA){
                for(R_xlen_t w = 0; w < nw; w++){
                	windxi = as<arma::uvec>(windx[w]); 
                    // if(sparse){
                    	varw = n * spsubmatmtp(ldm, g, windxi);
                    // }else{
                    	// varw = n * arma::as_scalar(g(windxi).t() * tXX_ma(windxi, windxi) * g(windxi));
                    // }
                    varw /= (n - 1);
                    wgvei[w] += (varw / vara_);
                    if((varw / vara_) >= wppa){
                        wppai[w] += 1;
                    }
                }
            }

		}
		
		// print iteration details
		if(verbose){
			if((iter + 1) % outfreq == 0){
				timer_loop.step("end");
				NumericVector tdiff(timer_loop);
			    tt0 = (tdiff[1] - tdiff[0]) / (iter + 1) / 1e9;
			    tt = floor(tt0 * (niter - iter));
				hor = floor(tt / 3600);
				min = floor(tt % 3600 / 60);
				sec = floor(tt % 3600 % 60);
				Rcpp::Rcout << " " << iter + 1 << " ";
				Rcpp::Rcout << std::fixed << lambda << " ";
				Rcpp::Rcout << std::fixed << sumvg << " ";
				Rcpp::Rcout << std::fixed << vara_ << " ";
				Rcpp::Rcout << std::fixed << vare_ << " ";
				Rcpp::Rcout << std::fixed << (vara_ / (vara_ + vare_)) << " ";
				Rprintf("%02dh%02dm%02ds \n", hor, min, sec);
			}
		}
	}

	// get average of all parameters
    g = g_store / count;
    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    sumvg = mean(sumvg_store);
    double sumvgsd = sd(sumvg_store);
    vara_ = mean(vara_store);
    double varasd = sd(vara_store);
    vare_ = mean(vare_store);
    double varesd = sd(vare_store);
    hsq = mean(hsq_store);
    double hsqsd = sd(hsq_store);

    if(verbose){
        Rcpp::Rcout << "Posterior parameters:" << std::endl;
        Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << " ± " << std::fixed << sumvgsd << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << " ± " << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << " ± " << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated hsq " << std::fixed << hsq << " ± " << std::fixed << hsqsd << std::endl;
    }

    timer.step("end");
	NumericVector res(timer);
    tt = floor((res[1] - res[0]) / 1e9);
	hor = floor(tt / 3600);
	min = floor(tt % 3600 / 60);
	sec = floor(tt % 3600 % 60);
	if(verbose){Rprintf("Finished within total run time: %02dh%02dm%02ds \n", hor, min, sec);}
	  if(WPPA){
        return List::create(
            Named("vara") = vara_, 
            Named("vare") = vare_,
            Named("g") = g,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("vara") = vara_, 
            Named("vare") = vare_,
            Named("g") = g
        );
    }
}

// [[Rcpp::export]]
Rcpp::List SBayesR_spa(
	const NumericMatrix sumstat, 
	const arma::sp_mat & ldm,
	const Nullable<NumericVector> pi = R_NilValue,
	const Nullable<NumericVector> fold = R_NilValue,
	const int niter = 50000,
	const int nburn = 20000,
	const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vara = R_NilValue,
    const Nullable<double> dfvara = R_NilValue,
    const Nullable<double> s2vara = R_NilValue,
    const Nullable<double> vare = R_NilValue,
    const Nullable<double> dfvare = R_NilValue,
    const Nullable<double> s2vare = R_NilValue,
	const int outfreq = 100,
	const bool fixpi = false,
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
	// arma::sp_mat ldm;
	// // MSpMat ldm = as<MSpMat>(ldm);
	// arma::mat tXX_ma;
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
	// 	tXX_ma = as<arma::mat>(ldm);
	// 	sparse = FALSE;
	// 	if(sumstat.nrow() != tXX_ma.n_rows){
	// 		throw Rcpp::exception("Number of SNPs not equals.");
	// 	}
	// 	m = tXX_ma.n_rows;
	// 	// for(int i = 0; i < m; i++){
	// 	// 	tXX_ma(i, i) = vx[i] * vx[i];
	// 	// 	for(int j = (i + 1); j < m; j++){
	// 	// 		tXX_ma(i, j) = tXX_ma(j, i) = tXX_ma(j, i) * vx[i] * vx[j];
	// 	// 	}
	// 	// }
	// }

	int n = mean(na_omit(sumstat(_, 3)));
	bool WPPA = false;
	int count = 0, indistflag, NnzSnp;
	double xx, gi, rhs, lhs, logdetV, uhat, acceptProb, r, v, varg, hsq;
	double vara_, vargi, dfvara_, s2vara_, vare_, dfvare_, s2vare_, sumvg, sumvargi;
	NumericVector pi_, fold_;
	IntegerVector snptracker(m);
	NumericVector nzrate(m);
	arma::vec g(m); g.zeros();
	arma::vec g_store(m); g_store.zeros();
	arma::vec xpx(m); xpx.zeros();
	arma::vec vx(m); vx.zeros();
	// if(sparse){
		for(int i = 0; i < m; i++){
			vx[i] = ldm(i, i);
			xpx[i] = vx[i] * n;
		}
	// }else{
	// 	vx = tXX_ma.diag();
	// 	xpx = tXX_ma.diag() * n;
	// }
	arma::vec xy(m);
	arma::vec yyi(m); yyi.zeros();
	NumericVector sumvg_store(niter - nburn + 1), vara_store(niter - nburn + 1), vare_store(niter - nburn + 1), hsq_store(niter - nburn + 1);
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
	double vary = yy / n;
	if(pi.isNotNull()){
		pi_ = as<NumericVector>(pi);
	}else{
		pi_ = {0.95, 0.02, 0.02, 0.01};
	}
    if(pi_.length() < 2){
        throw Rcpp::exception("pi should be a vector.");
    }
	if(sum(pi_) != 1){
		throw Rcpp::exception("sum of pi should be 1.");
	}
	for(int i = 0; i < pi_.length(); i++){
		if(pi_[i] <= 0 || pi_[i] >= 1){
			throw Rcpp::exception("elements of pi should be at (0, 1)");
		}
	}
	if(fold.isNotNull()){
		fold_ = as<NumericVector>(fold);
	}else{
		fold_ = {0, 0.0001, 0.001, 0.01};
	}
	if(fold_.length() != pi_.length()){
		throw Rcpp::exception("length of pi and fold not equals.");
	}
	if(dfvara.isNotNull()){
		dfvara_ = as<double>(dfvara);
	}else{
		dfvara_ = 4;
	}
	if(dfvara_ <= 2){
        throw Rcpp::exception("dfvara should not be less than 2.");
    }
	if(vara.isNotNull()){
		vara_ = as<double>(vara);
	}else{
		vara_ = ((dfvara_ - 2) / dfvara_) * yy / (n - 1);
	}
	if(vare.isNotNull()){
		vare_ = as<double>(vare);
	}else{
		vare_ = yy / (n - 1);
	}
	if(dfvare.isNotNull()){
		dfvare_ = as<double>(dfvare);
	}else{
		dfvare_ = -2;
	}
	if(s2vara.isNotNull()){
		s2vara_ = as<double>(s2vara);
	}else{
		s2vara_ = (vara_ / ((1 - pi_[0]) * sum(xpx) / (n - 1))) * (dfvara_ - 2) / dfvara_;
	}
	if(s2vare.isNotNull()){
		s2vare_ = as<double>(s2vare);
	}else{
		s2vare_ = 0;
	}
	int n_fold = fold_.length();
	NumericVector stemp(n_fold);
	NumericMatrix pi_store(niter - nburn + 1, n_fold);
	NumericVector fold_snp_num(n_fold);
	NumericVector logpi(n_fold), s(n_fold);
	NumericVector vara_fold = (vara_ / ((1 - pi_[0]) * sum(xpx) / (n - 1))) * fold_;
	arma::vec r_hat = xy;
	// arma::vec varei(m); varei.fill(vare_);
    // vec vargx(m);
    // double vargfilter = yy / (n - 1);

	double varei;
	vec varediff(m);
	// if(sparse){
		for(int i = 0; i < m; i++){
			start = ldm.begin_col(i);
			end = ldm.end_col(i);
			double j = 0;
			for (; start != end; ++start){
				j++;
			}
			varediff[i] = (m - j) / m;
		}
	// }

	if(niter < nburn){
		throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
	}

    // for gwas
    R_xlen_t nw;
    int ni;
    double varw;
    arma::vec windindx_(m);
    Rcpp::List windx(m);
    arma::uvec windxi;
    NumericVector wppai(m);
    NumericVector wgvei(m);
    if(windindx.isNotNull()){
        windindx_ = as<arma::vec>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = wppai[seq(0, (nw - 1))];
        wgvei = wgvei[seq(0, (nw - 1))];
        for(R_xlen_t w = 0; w < nw; w++){
            windx[w] = find(windindx_ == (w+1));
        }
    }

	if(verbose){
		Rcpp::Rcout.precision(4);
		Rcpp::Rcout << "Prior parameters:" << std::endl;
		Rcpp::Rcout << "    Model fitted at [SBayesR]" << std::endl;
		Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
		Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
		Rcpp::Rcout << "    Pi at " << pi_ << std::endl;
		Rcpp::Rcout << "    Group fold " << fold_ << std::endl;
		Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
		Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << s2vara_ << std::endl;
		Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
		Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << s2vare_ << std::endl;
		Rcpp::Rcout << "    Phenotypic var " << std::fixed << yy / (n - 1) << std::endl;
		if(WPPA){
            Rcpp::Rcout << "    Number of windows " << nw << std::endl;
            Rcpp::Rcout << "    GVE threshold for windows " << wppa << std::endl;
        }
		Rcpp::Rcout << "MCMC started: " << std::endl;
		Rcpp::Rcout << " Iter" << "\t";
		Rcpp::Rcout << "NnzSnp" << "\t";
		for(int kk = 0; kk < n_fold; kk++){
			Rcpp::Rcout << "Pi" << kk + 1 << "\t";
		}
		Rcpp::Rcout << "SigmaSq" << "\t";
		Rcpp::Rcout << "GenVar" << "\t";
		Rcpp::Rcout << "ResVar" << "\t";
		Rcpp::Rcout << "hsq" << "\t";
		Rcpp::Rcout << "Timeleft" << std::endl;
	}

	MyTimer timer;
    timer.step("start");
    MyTimer timer_loop;
    timer_loop.step("start");
    double tt0; 
    int tt, hor, min, sec;

	// MCMC procedure
	for(int iter = 0; iter < niter; iter++){

		// log transform
		logpi = log(pi_);
		s[0] = logpi[0];
		sumvg = 0;
		sumvargi = 0;
		// vargx.fill(0);

		// loop on snps
		for(int i = 0; i < m; i++){
			xx = xpx[i];
			gi = g[i];
			// if(sparse){
			// 	vare_ = varei[i];
			// }

			// if(sparse){
				varei = varediff[i] * vara_ + vare_;
			// }else{
			// 	varei = vare_;
			// }

			// right hand
			rhs = r_hat[i];
			if(gi){rhs += xx * gi;}

			// left hand
			lhs = xx / varei;

			for(int j = 1; j < n_fold; j++){
				logdetV = log(vara_fold[j] * lhs + 1);
				uhat = rhs / (xx + varei / vara_fold[j]);
				s[j] = -0.5 * (logdetV - (rhs * uhat / varei)) + logpi[j];
			}

			for(int j = 0; j < n_fold; j++){
				double temp = 0.0;
				for(int k = 0; k < n_fold; k++){
					temp += exp(s[k] - s[j]);
				}
				stemp[j] = 1 / temp;
			}

			// group assignation
			acceptProb = 0;
			indistflag = 0;
			r = uniform_sample();

			for(int j = 0; j < n_fold; j++){
				acceptProb += stemp[j];
				if(r < acceptProb){
					indistflag = j;
					break;
				}
			}

			snptracker[i] = indistflag;
			if(indistflag == 0){

				// zero effect
				gi = 0;
			}else{
				v = xx + varei / vara_fold[indistflag];

				// gibbs sample snp effect
				gi = norm_sample(rhs / v, sqrt(varei / v));

				vargi = gi * gi;
				if((vargi * vx[i]) > vary){
					int ii = 0;
					while(((vargi * vx[i]) > vary)){
						gi = norm_sample(rhs / v, sqrt(varei / v));
						vargi = gi * gi;
						ii++;
						if(ii > 100){
							gi = 0;
							vargi = 0;
						}
					}
				}
				sumvargi += vargi / fold_[indistflag];
				sumvg += vargi * vx[i];
			}

			// X'e = X'(y - Xb)
			if(gi != g[i]){
				// if(sparse){
					start = ldm.begin_col(i);
	    			end = ldm.end_col(i);
	    			// double sumvare = 0;
	    			// for (InIterMat it(ldm, i); it; ++it){
	    			// 	sumvare += (g[it.index()] * (xy[it.index()] + r_hat[it.index()]));
	      	// 			r_hat[it.index()] -= (gi - g[i]) * (it.value());
        // 			}
	    			for (; start != end; ++start){
	    				// sumvare += (g[start.row()] * (xy[start.row()] + r_hat[start.row()]));
	      				r_hat[start.row()] -= (gi - g[i]) * n * (*start);
	   				}
	   				// varei[i] = (yy - sumvare) / (n - 1);
   				// }else{
   				// 	r_hat = r_hat - (gi - g[i]) * n * tXX_ma.col(i);
   				// }
   				// update snp effect
				g[i] = gi;
			}
		}

		for(int j = 0; j < n_fold; j++){
			fold_snp_num[j] = sum(snptracker == j);
		}

		NnzSnp = m - fold_snp_num[0];

        // if(sparse){
        // 	vargfilter = sum(vargx) / NnzSnp + stddev(vargx.elem(find(vargx))) * 5;
        // }
        varg = (sumvargi + s2vara_ * dfvara_) / (dfvara_ + NnzSnp);
		varg = invchisq_sample(NnzSnp + dfvara_, varg);

		// genetic variance (betaX'Xbeta)
		vara_ = sum(g % (xy - r_hat)) / (n - 1);
	
		// sample residual variance from inv-chisq distribution
		vare_ = (yy - sum(g % (xy + r_hat)) + s2vare_ * dfvare_) / (n + dfvare_);
		if(vare_ < 0)	vare_ = 0;
		vare_ = invchisq_sample(n + dfvare_, vare_);

		for(int j = 0; j < n_fold; j++){
			vara_fold[j] = varg * fold_[j];	
			// vara_fold[j] = vara_ * fold_[j];	
		}
		
		// update pi
		if(!fixpi){
			pi_ = rdirichlet_sample(n_fold, fold_snp_num + 1);
		}
		// Rcout << fold_snp_num << endl;

		if(iter >= nburn){
			count++;
			pi_store(iter - nburn, _) = pi_;
            sumvg_store[iter - nburn] = sumvg;
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            hsq_store[iter - nburn] = vara_ / (vara_ + vare_);
            g_store += g;
            for(int i = 0; i < m; i++){
                if(snptracker[i]){
                    nzrate[i] += 1;
                }
            }

            // record pve for each window
            if(WPPA){
                for(R_xlen_t w = 0; w < nw; w++){
                	windxi = as<arma::uvec>(windx[w]); 
                    // if(sparse){
                    	varw = n * spsubmatmtp(ldm, g, windxi);
                    // }else{
                    // 	varw = n * arma::as_scalar(g(windxi).t() * tXX_ma(windxi, windxi) * g(windxi));
                    // }
                    varw /= (n - 1);
                    wgvei[w] += (varw / vara_);
                    if((varw / vara_) >= wppa){
                        wppai[w] += 1;
                    }
                }
            }

		}
		
		// print iteration details
		if(verbose){
			if((iter + 1) % outfreq == 0){

				// Rcout << mean(vargx.elem(find(vargx))) << " ~ " << stddev(vargx.elem(find(vargx))) << " [" << vargx.elem(find(vargx)).min() << ", " << vargx.elem(find(vargx)).max() << "]" << endl;
		
				timer_loop.step("end");
				NumericVector tdiff(timer_loop);
			    tt0 = (tdiff[1] - tdiff[0]) / (iter + 1) / 1e9;
			    tt = floor(tt0 * (niter - iter));
				hor = floor(tt / 3600);
				min = floor(tt % 3600 / 60);
				sec = floor(tt % 3600 % 60);
				Rcpp::Rcout << " " << iter + 1 << " ";
				Rcpp::Rcout << NnzSnp << " ";
				for(int kk = 0; kk < n_fold; kk++){
					Rcpp::Rcout << std::fixed << pi_[kk] << " ";
				}
				Rcpp::Rcout << std::fixed << sumvg << " ";
				Rcpp::Rcout << std::fixed << vara_ << " ";
				Rcpp::Rcout << std::fixed << vare_ << " ";
				Rcpp::Rcout << std::fixed << (vara_ / (vara_ + vare_)) << " ";
				Rprintf("%02dh%02dm%02ds \n", hor, min, sec);
			}
		}
	}

	// get average of all parameters
    g = g_store / count;
    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    nzrate = nzrate / count;
    for(int i = 0; i < n_fold; i++){
        pi_[i] = mean(pi_store(_, i));
    }
    sumvg = mean(sumvg_store);
    double sumvgsd = sd(sumvg_store);
    vara_ = mean(vara_store);
    double varasd = sd(vara_store);
    vare_ = mean(vare_store);
    double varesd = sd(vare_store);
    hsq = mean(hsq_store);
    double hsqsd = sd(hsq_store);

    if(verbose){
        Rcpp::Rcout << "Posterior parameters:" << std::endl;
        for(int i = 0; i < n_fold; i++){
            Rcpp::Rcout << "    Pi" << i + 1 << " " << std::fixed << pi_[i] << " ± " << std::fixed << sd(pi_store(_, i)) << std::endl;;
        }
        Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << " ± " << std::fixed << sumvgsd << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << " ± " << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << " ± " << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated hsq " << std::fixed << hsq << " ± " << std::fixed << hsqsd << std::endl;
    }

    timer.step("end");
	NumericVector res(timer);
    tt = floor((res[1] - res[0]) / 1e9);
	hor = floor(tt / 3600);
	min = floor(tt % 3600 / 60);
	sec = floor(tt % 3600 % 60);
	if(verbose){Rprintf("Finished within total run time: %02dh%02dm%02ds \n", hor, min, sec);}
	  if(WPPA){
        return List::create(
            Named("pi") = pi_, 
            Named("vara") = vara_, 
            Named("vare") = vare_,
            Named("g") = g,
            Named("nzrate") = nzrate,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("pi") = pi_, 
            Named("vara") = vara_, 
            Named("vare") = vare_,
            Named("g") = g,
            Named("nzrate") = nzrate
        );
    }
}
