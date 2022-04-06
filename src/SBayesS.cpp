#include "hibayes.h"
#include "solver.h"

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
Rcpp::List SBayesS(
	arma::mat sumstat, 
	arma::sp_mat ldm,
    std::string model,
	arma::vec Pi,
	const int niter = 50000,
	const int nburn = 20000,
    const Nullable<arma::vec> fold = R_NilValue,
	const Nullable<arma::uvec> windindx = R_NilValue,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
	const int outfreq = 100,
	const int threads = 0,
	const bool verbose = true
){

    omp_setup(threads);

    int model_index = (model == "BayesRR" ? 1 : (model == "BayesA" ? 2 : (model == "BayesB" || model == "BayesBpi" ? 3 : (model == "BayesC" || model == "BayesCpi" ? 4 : (model == "BayesL" ? 5 : 6)))));
    if(sumstat.n_rows != ldm.n_rows){
		throw Rcpp::exception("Number of SNPs not equals.");
	}
    int m = ldm.n_rows;
    vec N_col = sumstat.col(3);
    int n = mean(N_col.elem(find_finite(N_col)));
    bool fixpi = false;
    if(model == "BayesB" || model == "BayesC")    fixpi = true;
    if(Pi.n_elem < 2)  throw Rcpp::exception("Pi should be a vector.");
    if(sum(Pi) != 1)   throw Rcpp::exception("sum of Pi should be 1.");
    if(Pi[0] == 1) throw Rcpp::exception("all markers have no effect size.");
    for(int i = 0; i < Pi.n_elem; i++){
        if(Pi[i] < 0 || Pi[i] > 1){
            throw Rcpp::exception("elements of Pi should be at the range of [0, 1]");
        }
    }
    vec fold_;
    if(fold.isNotNull()){
        fold_ = as<arma::vec>(fold);
    }else{
        if(model == "BayesR")    throw Rcpp::exception("'fold' should be provided for BayesR model.");
        fold_.resize(2);
    }
    if(fold_.n_elem != Pi.n_elem){
        throw Rcpp::exception("length of Pi and fold not equals.");
    }

    int n_records = (niter - nburn) / outfreq;

    int count = 0;
    int nzct = 0;
    int inc = 1;
    int n_fold = fold_.n_elem;
    int NnzSnp, indistflag;
    // double doc = 1.0;
    double xx, gi, gi_, rhs, lhs, logdetV, acceptProb, uhat, v;
    double vara_, dfvara_, s2vara_, vare_, dfvare_, s2vare_, vargi, hsq, s2varg_;
    vec snptracker;
    vec nzrate;
    if(model == "BayesRR" || model == "BayesA" || model == "BayesL"){
        NnzSnp = m;
        Pi[0] = 0; Pi[1] = 1;
        fixpi = true;
    }else{
        if(model != "BayesR" && Pi.n_elem != 2) throw Rcpp::exception("length of Pi should be 2, the first value is the proportion of non-effect markers.");
        nzrate.zeros(m); 
        snptracker.zeros(m); 
    }
    vec xy = zeros(m);
    vec r_hat = zeros(m);
	vec tmp = zeros(m);
	vec yyi = zeros(m);
    vec g = zeros(m);
    mat g_store = zeros(m, n_records);
    vec u = zeros(n);
    vec xpx = zeros(m);
    vec vx = zeros(m);
    // vec sumvg_store = zeros(niter - nburn + 1);
    vec vara_store = zeros(n_records);
    vec vare_store = zeros(n_records);
    vec hsq_store = zeros(n_records);
    mat pi_store;
    pi_store.resize(n_fold, n_records);

    #pragma omp parallel for
    for(int i = 0; i < m; i++){
		vx[i] = ldm(i, i);
		xpx[i] = vx[i] * n;
	}
    int count_y = 0;
    std::vector<bool> ifest(m); std::fill(ifest.begin(), ifest.end(), true);
    for(int k = 0; k < m; k++){
		if(Rcpp::NumericVector::is_na(sumstat(k, 1)) || Rcpp::NumericVector::is_na(sumstat(k, 2)) || Rcpp::NumericVector::is_na(sumstat(k, 3))){
			ifest[k] = false;
		}else{
            xy[k] = xpx[k] * sumstat(k, 1);
            r_hat[k] = xy[k];
			yyi[k] = xpx[k] * (sumstat(k, 1) * sumstat(k, 1) + (sumstat(k, 3) - 2) * sumstat(k, 2) * sumstat(k, 2));
			count_y++;
		}
	}
	if(count_y == 0){
		throw Rcpp::exception("Lack of SE.");
	}
    vec varediff = zeros(m);
    sp_mat::const_iterator start, end;
	for(int i = 0; i < m; i++){
		start = ldm.begin_col(i);
		end = ldm.end_col(i);
		double j = 0;
		for (; start != end; ++start){
			j++;
		}
		varediff[i] = (m - j) / m;
	}
	double yy = sum(yyi) / count_y;
    double vary = yy / (n - 1);
    double h2 = 0.5;
    double varei;
    
    if(dfvg.isNotNull()){
        dfvara_ = as<double>(dfvg);
    }else{
        dfvara_ = 4;
    }
    if(dfvara_ <= 2){
        throw Rcpp::exception("dfvg should not be less than 2.");
    }
    if(vg.isNotNull()){
        vara_ = as<double>(vg);
    }else{
        vara_ = ((dfvara_ - 2) / dfvara_) * vary * h2;
    }
    if(ve.isNotNull()){
        vare_ = as<double>(ve);
    }else{
        vare_ = vary * (1 - h2);
    }
    if(dfve.isNotNull()){
        dfvare_ = as<double>(dfve);
    }else{
        dfvare_ = -2;
    }
    if(s2vg.isNotNull()){
        s2vara_ = as<double>(s2vg);
    }else{
        s2vara_ = vara_ * (dfvara_ - 2) / dfvara_;
    }
    double varg = vara_ / ((1 - Pi[0]) * sum(vx));
    s2varg_ = s2vara_ / ((1 - Pi[0]) * sum(vx));
    if(s2ve.isNotNull()){
        s2vare_ = as<double>(s2ve);
    }else{
        s2vare_ = 0;
    }
    if(niter < nburn){
        throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
    }
    double R2 = (dfvara_ - 2) / dfvara_;
    double lambda2 = 2 * (1 - R2) / (R2) * sum(vx);
    double lambda = sqrt(lambda2); 
    double shape, shape0 = 1.1;
    double rate, rate0 = (shape0 - 1) / lambda2;
    vec vargL;
    if(model == "BayesL"){
        vargL.resize(m);
        vargL.fill(varg);
    }
    vec stemp = zeros(n_fold);
    vec fold_snp_num = zeros(n_fold);
    vec logpi = zeros(n_fold);
    vec s = zeros(n_fold);
    vec vara_fold = (vara_ / ((1 - Pi[0]) * sum(vx))) * fold_;
    
    // for gwas
    int nw;
    bool WPPA = false;
    uvec windindx_;
    vector<arma::uvec> windx;
    uvec windxi;
    vec wppai;
    if(windindx.isNotNull()){
        windindx_ = as<arma::uvec>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai.zeros(nw);
        for(int w = 0; w < nw; w++){
            windx.push_back(find(windindx_ == (w + 1)));
        }
    }

    if(verbose){
        Rcpp::Rcout.precision(4);
        Rcpp::Rcout << "Prior parameters:" << std::endl;
        Rcpp::Rcout << "    Model fitted at [" << (model == "BayesRR" ? "Bayes Ridge Regression" : model) << "]" << std::endl;
        Rcpp::Rcout << "    Number of observations " << n << std::endl;
        Rcpp::Rcout << "    Number of markers " << m << std::endl;
        Rcpp::Rcout << "    Number of markers used for analysis" << count_y << std::endl;
        for(int i = 0; i < Pi.n_elem; i++){
            if(i == 0){
                Rcpp::Rcout << "    π for markers in zero effect size " << Pi[i] << endl;
            }else{
                if(i == 1){
                    Rcpp::Rcout << "    π for markers in non-zero effect size " << Pi[i] << " ";
                }else{
                    Rcpp::Rcout << Pi[i] << " ";
                }
            }
        }
        Rcpp::Rcout << std::endl;
        if(model == "BayesR"){
            Rcpp::Rcout << "    Group fold ";
            for(int i = 0; i < fold_.n_elem; i++){
                Rcpp::Rcout << fold_[i] << " ";
            }
            Rcpp::Rcout << std::endl;
        }
        Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
        Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
        Rcpp::Rcout << "    Phenotypic var " << std::fixed << vary << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << std::fixed << s2vara_ << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << std::fixed << s2vare_ << std::endl;
        Rcpp::Rcout << "    Marker var " << std::fixed << varg << std::endl;
        Rcpp::Rcout << "    Inv-Chisq alpar " << std::fixed << dfvara_ << " " << std::fixed << s2varg_ << std::endl;
        if(WPPA){
            Rcpp::Rcout << "    Number of windows for GWAS analysis " << nw << std::endl;
        }
        Rcpp::Rcout << "MCMC started: " << std::endl;
        Rcpp::Rcout << " Iter" << "  ";
        Rcpp::Rcout << "NumNZSnp" << "  ";
        for(int i = 0; i < n_fold; i++){
            Rcpp::Rcout << "π" << i + 1 << "  ";
        }
        if(model == "BayesL")    Rcpp::Rcout << "Lambda" << "  ";
        Rcpp::Rcout << "Vg" << "  ";
        Rcpp::Rcout << "Ve" << "  ";
        Rcpp::Rcout << "h2" << "  ";
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

        switch(model_index){
            case 1:
                for(int i = 0; i < m; i++){
                    if(!ifest[i])   continue;
                    xx = xpx[i];
                    gi = g[i];
                    varei = varediff[i] * vara_ + vare_;
                    rhs = r_hat[i];
                    if(gi){rhs += xx * gi;}
                    lhs = xx / varei;
                    v = xx + varei / varg;
                    gi = norm_sample(rhs / v, sqrt(varei / v));
                    gi_ = (g[i] - gi) * n;
                    start = ldm.begin_col(i);
	    			end = ldm.end_col(i);
	    			for (; start != end; ++start){
	    				r_hat[start.row()] += gi_ * (*start);
	   				}
                    g[i] = gi;
                }
                varg = (ddot_(&m, g.memptr(), &inc, g.memptr(), &inc) + s2varg_ * dfvara_) / (dfvara_ + count_y);
                varg = invchisq_sample(count_y + dfvara_, varg);
                break;
            case 2:
                for(int i = 0; i < m; i++){
                    if(!ifest[i])   continue;
                    xx = xpx[i];
                    gi = g[i];
                    varei = varediff[i] * vara_ + vare_;
                    varg = (gi * gi + s2varg_ * dfvara_) / (dfvara_ + 1);
                    varg = invchisq_sample(1 + dfvara_, varg);
                    rhs = r_hat[i];
                    if(gi){rhs += xx * gi;}
                    lhs = xx / varei;
                    v = xx + varei / varg;
                    gi = norm_sample(rhs / v, sqrt(varei / v));
                    gi_ = (g[i] - gi) * n;
                    start = ldm.begin_col(i);
	    			end = ldm.end_col(i);
	    			for (; start != end; ++start){
	    				r_hat[start.row()] += gi_ * (*start);
	   				}
                    g[i] = gi;
                }
                break;
            case 3:
                logpi = log(Pi);
                s[0] = logpi[0];
                for(int i = 0; i < m; i++){
                    if(!ifest[i])   continue;
                    xx = xpx[i];
                    gi = g[i];
                    varei = varediff[i] * vara_ + vare_;
                    varg = (gi * gi + s2varg_ * dfvara_) / (dfvara_ + 1);
                    varg = invchisq_sample(1 + dfvara_, varg);
                    rhs = r_hat[i];
                    if(gi){rhs += xx * gi;}
                    lhs = xx / varei;
                    logdetV = log(varg * lhs + 1);
                    uhat = rhs / (xx + varei / varg);
                    s[1] = -0.5 * (logdetV - (rhs * uhat / varei)) + logpi[1];
                    acceptProb = 1 / sum(exp(s - s[0]));
                    indistflag = (uniform_sample()) < acceptProb ? 0 : 1;
                    snptracker[i] = indistflag;
                    if(indistflag == 0){
                        gi = 0;
                    }else{
                        v = xx + varei / varg;
                        gi = norm_sample(rhs / v, sqrt(varei / v));
                    }
                    if(gi != g[i]){
                        gi_ = (g[i] - gi) * n;
                        start = ldm.begin_col(i);
                        end = ldm.end_col(i);
                        for (; start != end; ++start){
                            r_hat[start.row()] += gi_ * (*start);
                        }
                        g[i] = gi;
                    }
                }
                fold_snp_num[1] = sum(snptracker);
                fold_snp_num[0] = m - fold_snp_num[1];
                NnzSnp = fold_snp_num[1];
                if(!fixpi)  Pi = rdirichlet_sample(n_fold, (fold_snp_num + 1));
                break;
            case 4:
                logpi = log(Pi);
                s[0] = logpi[0];
                vargi = 0;
                for(int i = 0; i < m; i++){
                    if(!ifest[i])   continue;
                    xx = xpx[i];
                    gi = g[i];
                    varei = varediff[i] * vara_ + vare_;
                    rhs = r_hat[i];
                    if(gi){rhs += xx * gi;}
                    lhs = xx / varei;
                    logdetV = log(varg * lhs + 1);
                    uhat = rhs / (xx + varei / varg);
                    s[1] = -0.5 * (logdetV - (rhs * uhat / varei)) + logpi[1];
                    acceptProb = 1 / sum(exp(s - s[0]));	
                    indistflag = (uniform_sample()) < acceptProb ? 0 : 1;
                    snptracker[i] = indistflag;

                    if(indistflag == 0){
                        gi = 0;
                    }else{
                        v = xx + varei / varg;
                        gi = norm_sample(rhs / v, sqrt(varei / v));
                        if((gi * gi * vx[i]) > vary){
                            int ii = 0;
                            while((gi * gi * vx[i]) > vary){
                                gi = norm_sample(rhs / v, sqrt(varei / v));
                                vargi = gi * gi;
                                ii++;
                                if(ii > 100){
                                    gi = 0;
                                }
                            }
                        }
                        vargi += gi * gi;
                    }
                    if(gi != g[i]){
                        gi_ = (g[i] - gi) * n;
                        start = ldm.begin_col(i);
                        end = ldm.end_col(i);
                        for (; start != end; ++start){
                            r_hat[start.row()] += gi_ * (*start);
                        }
                        g[i] = gi;
                    }
                }
                fold_snp_num[1] = sum(snptracker);
                fold_snp_num[0] = m - fold_snp_num[1];
                NnzSnp = fold_snp_num[1];
                varg = (vargi + s2varg_ * dfvara_) / (dfvara_ + NnzSnp);
                varg = invchisq_sample(NnzSnp + dfvara_, varg);
                if(!fixpi)  Pi = rdirichlet_sample(n_fold, (fold_snp_num + 1));
                break;
            case 5:
                for(int i = 0; i < m; i++){
                    if(!ifest[i])   continue;
                    xx = xpx[i];
                    gi = g[i];
                    varei = varediff[i] * vara_ + vare_;
                    rhs = r_hat[i];
                    if(gi){rhs += xx * gi;}
                    lhs = xx / varei;
                    v = xx + 1 / vargL[i];
                    gi = norm_sample(rhs / v, sqrt(varei / v));
                    if(abs(gi) < 1e-6){gi = 1e-6;}
                    vargi = 1 / rinvgaussian_sample(sqrt(varei) * lambda / abs(gi), lambda2);
                    if(vargi > 0)   vargL[i] = vargi;
                    if(gi != g[i]){
                        gi_ = (g[i] - gi) * n;
                        start = ldm.begin_col(i);
                        end = ldm.end_col(i);
                        for (; start != end; ++start){
                            r_hat[start.row()] += gi_ * (*start);
                        }
                        g[i] = gi;
                    }
                }
                shape = shape0 + count_y;
                rate = rate0 + sum(vargL) / 2;
                lambda2 = gamma_sample(shape, 1 / rate);
                lambda = sqrt(lambda2);
                break;
            case 6:
                logpi = log(Pi);
                s[0] = logpi[0];
                // sumvg = 0;
                varg = 0;
                for(int i = 0; i < m; i++){
                    if(!ifest[i])   continue;
                    xx = xpx[i];
			        gi = g[i];
                    varei = varediff[i] * vara_ + vare_;
                    rhs = r_hat[i];
                    if(gi){rhs += xx * gi;}
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
                    acceptProb = 0;
                    indistflag = 0;
                    double rval = uniform_sample();

                    for(int j = 0; j < n_fold; j++){
                        acceptProb += stemp[j];
                        if(rval < acceptProb){
                            indistflag = j;
                            break;
                        }
                    }
                    snptracker[i] = indistflag;
                    if(indistflag == 0){
                        gi = 0;
                    }else{
                        v = xx + varei / vara_fold[indistflag];
                        gi = norm_sample(rhs / v, sqrt(varei / v));
                        if((gi * gi * vx[i]) > vary){
                            int ii = 0;
                            while((gi * gi * vx[i]) > vary){
                                gi = norm_sample(rhs / v, sqrt(varei / v));
                                vargi = gi * gi;
                                ii++;
                                if(ii > 100){
                                    gi = 0;
                                }
                            }
                        }
                        varg += (gi * gi / fold_[indistflag]);
                    }
                    if(gi != g[i]){
                        gi_ = (g[i] - gi) * n;
                        start = ldm.begin_col(i);
                        end = ldm.end_col(i);
                        for (; start != end; ++start){
                            r_hat[start.row()] += gi_ * (*start);
                        }
                        g[i] = gi;
                    }
                }
                for(int j = 0; j < n_fold; j++){
                    fold_snp_num[j] = sum(snptracker == j);
                }
                NnzSnp = m - fold_snp_num[0];
                varg = (varg + s2varg_ * dfvara_) / (dfvara_ + NnzSnp);
                varg = invchisq_sample(NnzSnp + dfvara_, varg);
                for(int j = 0; j < n_fold; j++){
                    vara_fold[j] = varg * fold_[j]; 
                    // vara_fold[j] = vara_ * fold_[j]; 
                }
                if(!fixpi)  Pi = rdirichlet_sample(n_fold, fold_snp_num + 1);
                break;
        }

        // genetic variance
        // vara_ = sum(g * (xy - r_hat)) / (n - 1);
        tmp = (xy - r_hat);
		// vara_ = ddot_(&m, g.memptr(), &inc, tmp.memptr(), &inc) / (n - 1);
        vara_ = (ddot_(&m, g.memptr(), &inc, tmp.memptr(), &inc) + s2vara_ * dfvara_) / (n + dfvara_);
        vara_ = invchisq_sample(n + dfvara_, vara_);
    
        // sample residual variance from inv-chisq distribution
        // vare_ = (yy - sum(g * (xy + r_hat)) + s2vare_ * dfvare_) / (n + dfvare_);
        tmp = (xy + r_hat);
		vare_ = (yy - ddot_(&m, g.memptr(), &inc, tmp.memptr(), &inc) + s2vare_ * dfvare_) / (n + dfvare_);
		if(vare_ < 0)	vare_ = 0;
		vare_ = invchisq_sample(n + dfvare_, vare_);

        if(iter >= nburn){
            
            if(!snptracker.is_empty()){
                #pragma omp parallel for
                for(int i = 0; i < m; i++){
                    if(snptracker[i]){
                        nzrate[i] += 1;
                    }
                }
            }
            
            if(WPPA){
                for(int w = 0; w < nw; w++){
                    windxi = windx[w];
                    if(any(snptracker.elem(windxi))){
                        wppai[w] += 1;
                    }
                }
            }
            
            nzct++;
        }

        if(iter >= nburn && (iter + 1 - nburn) % outfreq == 0){
            
            if(!fixpi)  pi_store.col(count) = Pi;
            vara_store[count] = vara_;
            vare_store[count] = vare_;
            
            // daxpy_(&m, &doc, g.memptr(), &inc, g_store.memptr(), &inc);
            g_store.col(count) = g;

            hsq_store[count] = vara_ / (vara_ + vare_);

            count++;
        }
        
        // print iteration details
        if(verbose && (iter + 1) % outfreq == 0){
            timer_loop.step("end");
            NumericVector tdiff(timer_loop);
            tt0 = (tdiff[1] - tdiff[0]) / (iter + 1) / 1e9;
            tt = floor(tt0 * (niter - iter));
            hor = floor(tt / 3600);
            min = floor(tt % 3600 / 60);
            sec = floor(tt % 3600 % 60);
            Rcpp::Rcout << " " << iter + 1 << " ";
            Rcpp::Rcout << NnzSnp << " ";
            for(int i = 0; i < n_fold; i++){
                Rcpp::Rcout << std::fixed << Pi[i] << " ";
            }
            if(model == "BayesL")    Rcpp::Rcout << std::fixed << lambda << " ";
            Rcpp::Rcout << std::fixed << vara_ << " ";
            Rcpp::Rcout << std::fixed << vare_ << " ";
            Rcpp::Rcout << std::fixed << (vara_ / (vara_ + vare_)) << " ";
            Rprintf("%02dh%02dm%02ds \n", hor, min, sec);
        }

        if(count == n_records)  break;
    }

    List results;
    List MCMCsample;

    vara_ = mean(vara_store);
    double varasd = stddev(vara_store);
    vare_ = mean(vare_store);
    double varesd = stddev(vare_store);
    hsq = mean(hsq_store);
    double hsqsd = stddev(hsq_store);
    results["Vg"] = vara_;
    results["Ve"] = vare_;
    results["h2"] = hsq;
    MCMCsample["Vg"] = vara_store.t();
    MCMCsample["Ve"] = vare_store.t();
    MCMCsample["h2"] = hsq_store.t();
    
    g = conv_to<vec>::from(mean(g_store, 1));
    results["alpha"] = g;
    MCMCsample["alpha"] = g_store;

    vec pisd;
    if(!fixpi){
        Pi = conv_to<vec>::from(mean(pi_store, 1));
        pisd = conv_to<vec>::from(stddev(pi_store, 0, 1));
    }else{
        pisd = zeros(Pi.n_elem);
        pi_store.row(0).fill(Pi[0]);
        pi_store.row(1).fill(Pi[1]);
    }
    results["pi"] = Pi;
    MCMCsample["pi"] = pi_store;

    if(nzrate.is_empty()){
        nzrate.ones(m);
    }else{
        nzrate = nzrate / nzct;
        nzrate.elem(find(nzrate == 1)).fill((nzct - 1) / (double) nzct);
    }
    results["pip"] = nzrate;

    if(WPPA){
        wppai = wppai / nzct;
        wppai.elem(find(wppai == 1)).fill((nzct - 1) / (double) nzct);
        results["gwas"] = wppai;
    }

    results["MCMCsamples"] = MCMCsample;

    if(verbose){
        Rcpp::Rcout << "Posterior parameters:" << std::endl;
        for(int i = 0; i < Pi.n_elem; i++){
            if(i == 0){
                Rcpp::Rcout << "    π for markers in zero effect size " << std::fixed << Pi[i] << "±" << std::fixed << pisd[i] << endl;
            }else{
                if(i == 1){
                    Rcpp::Rcout << "    π for markers in non-zero effect size " << std::fixed << Pi[i] << "±" << std::fixed << pisd[i];
                }else{
                    Rcpp::Rcout << std::fixed << Pi[i] << "±" << std::fixed << pisd[i];
                }
                if(i != (Pi.n_elem - 1))  Rcpp::Rcout << ", ";
            }
        }
        Rcpp::Rcout << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << "±" << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << "±" << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated h2 " << std::fixed << hsq << "±" << std::fixed << hsqsd << std::endl;
    }

    timer.step("end");
    NumericVector res(timer);
    tt = floor((res[1] - res[0]) / 1e9);
    hor = floor(tt / 3600);
    min = floor(tt % 3600 / 60);
    sec = floor(tt % 3600 % 60);
    if(verbose){Rprintf("Finished within total run time: %02dh%02dm%02ds \n", hor, min, sec);}
    
    return results;
}
