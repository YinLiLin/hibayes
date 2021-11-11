#include "hibayes.h"
#include "solver.h"

template <class T>
bool xhasNA(T &X){
    bool hasna = false;
    for(int i = 0; i < X.ncol(); i++){
        for(int j = 0; j < X.nrow(); j++){
            if(T::is_na(X(j, i))){
                hasna = true;
                break;
            }
        }
    }
    return hasna;
}

bool yhasNA(NumericVector y){
    bool hasna = false;
    for(int i = 0; i < y.length(); i++){
        if(NumericVector::is_na(y[i])){
            hasna = true;
            break;
        }
    }
    return hasna;
}

Rcpp::List makeZ(
    const CharacterVector &R
)
{
    int n = R.length();
    std::vector<std::string> levels;
    std::vector<std::string> value = Rcpp::as<std::vector<std::string> >(R);
    stable_sort(value.begin(), value.end());
    value.erase(unique(value.begin(), value.end()), value.end());
    if (value.size() == n){
        throw Rcpp::exception("number of class of environmental random effects should be less than population size.");
    }
    if(value.size() == 1){
        throw Rcpp::exception("number of class of environmental random effects should be bigger than 1.");
    }
    std::map<string, int> val_map;
    for (int j = 0; j < value.size(); j++){
        levels.push_back(value[j]);
        val_map.insert(pair<string, int>(value[j], j));
    }
    map<string, int>::iterator iter;
    arma::sp_mat z(n, value.size());
    for (int j = 0; j < n; j++) {
        std::string v = Rcpp::as<string>(R[j]);
        iter = val_map.find(v);
        z(j, iter->second) = 1.0;
    }
    return List::create(Named("z") = z, Named("levels") = levels);
}

// --------------------------------------- BAYES ---------------------------------------

// [[Rcpp::export]]
Rcpp::List BayesRR(
    const NumericVector &y,
    const NumericMatrix &X,
    const Nullable<NumericMatrix> C = R_NilValue,
    const Nullable<CharacterMatrix> R = R_NilValue,
    const int niter = 50000,
    const int nburn = 20000,
    const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
    const int outfreq = 100,
    const bool verbose = true
){
    if(yhasNA(y)){
        throw Rcpp::exception("NAs are not allowed in y.");
    }
    // if(xhasNA(X)){
    //     throw Rcpp::exception("NAs are not allowed in X.");
    // }
    if(y.length() != X.nrow()){
        throw Rcpp::exception("Number of individuals not equals.");
    }
    double* dci;
    double* pC;
    int nc = 0;
    NumericMatrix C_;
    if(C.isNotNull()){
        C_ = as<NumericMatrix>(C);
        if(C_.nrow() != X.nrow()){
            throw Rcpp::exception("Number of individuals not match for covariates.");
        }
        if(xhasNA(C_)){
            throw Rcpp::exception("NAs are not allowed in covariates.");
        }
        pC = NUMERIC_POINTER(C_);
        nc = C_.ncol();
    }
    NumericVector beta(nc);
    NumericVector beta_store(nc);
    NumericVector cpc(nc);
    if(nc){
        for(int i = 0; i < nc; i++){
            NumericVector ci = C_(_, i);
            cpc[i] = sum(ci * ci);
        }
    }
    
    int nr = 0;
    vec vr;
    vec vrsd;
    mat vr_store;
    vector<sp_mat> Z;
    vector<sp_mat> ZZ;
    vector<vec> estR;
    vector<vec> estR_store;
    vector< vector<string> > Z_levels;
    if(R.isNotNull()){
        CharacterMatrix R_ = as<CharacterMatrix>(R);
        if(R_.nrow() != X.nrow()){
            throw Rcpp::exception("Number of individuals not match for environmental random effects.");
        }
        if(xhasNA(R_)){
            throw Rcpp::exception("NAs are not allowed in environmental random effects.");
        }
        nr = R_.ncol();
        vr.resize(nr);
        vr_store.resize(niter - nburn + 1, nr);
        vr.fill(var(y) / nr);
        for(int i = 0; i < nr; i++){
            CharacterVector Ri = R_(_, i);
            List Z_info = makeZ(Ri);
            sp_mat z = Z_info[0];
            vector<string> z_levels = Z_info[1];
            Z.push_back(z);
            ZZ.push_back(z.t() * z);
            Z_levels.push_back(z_levels);
            vec estRi = zeros(z_levels.size());
            estR.push_back(estRi);
            estR_store.push_back(estRi);
        }
    }

    int n = y.length();
    int m = X.ncol();
    bool WPPA = false;
    int count = 0;
    int inc = 1;
    double doc = 1.0;
    double xx, oldgi, gi, gi_, rhs, v;
    double vara_, dfvara_, s2vara_, vare_, dfvare_, s2vare_, sumvg, hsq;
    NumericVector g(m);
    NumericVector g_store(m);
    NumericVector xi(n);
    NumericVector u(n);
    double* du = NUMERIC_POINTER(u);
    NumericVector xpx(m), vx(m);
    NumericVector mu_store(niter - nburn + 1), sumvg_store(niter - nburn + 1), vara_store(niter - nburn + 1), vare_store(niter - nburn + 1), hsq_store(niter - nburn + 1);
    for(int i = 0; i < m; i++){
        xi = X(_, i);
        xpx[i] = sum(xi * xi);
        vx[i] = var(xi);
    }
    double vary = var(y);
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
        vara_ = ((dfvara_ - 2) / dfvara_) * vary;
    }
    if(ve.isNotNull()){
        vare_ = as<double>(ve);
    }else{
        vare_ = vary;
    }
    if(dfve.isNotNull()){
        dfvare_ = as<double>(dfve);
    }else{
        dfvare_ = -2;
    }
    if(s2vg.isNotNull()){
        s2vara_ = as<double>(s2vg);
    }else{
        s2vara_ = (vara_ / (sum(vx))) * (dfvara_ - 2) / dfvara_;
    }
    if(s2ve.isNotNull()){
        s2vare_ = as<double>(s2ve);
    }else{
        s2vare_ = 0;
    }
    if(niter < nburn){
        throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
    }

    // for gwas
    R_xlen_t nw;
    int ni;
    double varw;
    NumericVector windindx_(m);
    Rcpp::List windx(m);
    NumericVector wu(n);
    double* dwu = NUMERIC_POINTER(wu);
    IntegerVector windxi;
    NumericVector wppai(m);
    NumericVector wgvei(m);
    if(windindx.isNotNull()){
        windindx_ = as<IntegerVector>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = wppai[seq(0, (nw - 1))];
        wgvei = wgvei[seq(0, (nw - 1))];
        for(R_xlen_t w = 0; w < nw; w++){
            windx[w] = which_c(windindx_, (w+1), 5);
        }
    }

    double varg = (vara_ / (sum(vx)));

    if(verbose){
        Rcpp::Rcout.precision(4);
        Rcpp::Rcout << "Prior parameters:" << std::endl;
        Rcpp::Rcout << "    Model fitted at [Bayes Ridge Regression]" << std::endl;
        Rcpp::Rcout << "    Number of observations " << n << std::endl;
        Rcpp::Rcout << "    Number of covariates " << nc << std::endl;
        Rcpp::Rcout << "    Number of envir-random effects " << nr << std::endl;
        Rcpp::Rcout << "    Number of markers " << m << std::endl;
        Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
        Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
        if(nr){
            Rcpp::Rcout << "    Environmental var ";
            for(int i = 0; i < nr; i++){
                Rcpp::Rcout << std::fixed << vr[i] << " ";
            }
            Rcpp::Rcout << std::endl;
        }
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << std::fixed << s2vara_ << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << std::fixed << s2vare_ << std::endl;
        Rcpp::Rcout << "    Phenotypic var " << std::fixed << vary << std::endl;
        if(WPPA){
            Rcpp::Rcout << "    Number of windows " << nw << std::endl;
            Rcpp::Rcout << "    GVE threshold for windows " << wppa << std::endl;
        }
        Rcpp::Rcout << "MCMC started: " << std::endl;
        Rcpp::Rcout << " Iter" << "\t";
        for(int i = 0; i < nr; i++){
            Rcpp::Rcout << "Vr" << (i + 1) << "\t";
        }
        Rcpp::Rcout << "SigmaSq" << "\t";
        Rcpp::Rcout << "Vg" << "\t";
        Rcpp::Rcout << "Ve" << "\t";
        Rcpp::Rcout << "h2" << "\t";
        Rcpp::Rcout << "Timeleft" << std::endl;
    }

    MyTimer timer;
    timer.step("start");
    MyTimer timer_loop;
    timer_loop.step("start");
    double tt0; 
    int tt, hor, min, sec;
    double mu_, mu = mean(y);
    NumericVector one(n); one.fill(1);
    double* done = NUMERIC_POINTER(one);
    NumericVector yadj = y - mu;
    double* dyadj = NUMERIC_POINTER(yadj);
    double* dxi;
    double* pX = NUMERIC_POINTER(X);

    // MCMC procedure
    for(int iter = 0; iter < niter; iter++){

        // sample intercept
        mu_ = - norm_sample(mean(yadj), sqrt(vare_ / n));
        mu -= mu_;
		daxpy_(&n, &mu_, done, &inc, dyadj, &inc);

        for(int i = 0; i < nc; i++){
            dci = pC+(long long)i*n;
			oldgi = beta[i];
            v = cpc[i];
			rhs = ddot_(&n, dci, &inc, dyadj, &inc);
			rhs += v * oldgi;
			gi = norm_sample(rhs / v, sqrt(vare_ / v));
			gi_ = oldgi - gi;
			daxpy_(&n, &gi_, dci, &inc, dyadj, &inc);
            beta[i] = gi;
        }

        for(int i = 0; i < nr; i++){
            sp_mat LHS = ZZ[i];
            LHS.diag() += vare_ / vr[i];
            vec arma_yadj(yadj);
            vec estR_tmp(estR[i]);
            vec RHS = Z[i].t() * arma_yadj;
            Gibbs(LHS, estR_tmp, RHS);
            estR[i] -= estR_tmp;
            vec diff = Z[i].t() * estR[i];
            daxpy_(&n, &doc, diff.memptr(), &inc, dyadj, &inc);
            int q = estR_tmp.n_elem;
            vr[i] = (ddot_(&q, estR_tmp.memptr(), &inc, estR_tmp.memptr(), &inc) + s2vara_ * dfvara_) / (dfvara_ + q);
            vr[i] = invchisq_sample(q + dfvara_, vr[i]);
            estR[i] = estR_tmp;
        }

        // loop on snps
        for(int i = 0; i < m; i++){

            dxi = pX+(long long)i*n;
            // xi = X[i];   // for list
            // dxi = NUMERIC_POINTER(xi);
            
			xx = xpx[i];
			oldgi = g[i];

			// right hand
			rhs = ddot_(&n, dxi, &inc, dyadj, &inc);
			rhs += xx * oldgi;

			// left hand
			v = xx + vare_ / varg;

			// sample snp effect
			gi = norm_sample(rhs / v, sqrt(vare_ / v));

            // update residual
			gi_ = oldgi - gi;
			daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);

            // update random effect
            gi_ *= -1;
            daxpy_(&n, &gi_, dxi, &inc, du, &inc);

			// update snp effect
			g[i] = gi;

        }

        // sample variance
        varg = (ddot_(&m, g.begin(), &inc, g.begin(), &inc) + s2vara_ * dfvara_) / (dfvara_ + m);
        varg = invchisq_sample(m + dfvara_, varg);

        sumvg = sum(vx * varg);

        // genetic variance
        vara_ = var(u);
    
        // sample residual variance from inv-chisq distribution
        vare_ = (ddot_(&n, dyadj, &inc, dyadj, &inc) + s2vare_ * dfvare_) / (n + dfvare_);
        vare_ = invchisq_sample(n + dfvare_, vare_);

        if(iter >= nburn){
            count++;
            mu_store[iter - nburn] = mu;
            sumvg_store[iter - nburn] = sumvg;
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            hsq_store[iter - nburn] = vara_ / (vara_ + vare_);
            daxpy_(&m, &doc, g.begin(), &inc, g_store.begin(), &inc);
            // g_store += g;
            if(nc) daxpy_(&nc, &doc, beta.begin(), &inc, beta_store.begin(), &inc);
            if(nr){
                vr_store.row(iter - nburn) = vr;
                for(int i = 0; i < nr; i++){
                    estR_store[i] += estR[i];
                }
            }

            // record pve for each window
            if(WPPA){
                for(R_xlen_t w = 0; w < nw; w++){
                    wu.fill(0);
                    windxi = windx[w];
                    for(int i = 0; i < windxi.size(); i++){
                        ni = windxi[i];
                        gi_ = g[ni];
                        if(gi_){
                            dxi = pX+(long long)ni*n;
                            daxpy_(&n, &gi_, dxi, &inc, dwu, &inc);
                        }
                        // if(g[windxi[i]]){wu += X(_, windxi[i]) * g[windxi[i]];}
                    }
                    varw = var(wu);
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
                for(int i = 0; i < nr; i++){
                    Rcpp::Rcout << std::fixed << vr[i] << " ";
                }
                Rcpp::Rcout << std::fixed << sumvg << " ";
                Rcpp::Rcout << std::fixed << vara_ << " ";
                Rcpp::Rcout << std::fixed << vare_ << " ";
                Rcpp::Rcout << std::fixed << (vara_ / (vara_ + vare_)) << " ";
                Rprintf("%02dh%02dm%02ds \n", hor, min, sec);
            }
        }
    }

    g = g_store / count;
    if(nc) beta = beta_store / count;
    DataFrame r;
    if(nr){
        vr = mean(vr_store);
        vrsd = stddev(vr_store);
        vector<string> r_levers;
        vector<double> estr;
        for(int i = 0; i < nr; i++){
            for(int j = 0; j < (estR_store[i]).n_elem; j++){
                r_levers.push_back(Z_levels[i][j]);
                estr.push_back(estR_store[i][j] / count);
            }
        }
        r = DataFrame::create(r_levers, estr);
    }
    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    double Mu = mean(mu_store);
    double musd = sd(mu_store);
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
        Rcpp::Rcout << "    Mu " << std::fixed << Mu << " ± " << std::fixed << musd << std::endl;
        if(nr){
            Rcpp::Rcout << "    Environmental var ";
            for(int i = 0; i < nr; i++){
                Rcpp::Rcout << std::fixed << vr[i] << " ± " << std::fixed << vrsd[i];
                if(i != (nr - 1))  Rcpp::Rcout << ", ";
            }
            Rcpp::Rcout << std::endl;
        }
        Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << " ± " << std::fixed << sumvgsd << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << " ± " << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << " ± " << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated hsq " << std::fixed << hsq << " ± " << std::fixed << hsqsd << std::endl;
        if(nc){
            Rcpp::Rcout << "    Estimated effects for covariates: ";
            for(int i = 0; i < nc; i++){
                Rcpp::Rcout << std::fixed << beta[i] << " ";
            }
            Rcpp::Rcout << std::endl;
        }
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
            Named("mu") = Mu, 
            Named("beta") = beta,
            Named("r") = r, 
            Named("vr") = vr, 
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("g") = g,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("mu") = Mu,
            Named("beta") = beta,
            Named("r") = r, 
            Named("vr") = vr, 
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("g") = g
        );
    }
}

// [[Rcpp::export]]
Rcpp::List BayesA(
    const NumericVector &y,
    const NumericMatrix &X,
    const Nullable<NumericMatrix> C = R_NilValue,
    const int niter = 50000,
    const int nburn = 20000,
    const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
    const int outfreq = 100,
    const bool verbose = true
){
    if(yhasNA(y)){
        throw Rcpp::exception("NAs are not allowed in y.");
    }
    // if(xhasNA(X)){
    //     throw Rcpp::exception("NAs are not allowed in X.");
    // }
    if(y.length() != X.nrow()){
        throw Rcpp::exception("Number of individuals not equals.");
    }
    double* dci;
    double* pC;
    int nc = 0;
    NumericMatrix C_;
    if(C.isNotNull()){
        C_ = as<NumericMatrix>(C);
        if(C_.nrow() != X.nrow()){
            throw Rcpp::exception("Number of individuals not match for covariates.");
        }
        if(xhasNA(C_)){
            throw Rcpp::exception("NAs are not allowed in covariates.");
        }
        pC = NUMERIC_POINTER(C_);
        nc = C_.ncol();
    }
    NumericVector beta(nc);
    NumericVector beta_store(nc);
    NumericVector cpc(nc);
    if(nc){
        for(int i = 0; i < nc; i++){
            NumericVector ci = C_(_, i);
            cpc[i] = sum(ci * ci);
        }
    }

    int nr = 0;

    int n = y.length();
    int m = X.ncol();
    bool WPPA = false;
    int count = 0;
    int inc = 1;
    double doc = 1.0;
    double xx, oldgi, gi, gi_, rhs, v;
    double vara_, dfvara_, s2vara_, vare_, dfvare_, s2vare_, sumvg, hsq;
    NumericVector g(m);
    NumericVector g_store(m);
    NumericVector xi(n);
    NumericVector u(n);
    double* du = NUMERIC_POINTER(u);
    NumericVector xpx(m), vx(m);
    NumericVector mu_store(niter - nburn + 1), sumvg_store(niter - nburn + 1), vara_store(niter - nburn + 1), vare_store(niter - nburn + 1), hsq_store(niter - nburn + 1);
    for(int i = 0; i < m; i++){
        xi = X(_, i);
        xpx[i] = sum(xi * xi);
        vx[i] = var(xi);
    }
    double vary = var(y);
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
        vara_ = ((dfvara_ - 2) / dfvara_) * vary;
    }
    if(ve.isNotNull()){
        vare_ = as<double>(ve);
    }else{
        vare_ = vary;
    }
    if(dfve.isNotNull()){
        dfvare_ = as<double>(dfve);
    }else{
        dfvare_ = -2;
    }
    if(s2vg.isNotNull()){
        s2vara_ = as<double>(s2vg);
    }else{
        s2vara_ = (vara_ / (sum(vx))) * (dfvara_ - 2) / dfvara_;
    }
    if(s2ve.isNotNull()){
        s2vare_ = as<double>(s2ve);
    }else{
        s2vare_ = 0;
    }
    if(niter < nburn){
        throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
    }

    // for gwas
    R_xlen_t nw;
    int ni;
    double varw;
    NumericVector windindx_(m);
    Rcpp::List windx(m);
    NumericVector wu(n);
    double* dwu = NUMERIC_POINTER(wu);
    IntegerVector windxi;
    NumericVector wppai(m);
    NumericVector wgvei(m);
    if(windindx.isNotNull()){
        windindx_ = as<IntegerVector>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = wppai[seq(0, (nw - 1))];
        wgvei = wgvei[seq(0, (nw - 1))];
        for(R_xlen_t w = 0; w < nw; w++){
            windx[w] = which_c(windindx_, (w+1), 5);
        }
    }

    double varg = (vara_ / (sum(vx)));

    if(verbose){
        Rcpp::Rcout.precision(4);
        Rcpp::Rcout << "Prior parameters:" << std::endl;
        Rcpp::Rcout << "    Model fitted at [BayesA]" << std::endl;
        Rcpp::Rcout << "    Number of observations " << n << std::endl;
        Rcpp::Rcout << "    Number of covariates " << nc << std::endl;
        Rcpp::Rcout << "    Number of envir-random effects " << nr << std::endl;
        Rcpp::Rcout << "    Number of markers " << m << std::endl;
        Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
        Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << std::fixed << s2vara_ << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << std::fixed << s2vare_ << std::endl;
        Rcpp::Rcout << "    Phenotypic var " << std::fixed << vary << std::endl;
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
    double mu_, mu = mean(y);
    NumericVector one(n); one.fill(1);
    double* done = NUMERIC_POINTER(one);
    NumericVector yadj = y - mu;
    double* dyadj = NUMERIC_POINTER(yadj);
    double* dxi;
    double* pX = NUMERIC_POINTER(X);

    // MCMC procedure
    for(int iter = 0; iter < niter; iter++){

        // sample intercept
        mu_ = - norm_sample(mean(yadj), sqrt(vare_ / n));
        mu -= mu_;
        daxpy_(&n, &mu_, done, &inc, dyadj, &inc);
        sumvg = 0;

        for(int i = 0; i < nc; i++){
            dci = pC+(long long)i*n;
			oldgi = beta[i];
            v = cpc[i];
			rhs = ddot_(&n, dci, &inc, dyadj, &inc);
			rhs += v * oldgi;
			gi = norm_sample(rhs / v, sqrt(vare_ / v));
			gi_ = oldgi - gi;
			daxpy_(&n, &gi_, dci, &inc, dyadj, &inc);
            beta[i] = gi;
        }

        // loop on snps
        for(int i = 0; i < m; i++){

            dxi = pX+(long long)i*n;
            // xi = X[i];   // for list
            // dxi = NUMERIC_POINTER(xi);
            
            xx = xpx[i];
            oldgi = g[i];

            // sample variance
            varg = (oldgi * oldgi + s2vara_ * dfvara_) / (dfvara_ + 1);
            varg = invchisq_sample(1 + dfvara_, varg);

            // right hand
            rhs = ddot_(&n, dxi, &inc, dyadj, &inc);
            rhs += xx * oldgi;

            // left hand
            v = xx + vare_ / varg;

            // sample snp effect
            gi = norm_sample(rhs / v, sqrt(vare_ / v));

            // update residual
            gi_ = oldgi - gi;
            daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);

            // update random effect
            gi_ *= -1;
            daxpy_(&n, &gi_, dxi, &inc, du, &inc);

            // update snp effect
            g[i] = gi;
            sumvg += (gi * gi) * vx[i];

        }

        // genetic variance
        vara_ = var(u);
    
        // sample residual variance from inv-chisq distribution
        vare_ = (ddot_(&n, dyadj, &inc, dyadj, &inc) + s2vare_ * dfvare_) / (n + dfvare_);
        vare_ = invchisq_sample(n + dfvare_, vare_);

        if(iter >= nburn){
            count++;
            mu_store[iter - nburn] = mu;
            sumvg_store[iter - nburn] = sumvg;
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            hsq_store[iter - nburn] = vara_ / (vara_ + vare_);
            daxpy_(&m, &doc, g.begin(), &inc, g_store.begin(), &inc);
            // g_store += g;
            if(nc) daxpy_(&nc, &doc, beta.begin(), &inc, beta_store.begin(), &inc);
            
            // record pve for each window
            if(WPPA){
                for(R_xlen_t w = 0; w < nw; w++){
                    wu.fill(0);
                    windxi = windx[w];
                    for(int i = 0; i < windxi.size(); i++){
                        ni = windxi[i];
                        gi_ = g[ni];
                        if(gi_){
                            dxi = pX+(long long)ni*n;
                            daxpy_(&n, &gi_, dxi, &inc, dwu, &inc);
                        }
                        // if(g[windxi[i]]){wu += X(_, windxi[i]) * g[windxi[i]];}
                    }
                    varw = var(wu);
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

    g = g_store / count;
    beta = beta_store / count;
    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    double Mu = mean(mu_store);
    double musd = sd(mu_store);
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
        Rcpp::Rcout << "    Mu " << std::fixed << Mu << " ± " << std::fixed << musd << std::endl;
        Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << " ± " << std::fixed << sumvgsd << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << " ± " << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << " ± " << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated hsq " << std::fixed << hsq << " ± " << std::fixed << hsqsd << std::endl;
        if(nc){
            Rcpp::Rcout << "    Estimated effects for covariates: ";
            for(int i = 0; i < nc; i++){
                Rcpp::Rcout << std::fixed << beta[i] << " ";
            }
            Rcpp::Rcout << std::endl;
        }
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
            Named("mu") = Mu,
            Named("beta") = beta, 
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("g") = g,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("mu") = Mu, 
            Named("beta") = beta,
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("g") = g
        );
    }
}

// [[Rcpp::export]]
Rcpp::List BayesBpi(
    const NumericVector &y,
    const NumericMatrix &X,
    const Nullable<NumericMatrix> C = R_NilValue,
    const double pi = 0.95,
    const int niter = 50000,
    const int nburn = 20000,
    const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
    const int outfreq = 100,
    const bool fixpi = false,
    const bool verbose = true
){
    if(yhasNA(y)){
        throw Rcpp::exception("NAs are not allowed in y.");
    }
    // if(xhasNA(X)){
    //     throw Rcpp::exception("NAs are not allowed in X.");
    // }
    if(y.length() != X.nrow()){
        throw Rcpp::exception("Number of individuals not equals.");
    }
    double* dci;
    double* pC;
    int nc = 0;
    NumericMatrix C_;
    if(C.isNotNull()){
        C_ = as<NumericMatrix>(C);
        if(C_.nrow() != X.nrow()){
            throw Rcpp::exception("Number of individuals not match for covariates.");
        }
        if(xhasNA(C_)){
            throw Rcpp::exception("NAs are not allowed in covariates.");
        }
        pC = NUMERIC_POINTER(C_);
        nc = C_.ncol();
    }
    NumericVector beta(nc);
    NumericVector beta_store(nc);
    NumericVector cpc(nc);
    if(nc){
        for(int i = 0; i < nc; i++){
            NumericVector ci = C_(_, i);
            cpc[i] = sum(ci * ci);
        }
    }

    int nr = 0;

    int n = y.length();
    int m = X.ncol();
    bool WPPA = false;
    int NnzSnp, indistflag, count = 0;
    int inc = 1;
    double doc = 1.0;
    double xx, oldgi, gi, gi_, rhs, lhs, logdetV, acceptProb, uhat, r, v;
    double vara_, dfvara_, s2vara_, vare_, dfvare_, s2vare_, sumvg, hsq;
    NumericVector pi_;
    NumericVector snptracker(m);
    NumericVector nzrate(m);
    NumericVector g(m);
    NumericVector g_store(m);
    NumericVector xi(n);
    NumericVector u(n);
    double* du = NUMERIC_POINTER(u);
    NumericVector xpx(m), vx(m);
    NumericVector mu_store(niter - nburn + 1), sumvg_store(niter - nburn + 1), vara_store(niter - nburn + 1), vare_store(niter - nburn + 1), hsq_store(niter - nburn + 1);
    for(int i = 0; i < m; i++){
        xi = X(_, i);
        xpx[i] = sum(xi * xi);
        vx[i] = var(xi);
    }
    double vary = var(y);
    if(pi <= 0 || pi >= 1){
        throw Rcpp::exception("pi should be at 0 < pi < 1.");
    }else{
        pi_ = {pi, 1 - pi};
    }
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
        vara_ = ((dfvara_ - 2) / dfvara_) * vary;
    }
    if(ve.isNotNull()){
        vare_ = as<double>(ve);
    }else{
        vare_ = vary;
    }
    if(dfve.isNotNull()){
        dfvare_ = as<double>(dfve);
    }else{
        dfvare_ = -2;
    }
    if(s2vg.isNotNull()){
        s2vara_ = as<double>(s2vg);
    }else{
        s2vara_ = (vara_ / ((1 - pi_[0]) * sum(vx))) * (dfvara_ - 2) / dfvara_;
    }
    if(s2ve.isNotNull()){
        s2vare_ = as<double>(s2ve);
    }else{
        s2vare_ = 0;
    }
    if(niter < nburn){
        throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
    }

    int n_fold = 2;
    NumericMatrix pi_store(niter - nburn + 1, n_fold);
    NumericVector fold_snp_num(n_fold);
    NumericVector logpi(n_fold), s(n_fold);

    // for gwas
    R_xlen_t nw;
    int ni;
    double varw;
    NumericVector windindx_(m);
    Rcpp::List windx(m);
    NumericVector wu(n);
    double* dwu = NUMERIC_POINTER(wu);
    IntegerVector windxi;
    NumericVector wppai(m);
    NumericVector wgvei(m);
    if(windindx.isNotNull()){
        windindx_ = as<IntegerVector>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = wppai[seq(0, (nw - 1))];
        wgvei = wgvei[seq(0, (nw - 1))];
        for(R_xlen_t w = 0; w < nw; w++){
            windx[w] = which_c(windindx_, (w+1), 5);
        }
    }

    double varg = (vara_ / ((1 - pi_[0]) * sum(vx)));

    if(verbose){
        Rcpp::Rcout.precision(4);
        Rcpp::Rcout << "Prior parameters:" << std::endl;
        if(fixpi){
            Rcpp::Rcout << "    Model fitted at [BayesB]" << std::endl;
        }else{
            Rcpp::Rcout << "    Model fitted at [BayesBπ]" << std::endl;
        }
        Rcpp::Rcout << "    Number of observations " << n << std::endl;
        Rcpp::Rcout << "    Number of covariates " << nc << std::endl;
        Rcpp::Rcout << "    Number of envir-random effects " << nr << std::endl;
        Rcpp::Rcout << "    Number of markers " << m << std::endl;
        Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
        Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
        Rcpp::Rcout << "    Pi at " << pi_ << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << s2vara_ << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << s2vare_ << std::endl;
        Rcpp::Rcout << "    Phenotypic var " << std::fixed << vary << std::endl;
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
    double mu_, mu = mean(y);
    NumericVector one(n); one.fill(1);
    double* done = NUMERIC_POINTER(one);
    NumericVector yadj = y - mu;
    double* dyadj = NUMERIC_POINTER(yadj);
    double* dxi;
    double* pX = NUMERIC_POINTER(X);

    // MCMC procedure
    for(int iter = 0; iter < niter; iter++){

        // sample intercept
        mu_ = - norm_sample(mean(yadj), sqrt(vare_ / n));
        mu -= mu_;
        daxpy_(&n, &mu_, done, &inc, dyadj, &inc);
        
        for(int i = 0; i < nc; i++){
            dci = pC+(long long)i*n;
			oldgi = beta[i];
            v = cpc[i];
			rhs = ddot_(&n, dci, &inc, dyadj, &inc);
			rhs += v * oldgi;
			gi = norm_sample(rhs / v, sqrt(vare_ / v));
			gi_ = oldgi - gi;
			daxpy_(&n, &gi_, dci, &inc, dyadj, &inc);
            beta[i] = gi;
        }

        // log transform
        logpi = log(pi_);
        s[0] = logpi[0];
        sumvg = 0;

        // loop on snps
        for(int i = 0; i < m; i++){

            dxi = pX+(long long)i*n;
            // xi = X[i];   // for list
            // dxi = NUMERIC_POINTER(xi);
            
            xx = xpx[i];
            oldgi = g[i];

            // sample variance
            varg = (oldgi * oldgi + s2vara_ * dfvara_) / (dfvara_ + 1);
            varg = invchisq_sample(1 + dfvara_, varg);

            // right hand
            rhs = ddot_(&n, dxi, &inc, dyadj, &inc);
            if(oldgi){rhs += xx * oldgi;}

            // left hand
            lhs = xx / vare_;

            logdetV = log(varg * lhs + 1);
            uhat = rhs / (xx + vare_ / varg);
            s[1] = -0.5 * (logdetV - (rhs * uhat / vare_)) + logpi[1];
            acceptProb = 1 / sum(exp(s - s[0]));

            // group assignation
            r = uniform_sample();
            indistflag = r < acceptProb ? 0 : 1;
            snptracker[i] = indistflag;

            if(indistflag){

                v = xx + vare_ / varg;

                // gibbs sample snp effect
                gi = norm_sample(rhs / v, sqrt(vare_ / v));

                // update residual
                gi_ = oldgi - gi;
                daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);

                // update random effect
                gi_ *= -1;
                daxpy_(&n, &gi_, dxi, &inc, du, &inc);

                sumvg += (gi * gi) * vx[i];

            }else{

                // zero effect
                gi = 0;
                if(oldgi){
                    gi_ = oldgi;

                    // update residual
                    daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);

                    // update random effect
                    gi_ *= -1;
                    daxpy_(&n, &gi_, dxi, &inc, du, &inc);
                }
            }

            // update snp effect
            g[i] = gi;

        }

        fold_snp_num[1] = sum(snptracker);
        fold_snp_num[0] = m - fold_snp_num[1];
        NnzSnp = fold_snp_num[1];

        // genetic variance
        vara_ = var(u);
    
        // sample residual variance from inv-chisq distribution
        vare_ = (ddot_(&n, dyadj, &inc, dyadj, &inc) + s2vare_ * dfvare_) / (n + dfvare_);
        vare_ = invchisq_sample(n + dfvare_, vare_);

        if(!fixpi){
            pi_ = rdirichlet_sample(n_fold, fold_snp_num + 1);
        }

        if(iter >= nburn){
            count++;
            mu_store[iter - nburn] = mu;
            pi_store(iter - nburn, _) = pi_;
            sumvg_store[iter - nburn] = sumvg;
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            hsq_store[iter - nburn] = vara_ / (vara_ + vare_);
            // g_store += g;
            // nzrate += snptracker;
            daxpy_(&m, &doc, g.begin(), &inc, g_store.begin(), &inc);
            daxpy_(&m, &doc, snptracker.begin(), &inc, nzrate.begin(), &inc);
            if(nc) daxpy_(&nc, &doc, beta.begin(), &inc, beta_store.begin(), &inc);

            // record pve for each window
            if(WPPA){
                for(R_xlen_t w = 0; w < nw; w++){
                    wu.fill(0);
                    windxi = windx[w];
                    for(int i = 0; i < windxi.size(); i++){
                        ni = windxi[i];
                        gi_ = g[ni];
                        if(gi_){
                            dxi = pX+(long long)ni*n;
                            daxpy_(&n, &gi_, dxi, &inc, dwu, &inc);
                        }
                        // if(g[windxi[i]]){wu += X(_, windxi[i]) * g[windxi[i]];}
                    }
                    varw = var(wu);
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

    g = g_store / count;
    beta = beta_store / count;
    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    nzrate = nzrate / count;
    double Mu = mean(mu_store);
    double musd = sd(mu_store);
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
        Rcpp::Rcout << "    Mu " << std::fixed << Mu << " ± " << std::fixed << musd << std::endl;
        for(int i = 0; i < n_fold; i++){
            Rcpp::Rcout << "    Pi" << i + 1 << " " << std::fixed << pi_[i] << " ± " << std::fixed << sd(pi_store(_, i)) << std::endl;;
        }
        Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << " ± " << std::fixed << sumvgsd << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << " ± " << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << " ± " << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated hsq " << std::fixed << hsq << " ± " << std::fixed << hsqsd << std::endl;
        if(nc){
            Rcpp::Rcout << "    Estimated effects for covariates: ";
            for(int i = 0; i < nc; i++){
                Rcpp::Rcout << std::fixed << beta[i] << " ";
            }
            Rcpp::Rcout << std::endl;
        }
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
            Named("mu") = Mu, 
            Named("beta") = beta,
            Named("pi") = pi_, 
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("g") = g,
            Named("nzrate") = nzrate,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("mu") = Mu, 
            Named("beta") = beta,
            Named("pi") = pi_, 
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("g") = g,
            Named("nzrate") = nzrate
        );
    }
}

// [[Rcpp::export]]
Rcpp::List BayesB(
    const NumericVector &y,
    const NumericMatrix &X,
    const Nullable<NumericMatrix> C = R_NilValue,
    const double pi = 0.95,
    const int niter = 50000,
    const int nburn = 20000,
    const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
    const int outfreq = 100,
    const bool verbose = true
){
	Rcpp::List res = BayesBpi(y, X, C, pi, niter, nburn, windindx, wppa, vg, dfvg, s2vg, ve, dfve, s2ve, outfreq, true, verbose);
	return res;
}

// [[Rcpp::export]]
Rcpp::List BayesCpi(
    const NumericVector &y,
    const NumericMatrix &X,
    const Nullable<NumericMatrix> C = R_NilValue,
    const double pi = 0.95,
    const int niter = 50000,
    const int nburn = 20000,
    const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
    const int outfreq = 100,
    const bool fixpi = false,
    const bool verbose = true
){
    if(yhasNA(y)){
        throw Rcpp::exception("NAs are not allowed in y.");
    }
    // if(xhasNA(X)){
    //     throw Rcpp::exception("NAs are not allowed in X.");
    // }
    if(y.length() != X.nrow()){
        throw Rcpp::exception("Number of individuals not equals.");
    }
    double* dci;
    double* pC;
    int nc = 0;
    NumericMatrix C_;
    if(C.isNotNull()){
        C_ = as<NumericMatrix>(C);
        if(C_.nrow() != X.nrow()){
            throw Rcpp::exception("Number of individuals not match for covariates.");
        }
        if(xhasNA(C_)){
            throw Rcpp::exception("NAs are not allowed in covariates.");
        }
        pC = NUMERIC_POINTER(C_);
        nc = C_.ncol();
    }
    NumericVector beta(nc);
    NumericVector beta_store(nc);
    NumericVector cpc(nc);
    if(nc){
        for(int i = 0; i < nc; i++){
            NumericVector ci = C_(_, i);
            cpc[i] = sum(ci * ci);
        }
    }

    int nr = 0;

    int n = y.length();
    int m = X.ncol();
    bool WPPA = false;
    int NnzSnp, indistflag, count = 0;
    int inc = 1;
    double doc = 1.0;
    double xx, oldgi, gi, gi_, rhs, lhs, logdetV, acceptProb, uhat, r, v;
    double vara_, dfvara_, s2vara_, vare_, dfvare_, s2vare_, sumvg, vargi, hsq;
    NumericVector pi_;
    NumericVector snptracker(m);
    NumericVector nzrate(m);
    NumericVector g(m);
    NumericVector g_store(m);
    NumericVector xi(n);
    NumericVector u(n);
    double* du = NUMERIC_POINTER(u);
    NumericVector xpx(m), vx(m);
    NumericVector mu_store(niter - nburn + 1), sumvg_store(niter - nburn + 1), vara_store(niter - nburn + 1), vare_store(niter - nburn + 1), hsq_store(niter - nburn + 1);
    for(int i = 0; i < m; i++){
        xi = X(_, i);
        xpx[i] = sum(xi * xi);
        vx[i] = var(xi);
    }
    double vary = var(y);
    if(pi <= 0 || pi >= 1){
        throw Rcpp::exception("pi should be at 0 < pi < 1.");
    }else{
        pi_ = {pi, 1 - pi};
    }
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
        vara_ = ((dfvara_ - 2) / dfvara_) * vary;
    }
    if(ve.isNotNull()){
        vare_ = as<double>(ve);
    }else{
        vare_ = vary;
    }
    if(dfve.isNotNull()){
        dfvare_ = as<double>(dfve);
    }else{
        dfvare_ = -2;
    }
    if(s2vg.isNotNull()){
        s2vara_ = as<double>(s2vg);
    }else{
        s2vara_ = (vara_ / ((1 - pi_[0]) * sum(vx))) * (dfvara_ - 2) / dfvara_;
    }
    if(s2ve.isNotNull()){
        s2vare_ = as<double>(s2ve);
    }else{
        s2vare_ = 0;
    }
    if(niter < nburn){
        throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
    }

    int n_fold = 2;
    NumericMatrix pi_store(niter - nburn + 1, n_fold);
    NumericVector fold_snp_num(n_fold);
    NumericVector logpi(n_fold), s(n_fold);

    // for gwas
    R_xlen_t nw;
    int ni;
    double varw;
    NumericVector windindx_(m);
    Rcpp::List windx(m);
    NumericVector wu(n);
    double* dwu = NUMERIC_POINTER(wu);
    IntegerVector windxi;
    NumericVector wppai(m);
    NumericVector wgvei(m);
    if(windindx.isNotNull()){
        windindx_ = as<IntegerVector>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = wppai[seq(0, (nw - 1))];
        wgvei = wgvei[seq(0, (nw - 1))];
        for(R_xlen_t w = 0; w < nw; w++){
            windx[w] = which_c(windindx_, (w+1), 5);
        }
    }

    double varg = (vara_ / ((1 - pi_[0]) * sum(vx)));

    if(verbose){
        Rcpp::Rcout.precision(4);
        Rcpp::Rcout << "Prior parameters:" << std::endl;
        if(fixpi){
            Rcpp::Rcout << "    Model fitted at [BayesC]" << std::endl;
        }else{
            Rcpp::Rcout << "    Model fitted at [BayesCπ]" << std::endl;
        }
        Rcpp::Rcout << "    Number of observations " << n << std::endl;
        Rcpp::Rcout << "    Number of covariates " << nc << std::endl;
        Rcpp::Rcout << "    Number of envir-random effects " << nr << std::endl;
        Rcpp::Rcout << "    Number of markers " << m << std::endl;
        Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
        Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
        Rcpp::Rcout << "    Pi at " << pi_ << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << s2vara_ << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << s2vare_ << std::endl;
        Rcpp::Rcout << "    Phenotypic var " << std::fixed << vary << std::endl;
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
    double mu_, mu = mean(y);
    NumericVector one(n); one.fill(1);
    double* done = NUMERIC_POINTER(one);
    NumericVector yadj = y - mu;
    double* dyadj = NUMERIC_POINTER(yadj);
    double* dxi;
    double* pX = NUMERIC_POINTER(X);

    // MCMC procedure
    for(int iter = 0; iter < niter; iter++){

        // sample intercept
        mu_ = - norm_sample(mean(yadj), sqrt(vare_ / n));
        mu -= mu_;
        daxpy_(&n, &mu_, done, &inc, dyadj, &inc);
        
        for(int i = 0; i < nc; i++){
            dci = pC+(long long)i*n;
			oldgi = beta[i];
            v = cpc[i];
			rhs = ddot_(&n, dci, &inc, dyadj, &inc);
			rhs += v * oldgi;
			gi = norm_sample(rhs / v, sqrt(vare_ / v));
			gi_ = oldgi - gi;
			daxpy_(&n, &gi_, dci, &inc, dyadj, &inc);
            beta[i] = gi;
        }

        // log transform
        logpi = log(pi_);
        s[0] = logpi[0];
        sumvg = 0;
        vargi = 0;

        // loop on snps
        for(int i = 0; i < m; i++){

            dxi = pX+(long long)i*n;
            // xi = X[i];   // for list
            // dxi = NUMERIC_POINTER(xi);
            
            xx = xpx[i];
            oldgi = g[i];

            // right hand
            rhs = ddot_(&n, dxi, &inc, dyadj, &inc);
            if(oldgi){rhs += xx * oldgi;}

            // left hand
            lhs = xx / vare_;

            logdetV = log(varg * lhs + 1);
            uhat = rhs / (xx + vare_ / varg);
            s[1] = -0.5 * (logdetV - (rhs * uhat / vare_)) + logpi[1];
            acceptProb = 1 / sum(exp(s - s[0]));

            // group assignation
            r = uniform_sample();
            indistflag = r < acceptProb ? 0 : 1;
            snptracker[i] = indistflag;

            if(indistflag){

                v = xx + vare_ / varg;

                // gibbs sample snp effect
                gi = norm_sample(rhs / v, sqrt(vare_ / v));

                // update residual
                gi_ = oldgi - gi;
                daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);

                // update random effect
                gi_ *= -1;
                daxpy_(&n, &gi_, dxi, &inc, du, &inc);

                sumvg += (gi * gi) * vx[i];
                vargi += (gi * gi);

            }else{

                // zero effect
                gi = 0;
                if(oldgi){
                    gi_ = oldgi;

                    // update residual
                    daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);

                    // update random effect
                    gi_ *= -1;
                    daxpy_(&n, &gi_, dxi, &inc, du, &inc);
                }
            }

            // update snp effect
            g[i] = gi;

        }

        fold_snp_num[1] = sum(snptracker);
        fold_snp_num[0] = m - fold_snp_num[1];
        NnzSnp = fold_snp_num[1];

        varg = (vargi + s2vara_ * dfvara_) / (dfvara_ + NnzSnp);
        varg = invchisq_sample(NnzSnp + dfvara_, varg);

        // genetic variance
        vara_ = var(u);
    
        // sample residual variance from inv-chisq distribution
        vare_ = (ddot_(&n, dyadj, &inc, dyadj, &inc) + s2vare_ * dfvare_) / (n + dfvare_);
        vare_ = invchisq_sample(n + dfvare_, vare_);

        if(!fixpi){
            pi_ = rdirichlet_sample(n_fold, fold_snp_num + 1);
        }

        if(iter >= nburn){
            count++;
            mu_store[iter - nburn] = mu;
            pi_store(iter - nburn, _) = pi_;
            sumvg_store[iter - nburn] = sumvg;
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            hsq_store[iter - nburn] = vara_ / (vara_ + vare_);
            // g_store += g;
            // nzrate += snptracker;
            daxpy_(&m, &doc, g.begin(), &inc, g_store.begin(), &inc);
            daxpy_(&m, &doc, snptracker.begin(), &inc, nzrate.begin(), &inc);
            if(nc) daxpy_(&nc, &doc, beta.begin(), &inc, beta_store.begin(), &inc);


            // record pve for each window
            if(WPPA){
                for(R_xlen_t w = 0; w < nw; w++){
                    wu.fill(0);
                    windxi = windx[w];
                    for(int i = 0; i < windxi.size(); i++){
                        ni = windxi[i];
                        gi_ = g[ni];
                        if(gi_){
                            dxi = pX+(long long)ni*n;
                            daxpy_(&n, &gi_, dxi, &inc, dwu, &inc);
                        }
                        // if(g[windxi[i]]){wu += X(_, windxi[i]) * g[windxi[i]];}
                    }
                    varw = var(wu);
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

    g = g_store / count;
    beta = beta_store / count;
    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    nzrate = nzrate / count;
    double Mu = mean(mu_store);
    double musd = sd(mu_store);
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
        Rcpp::Rcout << "    Mu " << std::fixed << Mu << " ± " << std::fixed << musd << std::endl;
        for(int i = 0; i < n_fold; i++){
            Rcpp::Rcout << "    Pi" << i + 1 << " " << std::fixed << pi_[i] << " ± " << std::fixed << sd(pi_store(_, i)) << std::endl;;
        }
        Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << " ± " << std::fixed << sumvgsd << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << " ± " << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << " ± " << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated hsq " << std::fixed << hsq << " ± " << std::fixed << hsqsd << std::endl;
        if(nc){
            Rcpp::Rcout << "    Estimated effects for covariates: ";
            for(int i = 0; i < nc; i++){
                Rcpp::Rcout << std::fixed << beta[i] << " ";
            }
            Rcpp::Rcout << std::endl;
        }
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
            Named("mu") = Mu, 
            Named("beta") = beta,
            Named("pi") = pi_, 
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("g") = g,
            Named("nzrate") = nzrate,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("mu") = Mu, 
            Named("beta") = beta,
            Named("pi") = pi_, 
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("g") = g,
            Named("nzrate") = nzrate
        );
    }
}

// [[Rcpp::export]]
Rcpp::List BayesC(
    const NumericVector &y,
    const NumericMatrix &X,
    const Nullable<NumericMatrix> C = R_NilValue,
    const double pi = 0.95,
    const int niter = 50000,
    const int nburn = 20000,
    const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
    const int outfreq = 100,
    const bool verbose = true
){
    Rcpp::List res = BayesCpi(y, X, C, pi, niter, nburn, windindx, wppa, vg, dfvg, s2vg, ve, dfve, s2ve, outfreq, true, verbose);
    return res;
}

// [[Rcpp::export]]
Rcpp::List BayesLASSO(
    const NumericVector &y,
    const NumericMatrix &X,
    const Nullable<NumericMatrix> C = R_NilValue,
    const int niter = 50000,
    const int nburn = 20000,
    const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
    const int outfreq = 100,
    const bool verbose = true
){
    if(yhasNA(y)){
        throw Rcpp::exception("NAs are not allowed in y.");
    }
    // if(xhasNA(X)){
    //     throw Rcpp::exception("NAs are not allowed in X.");
    // }
    if(y.length() != X.nrow()){
        throw Rcpp::exception("Number of individuals not equals.");
    }
    double* dci;
    double* pC;
    int nc = 0;
    NumericMatrix C_;
    if(C.isNotNull()){
        C_ = as<NumericMatrix>(C);
        if(C_.nrow() != X.nrow()){
            throw Rcpp::exception("Number of individuals not match for covariates.");
        }
        if(xhasNA(C_)){
            throw Rcpp::exception("NAs are not allowed in covariates.");
        }
        pC = NUMERIC_POINTER(C_);
        nc = C_.ncol();
    }
    NumericVector beta(nc);
    NumericVector beta_store(nc);
    NumericVector cpc(nc);
    if(nc){
        for(int i = 0; i < nc; i++){
            NumericVector ci = C_(_, i);
            cpc[i] = sum(ci * ci);
        }
    }

    int nr = 0;

    int n = y.length();
    int m = X.ncol();
    bool WPPA = false;
    int count = 0;
    int inc = 1;
    double doc = 1.0;
    double xx, oldgi, gi, gi_, rhs, v;
    double vara_, dfvara_, s2vara_, vare_, dfvare_, s2vare_, sumvg, vargi, hsq;
    NumericVector g(m);
    NumericVector g_store(m);
    NumericVector xi(n);
    NumericVector u(n);
    double* du = NUMERIC_POINTER(u);
    NumericVector xpx(m), vx(m);
    NumericVector mu_store(niter - nburn + 1), sumvg_store(niter - nburn + 1), vara_store(niter - nburn + 1), vare_store(niter - nburn + 1), hsq_store(niter - nburn + 1);
    for(int i = 0; i < m; i++){
        xi = X(_, i);
        xpx[i] = sum(xi * xi);
        vx[i] = var(xi);
    }
    double vary = var(y);
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
        vara_ = ((dfvara_ - 2) / dfvara_) * vary;
    }
    if(ve.isNotNull()){
        vare_ = as<double>(ve);
    }else{
        vare_ = vary;
    }
    if(dfve.isNotNull()){
        dfvare_ = as<double>(dfve);
    }else{
        dfvare_ = -2;
    }
    if(s2vg.isNotNull()){
        s2vara_ = as<double>(s2vg);
    }else{
        s2vara_ = (vara_ / (sum(vx))) * (dfvara_ - 2) / dfvara_;
    }
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

    // for gwas
    R_xlen_t nw;
    int ni;
    double varw;
    NumericVector windindx_(m);
    Rcpp::List windx(m);
    NumericVector wu(n);
    double* dwu = NUMERIC_POINTER(wu);
    IntegerVector windxi;
    NumericVector wppai(m);
    NumericVector wgvei(m);
    if(windindx.isNotNull()){
        windindx_ = as<IntegerVector>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = wppai[seq(0, (nw - 1))];
        wgvei = wgvei[seq(0, (nw - 1))];
        for(R_xlen_t w = 0; w < nw; w++){
            windx[w] = which_c(windindx_, (w+1), 5);
        }
    }

    NumericVector varg(m); varg.fill(vara_ / sum(vx));

    if(verbose){
        Rcpp::Rcout.precision(4);
        Rcpp::Rcout << "Prior parameters:" << std::endl;
        Rcpp::Rcout << "    Model fitted at [BayesLASSO]" << std::endl;
        Rcpp::Rcout << "    Number of observations " << n << std::endl;
        Rcpp::Rcout << "    Number of covariates " << nc << std::endl;
        Rcpp::Rcout << "    Number of envir-random effects " << nr << std::endl;
        Rcpp::Rcout << "    Number of markers " << m << std::endl;
        Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
        Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << std::fixed << s2vara_ << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << std::fixed << s2vare_ << std::endl;
        Rcpp::Rcout << "    Phenotypic var " << std::fixed << vary << std::endl;
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
    double mu_, mu = mean(y);
    NumericVector one(n); one.fill(1);
    double* done = NUMERIC_POINTER(one);
    NumericVector yadj = y - mu;
    double* dyadj = NUMERIC_POINTER(yadj);
    double* dxi;
    double* pX = NUMERIC_POINTER(X);

    // MCMC procedure
    for(int iter = 0; iter < niter; iter++){

        // sample intercept
        mu_ = - norm_sample(mean(yadj), sqrt(vare_ / n));
        mu -= mu_;
        daxpy_(&n, &mu_, done, &inc, dyadj, &inc);
        sumvg = 0;

        for(int i = 0; i < nc; i++){
            dci = pC+(long long)i*n;
			oldgi = beta[i];
            v = cpc[i];
			rhs = ddot_(&n, dci, &inc, dyadj, &inc);
			rhs += v * oldgi;
			gi = norm_sample(rhs / v, sqrt(vare_ / v));
			gi_ = oldgi - gi;
			daxpy_(&n, &gi_, dci, &inc, dyadj, &inc);
            beta[i] = gi;
        }

        // loop on snps
        for(int i = 0; i < m; i++){

            dxi = pX+(long long)i*n;
            // xi = X[i];   // for list
            // dxi = NUMERIC_POINTER(xi);
            
            xx = xpx[i];
            oldgi = g[i];

            // right hand
            rhs = ddot_(&n, dxi, &inc, dyadj, &inc);
            rhs += xx * oldgi;

            // left hand
            v = xx + 1 / varg[i];

            // sample snp effect
            gi = norm_sample(rhs / v, sqrt(vare_ / v));
            if(abs(gi) < 1e-6){gi = 1e-6;}

            vargi = 1 / rinvgaussian_sample(sqrt(vare_) * lambda / abs(gi), lambda2);
            if(vargi < 0){
                // nothing
            }else{
                varg[i] = vargi;
            }

            // update residual
            gi_ = oldgi - gi;
            daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);

            // update random effect
            gi_ *= -1;
            daxpy_(&n, &gi_, dxi, &inc, du, &inc);

            // update snp effect
            g[i] = gi;
            sumvg += (gi * gi) * vx[i];

        }

        // genetic variance
        vara_ = var(u);

        shape = shape0 + m;
        rate = rate0 + sum(varg) / 2;
        lambda2 = gamma_sample(shape, 1 / rate);
        lambda = sqrt(lambda2);

        // sample residual variance from inv-chisq distribution
        vare_ = (ddot_(&n, dyadj, &inc, dyadj, &inc) + s2vare_ * dfvare_) / (n + dfvare_);
        vare_ = invchisq_sample(n + dfvare_, vare_);

        if(iter >= nburn){
            count++;
            mu_store[iter - nburn] = mu;
            sumvg_store[iter - nburn] = sumvg;
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            hsq_store[iter - nburn] = vara_ / (vara_ + vare_);
            // g_store += g;
            daxpy_(&m, &doc, g.begin(), &inc, g_store.begin(), &inc);
            if(nc) daxpy_(&nc, &doc, beta.begin(), &inc, beta_store.begin(), &inc);

            // record pve for each window
            if(WPPA){
                for(R_xlen_t w = 0; w < nw; w++){
                    wu.fill(0);
                    windxi = windx[w];
                    for(int i = 0; i < windxi.size(); i++){
                        ni = windxi[i];
                        gi_ = g[ni];
                        if(gi_){
                            dxi = pX+(long long)ni*n;
                            daxpy_(&n, &gi_, dxi, &inc, dwu, &inc);
                        }
                        // if(g[windxi[i]]){wu += X(_, windxi[i]) * g[windxi[i]];}
                    }
                    varw = var(wu);
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

    g = g_store / count;
    beta = beta_store / count;
    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    double Mu = mean(mu_store);
    double musd = sd(mu_store);
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
        Rcpp::Rcout << "    Mu " << std::fixed << Mu << " ± " << std::fixed << musd << std::endl;
        Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << " ± " << std::fixed << sumvgsd << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << " ± " << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << " ± " << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated hsq " << std::fixed << hsq << " ± " << std::fixed << hsqsd << std::endl;
        if(nc){
            Rcpp::Rcout << "    Estimated effects for covariates: ";
            for(int i = 0; i < nc; i++){
                Rcpp::Rcout << std::fixed << beta[i] << " ";
            }
            Rcpp::Rcout << std::endl;
        }
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
            Named("mu") = Mu, 
            Named("beta") = beta,
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("g") = g,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("mu") = Mu, 
            Named("beta") = beta,
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("g") = g
        );
    }
}

// [[Rcpp::export]]
Rcpp::List BayesR(
    const NumericVector &y,
    const NumericMatrix &X,
    const Nullable<NumericMatrix> C = R_NilValue,
    const Nullable<NumericVector> pi = R_NilValue,
    const Nullable<NumericVector> fold = R_NilValue,
    const int niter = 50000,
    const int nburn = 20000,
    const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
    const int outfreq = 100,
    const bool fixpi = false,
    const bool verbose = true
){
    if(yhasNA(y)){
        throw Rcpp::exception("NAs are not allowed in y.");
    }
    // if(xhasNA(X)){
    //     throw Rcpp::exception("NAs are not allowed in X.");
    // }
    if(y.length() != X.nrow()){
        throw Rcpp::exception("Number of individuals not equals.");
    }
    double* dci;
    double* pC;
    int nc = 0;
    NumericMatrix C_;
    if(C.isNotNull()){
        C_ = as<NumericMatrix>(C);
        if(C_.nrow() != X.nrow()){
            throw Rcpp::exception("Number of individuals not match for covariates.");
        }
        if(xhasNA(C_)){
            throw Rcpp::exception("NAs are not allowed in covariates.");
        }
        pC = NUMERIC_POINTER(C_);
        nc = C_.ncol();
    }
    NumericVector beta(nc);
    NumericVector beta_store(nc);
    NumericVector cpc(nc);
    if(nc){
        for(int i = 0; i < nc; i++){
            NumericVector ci = C_(_, i);
            cpc[i] = sum(ci * ci);
        }
    }

    int nr = 0;

    int n = y.length();
    int m = X.ncol();
    bool WPPA = false;
    int NnzSnp, indistflag, count = 0;
    int inc = 1;
    double doc = 1.0;
    double xx, oldgi, gi, gi_, rhs, lhs, temp, logdetV, acceptProb, uhat, r, v;
    double vara_, dfvara_, s2vara_, vare_, dfvare_, s2vare_, sumvg, varg, hsq;
    NumericVector pi_, fold_;
    NumericVector snptracker(m);
    NumericVector nzrate(m);
    NumericVector g(m);
    NumericVector g_store(m);
    NumericVector xi(n);
    NumericVector u(n);
    double* du = NUMERIC_POINTER(u);
    NumericVector xpx(m), vx(m);
    NumericVector mu_store(niter - nburn + 1), sumvg_store(niter - nburn + 1), vara_store(niter - nburn + 1), vare_store(niter - nburn + 1), hsq_store(niter - nburn + 1);
    for(int i = 0; i < m; i++){
        xi = X(_, i);
        xpx[i] = sum(xi * xi);
        vx[i] = var(xi);
    }
    double vary = var(y);
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
        vara_ = ((dfvara_ - 2) / dfvara_) * vary;
    }
    if(ve.isNotNull()){
        vare_ = as<double>(ve);
    }else{
        vare_ = vary;
    }
    if(dfve.isNotNull()){
        dfvare_ = as<double>(dfve);
    }else{
        dfvare_ = -2;
    }
    if(s2vg.isNotNull()){
        s2vara_ = as<double>(s2vg);
    }else{
        s2vara_ = (vara_ / ((1 - pi_[0]) * sum(vx))) * (dfvara_ - 2) / dfvara_;
    }
    if(s2ve.isNotNull()){
        s2vare_ = as<double>(s2ve);
    }else{
        s2vare_ = 0;
    }
    if(niter < nburn){
        throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
    }

    int n_fold = fold_.length();
    NumericVector stemp(n_fold);
    NumericMatrix pi_store(niter - nburn + 1, n_fold);
    NumericVector fold_snp_num(n_fold);
    NumericVector logpi(n_fold), s(n_fold);
    NumericVector vara_fold = (vara_ / ((1 - pi_[0]) * sum(vx))) * fold_;
    NumericVector vare_vara_fold(n_fold);

    // for gwas
    R_xlen_t nw;
    int ni;
    double varw;
    NumericVector windindx_(m);
    Rcpp::List windx(m);
    NumericVector wu(n);
    double* dwu = NUMERIC_POINTER(wu);
    IntegerVector windxi;
    NumericVector wppai(m);
    NumericVector wgvei(m);
    if(windindx.isNotNull()){
        windindx_ = as<IntegerVector>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = wppai[seq(0, (nw - 1))];
        wgvei = wgvei[seq(0, (nw - 1))];
        for(R_xlen_t w = 0; w < nw; w++){
            windx[w] = which_c(windindx_, (w+1), 5);
        }
    }
    if(verbose){
        Rcpp::Rcout.precision(4);
        Rcpp::Rcout << "Prior parameters:" << std::endl;
        Rcpp::Rcout << "    Model fitted at [BayesR]" << std::endl;
        Rcpp::Rcout << "    Number of observations " << n << std::endl;
        Rcpp::Rcout << "    Number of covariates " << nc << std::endl;
        Rcpp::Rcout << "    Number of envir-random effects " << nr << std::endl;
        Rcpp::Rcout << "    Number of markers " << m << std::endl;
        Rcpp::Rcout << "    Total number of iteration " << niter << std::endl;
        Rcpp::Rcout << "    Total number of burn-in " << nburn << std::endl;
        Rcpp::Rcout << "    Pi at " << pi_ << std::endl;
        Rcpp::Rcout << "    Group fold " << fold_ << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << s2vara_ << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << s2vare_ << std::endl;
        Rcpp::Rcout << "    Phenotypic var " << std::fixed << vary << std::endl;
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
    double mu_, mu = mean(y);
    NumericVector one(n); one.fill(1);
    double* done = NUMERIC_POINTER(one);
    NumericVector yadj = y - mu;
    double* dyadj = NUMERIC_POINTER(yadj);
    double* dxi;
    double* pX = NUMERIC_POINTER(X);

    // MCMC procedure
    for(int iter = 0; iter < niter; iter++){

        // sample intercept
        mu_ = - norm_sample(mean(yadj), sqrt(vare_ / n));
        mu -= mu_;
        daxpy_(&n, &mu_, done, &inc, dyadj, &inc);

        for(int i = 0; i < nc; i++){
            dci = pC+(long long)i*n;
			oldgi = beta[i];
            v = cpc[i];
			rhs = ddot_(&n, dci, &inc, dyadj, &inc);
			rhs += v * oldgi;
			gi = norm_sample(rhs / v, sqrt(vare_ / v));
			gi_ = oldgi - gi;
			daxpy_(&n, &gi_, dci, &inc, dyadj, &inc);
            beta[i] = gi;
        }

        // log transform
        logpi = log(pi_);
        s[0] = logpi[0];
        sumvg = 0;
        varg = 0;
        for(int j = 1; j < n_fold; j++){
            vare_vara_fold[j] = vare_ / vara_fold[j];
        }

        // loop on snps
        for(int i = 0; i < m; i++){

            dxi = pX+(long long)i*n;
            // xi = X[i];   // for list
            // dxi = NUMERIC_POINTER(xi);
            
            xx = xpx[i];
            oldgi = g[i];

            // right hand
            rhs = ddot_(&n, dxi, &inc, dyadj, &inc);
            if(oldgi){rhs += xx * oldgi;}

            // left hand
            lhs = xx / vare_;

			for(int j = 1; j < n_fold; j++){
				logdetV = log(vara_fold[j] * lhs + 1);
				uhat = rhs / (xx + vare_vara_fold[j]);
				s[j] = -0.5 * (logdetV - (rhs * uhat / vare_)) + logpi[j];
			}

            for(int j = 0; j < n_fold; j++){
                temp = 0.0;
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
            if(indistflag){

                v = xx + vare_vara_fold[indistflag];

                // gibbs sample snp effect
                gi = norm_sample(rhs / v, sqrt(vare_ / v));

                // update residual
                gi_ = oldgi - gi;
                daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);

                // update random effect
                gi_ *= -1;
                daxpy_(&n, &gi_, dxi, &inc, du, &inc);

                sumvg += (gi * gi) * vx[i];
                varg += gi * gi / fold_[indistflag];

            }else{

                // zero effect
                gi = 0;
                if(oldgi){
                    gi_ = oldgi;

                    // update residual
                    daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);

                    // update random effect
                    gi_ *= -1;
                    daxpy_(&n, &gi_, dxi, &inc, du, &inc);
                }
            }

            // update snp effect
            g[i] = gi;

        }

        for(int j = 0; j < n_fold; j++){
            fold_snp_num[j] = sum(snptracker == j);
        }
        NnzSnp = m - fold_snp_num[0];
        varg = (varg + s2vara_ * dfvara_) / (dfvara_ + NnzSnp);
        varg = invchisq_sample(NnzSnp + dfvara_, varg);

        // genetic variance
        vara_ = var(u);
    
        // sample residual variance from inv-chisq distribution
        vare_ = (ddot_(&n, dyadj, &inc, dyadj, &inc) + s2vare_ * dfvare_) / (n + dfvare_);
        vare_ = invchisq_sample(n + dfvare_, vare_);

        for(int j = 0; j < n_fold; j++){
            vara_fold[j] = varg * fold_[j]; 
            // vara_fold[j] = vara_ * fold_[j]; 
        }

        if(!fixpi){
            pi_ = rdirichlet_sample(n_fold, fold_snp_num + 1);
        }

        if(iter >= nburn){
            count++;
            mu_store[iter - nburn] = mu;
            pi_store(iter - nburn, _) = pi_;
            sumvg_store[iter - nburn] = sumvg;
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            hsq_store[iter - nburn] = vara_ / (vara_ + vare_);
            // g_store += g;
            daxpy_(&m, &doc, g.begin(), &inc, g_store.begin(), &inc);
            if(nc) daxpy_(&nc, &doc, beta.begin(), &inc, beta_store.begin(), &inc);

            for(int i = 0; i < m; i++){
                if(snptracker[i]){
                    nzrate[i] += 1;
                }
            }
            
            // record pve for each window
            if(WPPA){
                for(R_xlen_t w = 0; w < nw; w++){
                    wu.fill(0);
                    windxi = windx[w];
                    for(int i = 0; i < windxi.size(); i++){
                        ni = windxi[i];
                        gi_ = g[ni];
                        if(gi_){
                            dxi = pX+(long long)ni*n;
                            daxpy_(&n, &gi_, dxi, &inc, dwu, &inc);
                        }
                        // if(g[windxi[i]]){wu += X(_, windxi[i]) * g[windxi[i]];}
                    }
                    varw = var(wu);
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

    g = g_store / count;
    beta = beta_store / count;
    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    nzrate = nzrate / count;
    double Mu = mean(mu_store);
    double musd = sd(mu_store);
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
        Rcpp::Rcout << "    Mu " << std::fixed << Mu << " ± " << std::fixed << musd << std::endl;
        for(int i = 0; i < n_fold; i++){
            Rcpp::Rcout << "    Pi" << i + 1 << " " << std::fixed << pi_[i] << " ± " << std::fixed << sd(pi_store(_, i)) << std::endl;;
        }
        Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << " ± " << std::fixed << sumvgsd << std::endl;
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << " ± " << std::fixed << varasd << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << " ± " << std::fixed << varesd << std::endl;
        Rcpp::Rcout << "    Estimated hsq " << std::fixed << hsq << " ± " << std::fixed << hsqsd << std::endl;
        if(nc){
            Rcpp::Rcout << "    Estimated effects for covariates: ";
            for(int i = 0; i < nc; i++){
                Rcpp::Rcout << std::fixed << beta[i] << " ";
            }
            Rcpp::Rcout << std::endl;
        }
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
            Named("mu") = Mu,
            Named("beta") = beta, 
            Named("pi") = pi_, 
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("g") = g,
            Named("nzrate") = nzrate,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("mu") = Mu, 
            Named("beta") = beta,
            Named("pi") = pi_, 
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("g") = g,
            Named("nzrate") = nzrate
        );
    }
}
