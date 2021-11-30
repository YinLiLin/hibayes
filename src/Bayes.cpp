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

// [[Rcpp::export]]
Rcpp::List Bayes(
    arma::vec &y,
    arma::mat &X,
    std::string &model,
    arma::vec &pi,
    const Nullable<arma::mat> C = R_NilValue,
    const Nullable<CharacterMatrix> R = R_NilValue,
    const Nullable<arma::vec> fold = R_NilValue,
    const int niter = 50000,
    const int nburn = 20000,
    const Nullable<arma::vec> epsl_y_J = R_NilValue,
    const Nullable<arma::sp_mat> epsl_Gi = R_NilValue,
    const Nullable<arma::uvec> epsl_index = R_NilValue,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
    const Nullable<IntegerVector> windindx = R_NilValue,
    const double wppa = 0.01,
    const int outfreq = 100,
    const int threads = 0,
    const bool verbose = true
){

    omp_setup(threads);

    if(y.has_nan())   throw Rcpp::exception("NAs are not allowed in y.");
    // if(X.has_nan()){
    //     throw Rcpp::exception("NAs are not allowed in genotype.");
    // }
    if(y.n_elem != X.n_rows)  throw Rcpp::exception("Number of individuals not equals.");
    int model_index = (model == "BayesRR" ? 1 : (model == "BayesA" ? 2 : (model == "BayesB" || model == "BayesBpi" ? 3 : (model == "BayesC" || model == "BayesCpi" ? 4 : (model == "BayesL" ? 5 : 6)))));
    bool fixpi = false;
    if(model == "BayesB" || model == "BayesC")    fixpi = true;
    if(pi.n_elem < 2)  throw Rcpp::exception("pi should be a vector.");
    if(sum(pi) != 1)   throw Rcpp::exception("sum of pi should be 1.");
    if(pi[0] == 1) throw Rcpp::exception("all markers have no effect size.");
    for(int i = 0; i < pi.n_elem; i++){
        if(pi[i] < 0 || pi[i] > 1){
            throw Rcpp::exception("elements of pi should be at the range of [0, 1]");
        }
    }
    vec fold_;
    if(fold.isNotNull()){
        fold_ = as<arma::vec>(fold);
    }else{
        if(model == "BayesR")    throw Rcpp::exception("'fold' should be provided for BayesR model.");
        fold_.resize(2);
    }
    if(fold_.n_elem != pi.n_elem){
        throw Rcpp::exception("length of pi and fold not equals.");
    }

    double vary = var(y);
    double h2 = 0.5;
    
    double* dci;
    int nc = 0;
    mat C_;
    vec beta;
    vec betase;
    mat beta_store;
    if(C.isNotNull()){
        C_ = as<arma::mat>(C);
        if(C_.n_rows != X.n_rows)   throw Rcpp::exception("Number of individuals does not match for covariates.");
        if(C_.has_nan())  throw Rcpp::exception("Individuals with phenotypic value should not have missing covariates.");
        nc = C_.n_cols;
        beta.resize(nc);
        betase.resize(nc);
        beta_store.resize(niter - nburn + 1, nc);
    }
    vec cpc;
    if(nc){
        cpc.resize(nc);
        #pragma omp parallel for
        for(int i = 0; i < nc; i++){
            cpc[i] = dot(C_.col(i), C_.col(i));
        }
    }
    
    int nr = 0;
    vec vr;
    vec vrtmp;
    vec vrsd;
    sp_mat r_LHS;
    vec estR_tmp;
    vec r_RHS;
    vec diff;
    int qr;
    int dfr = 5;
    double s2r;
    mat vr_store;
    vector<sp_mat> Z;
    vector<sp_mat> ZZ;
    vector<vec> estR;
    vector<vec> estR_store;
    vector< vector<string> > Z_levels;
    if(R.isNotNull()){
        CharacterMatrix R_ = as<CharacterMatrix>(R);
        if(R_.nrow() != X.n_rows)   throw Rcpp::exception("Number of individuals does not match for environmental random effects.");
        if(xhasNA(R_))  throw Rcpp::exception("Individuals with phenotypic value should not have missing environmental random effects.");
        nr = R_.ncol();
        vr.resize(nr);
        vrtmp.resize(nr);
        vr_store.resize(niter - nburn + 1, nr);
        vrtmp.fill(vary * (1 - h2) / (nr + 1));
        s2r = vrtmp[0] * (dfr - 2) / dfr;
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

    int ne = 0;
    double veps;
    double vepstmp;
    double vepsse;
    double JtJ;
    vec veps_store;
    sp_mat epsl_Gi_, epsl_Z, epsl_ZZ;
    uvec epsl_index_;
    vec epsl_estR;
    vec epsl_estR_store;
    vec epsl_y_J_;
    sp_mat epsl_LHS;
    vec epsl_yadj;
    vec epsl_estR_tmp;
    vec epsl_RHS;
    int qe = 0;
    double epsl_J_beta = 0;
    vec epsl_J_beta_store;
    if(epsl_index.isNotNull()){
        epsl_index_ = as<uvec>(epsl_index);
        epsl_index_ -= 1;
        ne = epsl_index_.n_elem;
    }
    if(ne){
        if(!epsl_Gi.isNotNull())    throw Rcpp::exception("variance-covariance matrix should be provided for epsilon term.");
        epsl_Gi_ = as<sp_mat>(epsl_Gi);
        epsl_y_J_ = as<vec>(epsl_y_J);
        if(epsl_Gi_.n_cols != epsl_Gi_.n_rows) throw Rcpp::exception("variance-covariance matrix should be in square.");
        JtJ = dot(epsl_y_J_, epsl_y_J_);
        veps_store.resize(niter - nburn + 1);
        epsl_J_beta_store.resize(niter - nburn + 1);
        epsl_Z.resize(ne, epsl_Gi_.n_cols);
        epsl_estR.zeros(epsl_Gi_.n_cols);
        epsl_estR_store.zeros(epsl_Gi_.n_cols);
        epsl_yadj.resize(ne);
        estR_tmp.resize(epsl_Gi_.n_cols);
        qe = epsl_Gi_.n_cols;
        for(int i = 0; i < ne; i++){epsl_Z(i, epsl_index_[i]) = 1.0;}
        epsl_ZZ = epsl_Z.t() * epsl_Z;
    }

    int n = y.n_elem;
    int m = X.n_cols;
    int count = 0;
    int inc = 1;
    int n_fold = fold_.n_elem;
    int NnzSnp, indistflag;
    double doc = 1.0;
    double xx, oldgi, gi, gi_, rhs, lhs, logdetV, acceptProb, uhat, v;
    double vara_, dfvara_, s2vara_, vare_, dfvare_, s2vare_, vargi, hsq, s2varg_;
    vec snptracker;
    vec nzrate;
    if(model == "BayesRR" || model == "BayesA" || model == "BayesL"){
        NnzSnp = m;
        pi[0] = 0; pi[1] = 1;
        fixpi = true;
    }else{
        if(model != "BayesR" && pi.n_elem != 2) throw Rcpp::exception("length of pi should be 2, the first value is the proportion of non-effect markers.");
        nzrate.zeros(m); 
        snptracker.zeros(m); 
    }
    vec g = zeros(m);
    vec g_store = zeros(m);
    vec u = zeros(n);
    vec xpx = zeros(m);
    vec vx = zeros(m);
    vec mu_store = zeros(niter - nburn + 1);
    // vec sumvg_store = zeros(niter - nburn + 1);
    vec vara_store = zeros(niter - nburn + 1);
    vec vare_store = zeros(niter - nburn + 1);
    vec hsq_store = zeros(niter - nburn + 1);
    mat pi_store;
    if(!fixpi)  pi_store.resize(niter - nburn + 1, n_fold);

    #pragma omp parallel for
    for(int i = 0; i < m; i++){
        vec Xi = X.col(i);
        xpx[i] = sum(square(Xi));
        vx[i] = var(Xi);
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
        vara_ = ((dfvara_ - 2) / dfvara_) * vary * h2;
    }
    vepstmp = vara_;
    if(ve.isNotNull()){
        vare_ = as<double>(ve);
    }else{
        vare_ = vary * (1 - h2) / (nr + 1);
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
    double varg = vara_ / ((1 - pi[0]) * sum(vx));
    s2varg_ = s2vara_ / ((1 - pi[0]) * sum(vx));
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
    vec vara_fold = (vara_ / ((1 - pi[0]) * sum(vx))) * fold_;
    vec vare_vara_fold = zeros(n_fold);

    // for gwas
    int nw;
    int ni;
    double varw;
    bool WPPA = false;
    NumericVector windindx_;
    vector<IntegerVector> windx;
    vec wu;
    IntegerVector windxi;
    NumericVector wppai;
    NumericVector wgvei;
    if(windindx.isNotNull()){
        windindx_ = as<IntegerVector>(windindx);
        WPPA = true;
        nw = max(windindx_);
        wppai = seq(0, (nw - 1));
        wgvei = seq(0, (nw - 1));
        for(int w = 0; w < nw; w++){
            windx.push_back(which_c(windindx_, (w+1), 5));
        }
        wu.zeros(n);
    }

    if(verbose){
        Rcpp::Rcout.precision(4);
        Rcpp::Rcout << "Prior parameters:" << std::endl;
        Rcpp::Rcout << "    Model fitted at [" << (model == "BayesRR" ? "Bayes Ridge Regression" : model) << "]" << std::endl;
        Rcpp::Rcout << "    Number of observations " << n << std::endl;
        if(epsl_index.isNotNull()){
            Rcpp::Rcout << "    Observations with genotype " << (n - ne) << std::endl;
            Rcpp::Rcout << "    Observations with imputed genotype " << ne << std::endl;
        }
        Rcpp::Rcout << "    Number of covariates " << nc << std::endl;
        Rcpp::Rcout << "    Number of envir-random effects " << nr << std::endl;
        Rcpp::Rcout << "    Number of markers " << m << std::endl;
        for(int i = 0; i < pi.n_elem; i++){
            if(i == 0){
                Rcpp::Rcout << "    π for markers in zero effect size " << pi[i] << endl;
            }else{
                if(i == 1){
                    Rcpp::Rcout << "    π for markers in non-zero effect size " << pi[i] << " ";
                }else{
                    Rcpp::Rcout << pi[i] << " ";
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
        if(nr){
            Rcpp::Rcout << "    Environmental var ";
            for(int i = 0; i < nr; i++){
                Rcpp::Rcout << std::fixed << vrtmp[i] << " ";
            }
            Rcpp::Rcout << std::endl;
        }
        Rcpp::Rcout << "    Genetic var " << std::fixed << vara_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq gpar " << std::fixed << dfvara_ << " " << std::fixed << s2vara_ << std::endl;
        Rcpp::Rcout << "    Residual var " << std::fixed << vare_ << std::endl;
        Rcpp::Rcout << "    Inv-Chisq epar " << std::fixed << dfvare_ << " " << std::fixed << s2vare_ << std::endl;
        Rcpp::Rcout << "    Marker var " << std::fixed << varg << std::endl;
        Rcpp::Rcout << "    Inv-Chisq alpar " << std::fixed << dfvara_ << " " << std::fixed << s2varg_ << std::endl;
        
        if(WPPA){
            Rcpp::Rcout << "    Number of windows " << nw << std::endl;
            Rcpp::Rcout << "    GVE threshold for windows " << wppa << std::endl;
        }
        Rcpp::Rcout << "MCMC started: " << std::endl;
        Rcpp::Rcout << " Iter" << "  ";
        Rcpp::Rcout << "NumNZSnp" << "  ";
        for(int i = 0; i < n_fold; i++){
            Rcpp::Rcout << "π" << i + 1 << "  ";
        }
        if(model == "BayesL")    Rcpp::Rcout << "Lambda" << "  ";
        for(int i = 0; i < nr; i++){
            Rcpp::Rcout << "Vr" << (i + 1) << "  ";
        }
        if(ne)  Rcpp::Rcout << "Vε" << "  ";
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
    double mu_, mu = mean(y);
    vec one(n); one.fill(1);
    vec yadj = y - mu;
    double* dyadj = yadj.memptr();
    double* dxi;
    // double* pX = NUMERIC_POINTER(X);

    // MCMC procedure
    for(int iter = 0; iter < niter; iter++){

        // sample intercept
        mu_ = - norm_sample(mean(yadj), sqrt(vare_ / n));
        mu -= mu_;
		daxpy_(&n, &mu_, one.memptr(), &inc, dyadj, &inc);

        for(int i = 0; i < nc; i++){
            dci = C_.colptr(i);
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
            r_LHS = ZZ[i];
            r_LHS.diag() += vare_ / vrtmp[i];
            estR_tmp = estR[i];
            r_RHS = Z[i].t() * yadj;
            r_RHS += ZZ[i] * estR_tmp;
            Gibbs(r_LHS, estR_tmp, r_RHS, vare_);
            estR[i] -= estR_tmp;
            diff = Z[i] * estR[i];
            daxpy_(&n, &doc, diff.memptr(), &inc, dyadj, &inc);
            qr = estR_tmp.n_elem;
            vrtmp[i] = (ddot_(&qr, estR_tmp.memptr(), &inc, estR_tmp.memptr(), &inc) + s2r * dfr) / (qr + dfr);
            vrtmp[i] = invchisq_sample(qr + dfr, vrtmp[i]);
            vr[i] = var(estR_tmp);
            estR[i] = estR_tmp;
        }

        if(ne){
            oldgi = epsl_J_beta;
            v = JtJ;
			rhs = ddot_(&n, epsl_y_J_.memptr(), &inc, dyadj, &inc);
			rhs += v * oldgi;
			gi = norm_sample(rhs / v, sqrt(vare_ / v));
			gi_ = oldgi - gi;
			daxpy_(&n, &gi_, epsl_y_J_.memptr(), &inc, dyadj, &inc);
            gi_ *= -1;
            daxpy_(&n, &gi_, epsl_y_J_.memptr(), &inc, u.memptr(), &inc);
            epsl_J_beta = gi;
            epsl_LHS = epsl_ZZ;
            epsl_LHS += (epsl_Gi_ * (vare_ / vepstmp));
            epsl_yadj = yadj.tail(ne);
            epsl_estR_tmp = epsl_estR;
            epsl_RHS = epsl_Z.t() * epsl_yadj;
            epsl_RHS += epsl_ZZ * epsl_estR_tmp;
            Gibbs(epsl_LHS, epsl_estR_tmp, epsl_RHS, vare_);
            epsl_estR -= epsl_estR_tmp;
            epsl_yadj = epsl_Z * epsl_estR;
            yadj.tail(ne) += epsl_yadj;
            u.tail(ne) -= epsl_yadj;
            vepstmp = as_scalar(epsl_estR_tmp.t() * epsl_Gi_ * epsl_estR_tmp);
            vepstmp += s2vara_ * dfvara_;
            vepstmp /= (dfvara_ + qe);
            vepstmp = invchisq_sample(qe + dfvara_, vepstmp);
            epsl_estR = epsl_estR_tmp;
            veps = var(epsl_estR);
        }

        switch(model_index){
            case 1:
                for(int i = 0; i < m; i++){
                    dxi = X.colptr(i);
                    xx = xpx[i];
                    oldgi = g[i];
                    rhs = ddot_(&n, dxi, &inc, dyadj, &inc);
                    rhs += xx * oldgi;
                    v = xx + vare_ / varg;
                    gi = norm_sample(rhs / v, sqrt(vare_ / v));
                    gi_ = oldgi - gi;
                    daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);
                    gi_ *= -1;
                    daxpy_(&n, &gi_, dxi, &inc, u.memptr(), &inc);
                    g[i] = gi;

                }
                varg = (ddot_(&m, g.memptr(), &inc, g.memptr(), &inc) + s2varg_ * dfvara_) / (dfvara_ + m);
                varg = invchisq_sample(m + dfvara_, varg);
                // sumvg = sum(vx * varg);
                break;
            case 2:
                for(int i = 0; i < m; i++){
                    dxi = X.colptr(i);
                    xx = xpx[i];
                    oldgi = g[i];
                    varg = (oldgi * oldgi + s2varg_ * dfvara_) / (dfvara_ + 1);
                    varg = invchisq_sample(1 + dfvara_, varg);
                    rhs = ddot_(&n, dxi, &inc, dyadj, &inc);
                    rhs += xx * oldgi;
                    v = xx + vare_ / varg;
                    gi = norm_sample(rhs / v, sqrt(vare_ / v));
                    gi_ = oldgi - gi;
                    daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);
                    gi_ *= -1;
                    daxpy_(&n, &gi_, dxi, &inc, u.memptr(), &inc);
                    g[i] = gi;
                    // sumvg += (gi * gi) * vx[i];
                }
                break;
            case 3:
                logpi = log(pi);
                s[0] = logpi[0];
                // sumvg = 0;
                for(int i = 0; i < m; i++){
                    dxi = X.colptr(i);
                    xx = xpx[i];
                    oldgi = g[i];
                    varg = (oldgi * oldgi + s2varg_ * dfvara_) / (dfvara_ + 1);
                    varg = invchisq_sample(1 + dfvara_, varg);
                    rhs = ddot_(&n, dxi, &inc, dyadj, &inc);
                    if(oldgi)   rhs += xx * oldgi;
                    lhs = xx / vare_;
                    logdetV = log(varg * lhs + 1);
                    uhat = rhs / (xx + vare_ / varg);
                    s[1] = -0.5 * (logdetV - (rhs * uhat / vare_)) + logpi[1];
                    acceptProb = 1 / sum(exp(s - s[0]));
                    indistflag = (uniform_sample()) < acceptProb ? 0 : 1;
                    snptracker[i] = indistflag;
                    if(indistflag){
                        v = xx + vare_ / varg;
                        gi = norm_sample(rhs / v, sqrt(vare_ / v));
                        gi_ = oldgi - gi;
                        daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);
                        gi_ *= -1;
                        daxpy_(&n, &gi_, dxi, &inc, u.memptr(), &inc);
                        // sumvg += (gi * gi) * vx[i];
                    }else{
                        gi = 0;
                        if(oldgi){
                            gi_ = oldgi;
                            daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);
                            gi_ *= -1;
                            daxpy_(&n, &gi_, dxi, &inc, u.memptr(), &inc);
                        }
                    }
                    g[i] = gi;
                }
                fold_snp_num[1] = sum(snptracker);
                fold_snp_num[0] = m - fold_snp_num[1];
                NnzSnp = fold_snp_num[1];
                if(!fixpi)  pi = rdirichlet_sampler(n_fold, (fold_snp_num + 1));
                break;
            case 4:
                logpi = log(pi);
                s[0] = logpi[0];
                // sumvg = 0;
                vargi = 0;
                for(int i = 0; i < m; i++){
                    dxi = X.colptr(i);
                    xx = xpx[i];
                    oldgi = g[i];
                    rhs = ddot_(&n, dxi, &inc, dyadj, &inc);
                    if(oldgi)   rhs += xx * oldgi;
                    lhs = xx / vare_;
                    logdetV = log(varg * lhs + 1);
                    uhat = rhs / (xx + vare_ / varg);
                    s[1] = -0.5 * (logdetV - (rhs * uhat / vare_)) + logpi[1];
                    acceptProb = 1 / sum(exp(s - s[0]));
                    indistflag = (uniform_sample()) < acceptProb ? 0 : 1;
                    snptracker[i] = indistflag;
                    if(indistflag){
                        v = xx + vare_ / varg;
                        gi = norm_sample(rhs / v, sqrt(vare_ / v));
                        gi_ = oldgi - gi;
                        daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);
                        gi_ *= -1;
                        daxpy_(&n, &gi_, dxi, &inc, u.memptr(), &inc);
                        // sumvg += (gi * gi) * vx[i];
                        vargi += (gi * gi);
                    }else{
                        gi = 0;
                        if(oldgi){
                            gi_ = oldgi;
                            daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);
                            gi_ *= -1;
                            daxpy_(&n, &gi_, dxi, &inc, u.memptr(), &inc);
                        }
                    }
                    g[i] = gi;
                }
                fold_snp_num[1] = sum(snptracker);
                fold_snp_num[0] = m - fold_snp_num[1];
                NnzSnp = fold_snp_num[1];
                varg = (vargi + s2varg_ * dfvara_) / (dfvara_ + NnzSnp);
                varg = invchisq_sample(NnzSnp + dfvara_, varg);
                if(!fixpi)  pi = rdirichlet_sampler(n_fold, (fold_snp_num + 1));
                break;
            case 5:
                for(int i = 0; i < m; i++){
                    dxi = X.colptr(i);
                    xx = xpx[i];
                    oldgi = g[i];
                    rhs = ddot_(&n, dxi, &inc, dyadj, &inc);
                    rhs += xx * oldgi;
                    v = xx + 1 / vargL[i];
                    gi = norm_sample(rhs / v, sqrt(vare_ / v));
                    if(abs(gi) < 1e-6)  gi = 1e-6;
                    vargi = 1 / rinvgaussian_sample(sqrt(vare_) * lambda / abs(gi), lambda2);
                    if(vargi >= 0)  vargL[i] = vargi;
                    gi_ = oldgi - gi;
                    daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);
                    gi_ *= -1;
                    daxpy_(&n, &gi_, dxi, &inc, u.memptr(), &inc);
                    g[i] = gi;
                    // sumvg += (gi * gi) * vx[i];
                }
                shape = shape0 + m;
                rate = rate0 + sum(vargL) / 2;
                lambda2 = gamma_sample(shape, 1 / rate);
                lambda = sqrt(lambda2);
                break;
            case 6:
                logpi = log(pi);
                s[0] = logpi[0];
                // sumvg = 0;
                varg = 0;
                for(int j = 1; j < n_fold; j++){
                    vare_vara_fold[j] = vare_ / vara_fold[j];
                }
                for(int i = 0; i < m; i++){
                    dxi = X.colptr(i);
                    xx = xpx[i];
                    oldgi = g[i];
                    rhs = ddot_(&n, dxi, &inc, dyadj, &inc);
                    if(oldgi){rhs += xx * oldgi;}
                    lhs = xx / vare_;
                    for(int j = 1; j < n_fold; j++){
                        logdetV = log(vara_fold[j] * lhs + 1);
                        uhat = rhs / (xx + vare_vara_fold[j]);
                        s[j] = -0.5 * (logdetV - (rhs * uhat / vare_)) + logpi[j];
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
                    if(indistflag){
                        v = xx + vare_vara_fold[indistflag];
                        gi = norm_sample(rhs / v, sqrt(vare_ / v));
                        gi_ = oldgi - gi;
                        daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);
                        gi_ *= -1;
                        daxpy_(&n, &gi_, dxi, &inc, u.memptr(), &inc);
                        // sumvg += (gi * gi) * vx[i];
                        varg += (gi * gi / fold_[indistflag]);
                    }else{
                        gi = 0;
                        if(oldgi){
                            gi_ = oldgi;
                            daxpy_(&n, &gi_, dxi, &inc, dyadj, &inc);
                            gi_ *= -1;
                            daxpy_(&n, &gi_, dxi, &inc, u.memptr(), &inc);
                        }
                    }
                    g[i] = gi;
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
                if(!fixpi)  pi = rdirichlet_sampler(n_fold, fold_snp_num + 1);
                break;
        }

        // genetic variance
        vara_ = var(u);
    
        // sample residual variance from inv-chisq distribution
        vare_ = (ddot_(&n, dyadj, &inc, dyadj, &inc) + s2vare_ * dfvare_) / (n + dfvare_);
        vare_ = invchisq_sample(n + dfvare_, vare_);

        if(iter >= nburn){
            count++;
            mu_store[iter - nburn] = mu;
            // sumvg_store[iter - nburn] = sumvg;
            if(!fixpi)  pi_store.row(iter - nburn) = pi.t();
            vara_store[iter - nburn] = vara_;
            vare_store[iter - nburn] = vare_;
            daxpy_(&m, &doc, g.memptr(), &inc, g_store.memptr(), &inc);
            // g_store += g;
            if(nc){
                beta_store.row(iter - nburn) = beta.t();
            }
            double vt = vara_ + vare_;
            if(nr){
                vt += sum(vr);
                vr_store.row(iter - nburn) = vr.t();
                for(int i = 0; i < nr; i++){
                    estR_store[i] += estR[i];
                }
            }
            if(ne){
                veps_store[iter - nburn] = veps;
                epsl_J_beta_store[iter - nburn] = epsl_J_beta;
                epsl_estR_store += epsl_estR;
            }
            hsq_store[iter - nburn] = vara_ / (vt);

            if(!snptracker.is_empty()){
                for(int i = 0; i < m; i++){
                    if(snptracker[i]){
                        nzrate[i] += 1;
                    }
                }
            }

            // record pve for each window
            if(WPPA){
                double vara_snp = vara_;
                if(ne){
                    vec u_minus = u - epsl_J_beta * epsl_y_J_;
                    u_minus.tail(ne) -= epsl_Z * epsl_estR;
                    vara_snp = var(u_minus);
                }
                for(int w = 0; w < nw; w++){
                    wu.fill(0);
                    windxi = windx[w];
                    for(int i = 0; i < windxi.size(); i++){
                        ni = windxi[i];
                        gi_ = g[ni];
                        if(gi_){
                            dxi = X.colptr(ni);
                            // dxi = pX+(long long)ni*n;
                            daxpy_(&n, &gi_, dxi, &inc, wu.memptr(), &inc);
                        }
                        // if(g[windxi[i]]){wu += X(_, windxi[i]) * g[windxi[i]];}
                    }
                    varw = var(wu);
                    wgvei[w] += (varw / vara_snp);
                    if((varw / vara_snp) >= wppa){
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
                for(int i = 0; i < n_fold; i++){
                    Rcpp::Rcout << std::fixed << pi[i] << " ";
                }
                if(model == "BayesL")    Rcpp::Rcout << std::fixed << lambda << " ";
                double vt = vara_ + vare_;
                for(int i = 0; i < nr; i++){
                    vt += vr[i];
                    Rcpp::Rcout << std::fixed << vr[i] << " ";
                }
                // Rcpp::Rcout << std::fixed << sumvg << " ";
                if(ne)  Rcpp::Rcout << std::fixed << veps << " ";
                Rcpp::Rcout << std::fixed << vara_ << " ";
                Rcpp::Rcout << std::fixed << vare_ << " ";
                Rcpp::Rcout << std::fixed << (vara_ / vt) << " ";
                Rprintf("%02dh%02dm%02ds \n", hor, min, sec);
            }
        }
    }
    if(nzrate.is_empty()){
        nzrate.ones(m);
    }else{
        nzrate = nzrate / count;
    }
    g = g_store / count;
    if(nc){
        beta = conv_to<vec>::from(mean(beta_store));
        betase = conv_to<vec>::from(stddev(beta_store));
    }
    if(ne){
        veps = mean(veps_store);
        vepsse = stddev(veps_store);
        epsl_J_beta = mean(epsl_J_beta_store);
        epsl_estR = epsl_estR_store / count;
    }
    List listr(2);
    if(nr){
        vr = conv_to<vec>::from(mean(vr_store));
        vrsd = conv_to<vec>::from(stddev(vr_store));
        vector<string> r_levers;
        vector<double> estr;
        for(int i = 0; i < nr; i++){
            for(int j = 0; j < (estR_store[i]).n_elem; j++){
                r_levers.push_back(Z_levels[i][j]);
                estr.push_back(estR_store[i][j] / count);
            }
        }
        listr[0] = wrap(r_levers.begin(), r_levers.end());
        listr[1] = wrap(estr.begin(), estr.end());
    }
    DataFrame r = listr;
    Rcpp::CharacterVector names(2);
    names[0] = "Levels";
    names[1] = "Estimation";
    if(nr) r.attr("names") = names;

    if(WPPA){
        wppai = wppai / count;
        wgvei = wgvei / count;
    }
    double Mu = mean(mu_store);
    double musd = stddev(mu_store);
    vec pise;
    if(!fixpi){
        pi= conv_to<vec>::from(mean(pi_store));
        pise = conv_to<vec>::from(stddev(pi_store));
    }else{
        pise = zeros(pi.n_elem);
    }
    // sumvg = mean(sumvg_store);
    // double sumvgsd = stddev(sumvg_store);
    vara_ = mean(vara_store);
    double varasd = stddev(vara_store);
    vare_ = mean(vare_store);
    double varesd = stddev(vare_store);
    hsq = mean(hsq_store);
    double hsqsd = stddev(hsq_store);

    if(verbose){
        Rcpp::Rcout << "Posterior parameters:" << std::endl;
        Rcpp::Rcout << "    Mu " << std::fixed << Mu << "±" << std::fixed << musd << std::endl;
        for(int i = 0; i < pi.n_elem; i++){
            if(i == 0){
                Rcpp::Rcout << "    π for markers in zero effect size " << std::fixed << pi[i] << "±" << std::fixed << pise[i] << endl;
            }else{
                if(i == 1){
                    Rcpp::Rcout << "    π for markers in non-zero effect size " << std::fixed << pi[i] << "±" << std::fixed << pise[i];
                }else{
                    Rcpp::Rcout << std::fixed << pi[i] << "±" << std::fixed << pise[i];
                }
                if(i != (pi.n_elem - 1))  Rcpp::Rcout << ", ";
            }
        }
        Rcpp::Rcout << std::endl;
        if(nc){
            Rcpp::Rcout << "    Effects for covariates ";
            for(int i = 0; i < nc; i++){
                Rcpp::Rcout << std::fixed << beta[i] << "±" << std::fixed << betase[i];
                if(i != (nc - 1))  Rcpp::Rcout << ", ";
            }
            Rcpp::Rcout << std::endl;
        }
        if(nr){
            Rcpp::Rcout << "    Environmental var ";
            for(int i = 0; i < nr; i++){
                Rcpp::Rcout << std::fixed << vr[i] << "±" << std::fixed << vrsd[i];
                if(i != (nr - 1))  Rcpp::Rcout << ", ";
            }
            Rcpp::Rcout << std::endl;
        }
        // Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << "±" << std::fixed << sumvgsd << std::endl;
        if(ne)  Rcpp::Rcout << "    ε var " << std::fixed << veps << "±" << std::fixed << vepsse << std::endl;
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
    if(WPPA){
        return List::create(
            Named("eps_J") = epsl_J_beta,
            Named("eps_R") = epsl_estR,
            Named("mu") = Mu, 
            Named("pi") = pi, 
            Named("beta") = beta,
            Named("r") = r, 
            Named("vr") = vr, 
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("alpha") = g,
            Named("modfreq") = nzrate,
            Named("wppa") = wppai,
            Named("wgve") = wgvei
        );
    }else{
        return List::create(
            Named("eps_J") = epsl_J_beta,
            Named("eps_R") = epsl_estR,
            Named("mu") = Mu,
            Named("pi") = pi, 
            Named("beta") = beta,
            Named("r") = r, 
            Named("vr") = vr, 
            Named("vg") = vara_, 
            Named("ve") = vare_,
            Named("alpha") = g,
            Named("modfreq") = nzrate
        );
    }
}
