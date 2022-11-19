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

template <typename T>
Rcpp::List Bayes(
    arma::vec &y,
    arma::mat &X,
    std::string model,
    arma::vec Pi,
    T &K,
    arma::uvec &K_index,
    const Nullable<arma::mat> C = R_NilValue,
    const Nullable<CharacterMatrix> R = R_NilValue,
    const Nullable<arma::vec> fold = R_NilValue,
    const int niter = 50000,
    const int nburn = 20000,
    const int thin = 5,
    const Nullable<arma::vec> epsl_y_J = R_NilValue,
    const Nullable<arma::sp_mat> epsl_Gi = R_NilValue,
    const Nullable<arma::uvec> epsl_index = R_NilValue,
    const Nullable<double> dfvr = R_NilValue,
    const Nullable<double> s2vr = R_NilValue,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
    const Nullable<arma::uvec> windindx = R_NilValue,
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
    int model_index = (model == "BayesRR" ? 1 : (model == "BayesA" ? 2 : (model == "BayesB" || model == "BayesBpi" ? 3 : (model == "BayesC" || model == "BayesCpi" || model == "BSLMM" ? 4 : (model == "BayesL" ? 5 : 6)))));
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

    int n = y.n_elem;
    int m = X.n_cols;
    double vary = var(y);
    double h2 = 0.5;

    int n_records = (niter - nburn) / thin;

    double* dci;
    int nc = 0;
    mat C_;
    vec beta;
    vec betasd;
    mat beta_store;
    if(C.isNotNull()){
        C_ = as<arma::mat>(C);
        if(C_.n_rows != X.n_rows)   throw Rcpp::exception("Number of individuals does not match for covariates.");
        if(C_.has_nan())  throw Rcpp::exception("Individuals with phenotypic value should not have missing covariates.");
        nc = C_.n_cols;
        beta.resize(nc);
        beta_store.resize(nc, n_records);
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
    vec vrsd;
    vec vrtmp;
    sp_mat r_LHS;
    vec estR;
    vec estR_tmp;
    vec R_idx;
    vec r_RHS;
    vec diff;
    int qr;
    double dfr;
    double s2r;
    if(dfvr.isNotNull()){
        dfr = as<double>(dfvr);
    }else{
        dfr = -1;
    }
    if(s2vr.isNotNull()){
        s2r = as<double>(s2vr);
    }else{
        s2r = 0;
    }
    mat vr_store;
    mat estR_store;
    vector<sp_mat> Z;
    vector<sp_mat> ZZ;
    vector< vector<string> > Z_levels;
    if(R.isNotNull()){
        CharacterMatrix R_ = as<CharacterMatrix>(R);
        if(R_.nrow() != X.n_rows)   throw Rcpp::exception("Number of individuals does not match for environmental random effects.");
        if(xhasNA(R_))  throw Rcpp::exception("Individuals with phenotypic value should not have missing environmental random effects.");
        nr = R_.ncol();
        vr.resize(nr);
        vrtmp.resize(nr);
        R_idx.resize(nr + 1); R_idx[0] = -1;
        vr_store.resize(nr, n_records);
        vrtmp.fill(vary * (1 - h2) / (nr + 1));
        // s2r = vrtmp[0] * (dfr - 2) / dfr;
        int n_levels = 0;
        for(int i = 0; i < nr; i++){
            CharacterVector Ri = R_(_, i);
            List Z_info = makeZ(Ri);
            sp_mat z = Z_info[0];
            vector<string> z_levels = Z_info[1];
            Z.push_back(z);
            ZZ.push_back(z.t() * z);
            Z_levels.push_back(z_levels);
            n_levels += z_levels.size();
            R_idx[i + 1] = n_levels - 1;
        }
        estR.zeros(n_levels);
        estR_store.zeros(n_levels, n_records);
    }

    int nk = 0;
    double va, vb;
    double vbtmp;
    double vasd, vbsd;
    vec va_store;
    vec vb_store;
    sp_mat k_Z, k_ZZ;
    vec k_estR;
    vec k_estR_store;
    T k_LHS;
    vec k_estR_tmp;
    vec k_RHS;
    // mat ghat;
    int qk = 0;
    if(!K_index.is_empty()){
        if(K.is_empty())    throw Rcpp::exception("variance-covariance matrix should be provided for K term.");
        if(K.n_cols != K.n_rows) throw Rcpp::exception("variance-covariance matrix should be in square.");
        nk = K.n_cols;
        K_index -= 1;
        va_store.resize(n_records);
        vb_store.resize(n_records);
        // ghat.resize(m, n_records);
        k_Z.resize(K_index.n_elem, nk);
        k_estR.zeros(nk);
        k_estR_store.zeros(nk);
        k_estR_tmp.zeros(nk);
        qk = nk;
        for(int i = 0; i < K_index.n_elem; i++){k_Z(i, K_index[i]) = 1.0;}
        k_ZZ = k_Z.t() * k_Z;
    }

    int ne = 0;
    double veps;
    double vepstmp;
    double vepssd;
    double JtJ;
    vec veps_store;
    sp_mat epsl_Gi_, epsl_Z, epsl_ZZ;
    uvec epsl_index_;
    vec epsl_estR;
    mat epsl_estR_store;
    vec epsl_y_J_;
    sp_mat epsl_LHS;
    vec epsl_yadj;
    vec epsl_estR_tmp;
    vec epsl_RHS;
    int qe = 0;
    double epsl_J_beta = 0;
    double epsl_J_betasd = 0;
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
        veps_store.resize(n_records);
        epsl_J_beta_store.resize(n_records);
        epsl_Z.resize(ne, epsl_Gi_.n_cols);
        epsl_estR.zeros(epsl_Gi_.n_cols);
        epsl_estR_tmp.zeros(epsl_Gi_.n_cols);
        epsl_estR_store.resize(epsl_Gi_.n_cols, n_records);
        epsl_yadj.resize(ne);
        qe = epsl_Gi_.n_cols;
        for(int i = 0; i < ne; i++){epsl_Z(i, epsl_index_[i]) = 1.0;}
        epsl_ZZ = epsl_Z.t() * epsl_Z;
    }

    int count = 0;
    int nzct = 0;
    int inc = 1;
    int n_fold = fold_.n_elem;
    int NnzSnp, indistflag;
    double doc = 1.0;
    double xx, oldgi, gi, gi_, rhs, lhs, logdetV, acceptProb, uhat, v;
    double vara_, dfvara_, s2vara_, vare_, dfvare_, s2vare_, vargi, hsq, s2varg_;
    // double sumvg;
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
    vec g = zeros(m);
    mat g_store = zeros(m, n_records);
    vec u = zeros(n);
    vec xpx = zeros(m);
    vec vx = zeros(m);
    vec mu_store = zeros(n_records);
    // vec sumvg_store = zeros(niter - nburn + 1);
    vec vara_store = zeros(n_records);
    vec vare_store = zeros(n_records);
    vec hsq_store = zeros(n_records);
    mat pi_store;
    pi_store.resize(n_fold, n_records);

    #pragma omp parallel for
    for(int i = 0; i < m; i++){
        vec Xi = X.col(i);
        xpx[i] = sum(square(Xi));
        vx[i] = var(Xi);
    }
    double sumvx = sum(vx);
    int nvar0 = sum(vx == 0);

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
    vbtmp = vara_;
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
    double varg = vara_ / ((1 - Pi[0]) * sumvx);
    s2varg_ = s2vara_ / ((1 - Pi[0]) * sumvx);
    if(s2ve.isNotNull()){
        s2vare_ = as<double>(s2ve);
    }else{
        s2vare_ = 0;
    }
    if(niter < nburn){
        throw Rcpp::exception("Number of total iteration ('niter') shold be larger than burn-in ('nburn').");
    }
    double R2 = (dfvara_ - 2) / dfvara_;
    double lambda2 = 2 * (1 - R2) / (R2) * sumvx;
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
    vec vara_fold = (vara_ / ((1 - Pi[0]) * sumvx)) * fold_;
    vec vare_vara_fold = zeros(n_fold);

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
        if(epsl_index.isNotNull()){
            Rcpp::Rcout << "    Observations with genotype " << (n - ne) << std::endl;
            Rcpp::Rcout << "    Observations with imputed genotype " << ne << std::endl;
        }
        Rcpp::Rcout << "    Number of covariates " << (nc + 1) << std::endl;
        Rcpp::Rcout << "    Number of envir-random effects " << nr << std::endl;
        Rcpp::Rcout << "    Number of markers " << m << std::endl;
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
        Rcpp::Rcout << "    Frequency of collecting " << thin << std::endl;
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
            Rcpp::Rcout << "    Number of windows for GWAS analysis " << nw << std::endl;
        }
        Rcpp::Rcout << "MCMC started: " << std::endl;
        Rcpp::Rcout << " Iter" << "  ";
        Rcpp::Rcout << "NumNZSnp" << "  ";
        for(int i = 0; i < n_fold; i++){
            Rcpp::Rcout << "π" << i + 1 << "  ";
        }
        if(model == "BayesL")    Rcpp::Rcout << "Lambda" << "  ";
        if(model == "BSLMM")    Rcpp::Rcout << "Va" << "  " << "Vb" << "  ";
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
        mu_ = - norm_sample(sum(yadj) / n, sqrt(vare_ / n));
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
            estR_tmp = estR.subvec(R_idx[i] + 1, R_idx[i + 1]);
            r_RHS = Z[i].t() * yadj;
            r_RHS += ZZ[i] * estR_tmp;
            Gibbs(r_LHS, estR_tmp, r_RHS, vare_);
            diff = Z[i] * (estR.subvec(R_idx[i] + 1, R_idx[i + 1]) - estR_tmp);
            daxpy_(&n, &doc, diff.memptr(), &inc, dyadj, &inc);
            qr = estR_tmp.n_elem;
            vrtmp[i] = (ddot_(&qr, estR_tmp.memptr(), &inc, estR_tmp.memptr(), &inc) + s2r * dfr) / (qr + dfr);
            vrtmp[i] = invchisq_sample(qr + dfr, vrtmp[i]);
            vr[i] = var(estR_tmp);
            // vr[i] = vrtmp[i];
            estR.subvec(R_idx[i] + 1, R_idx[i + 1]) = estR_tmp;
        }

        if(nk){
            k_LHS = (K * (vare_ / vbtmp));
            k_LHS += k_ZZ;
            k_RHS = k_Z.t() * yadj;
            vec zzk = k_ZZ * k_estR_tmp;
            daxpy_(&qk, &doc, zzk.memptr(), &inc, k_RHS.memptr(), &inc);
            Gibbs(k_LHS, k_estR_tmp, k_RHS, vare_);
            gi_ = -1;
            daxpy_(&qk, &gi_, k_estR_tmp.memptr(), &inc, k_estR.memptr(), &inc);
            diff = k_Z * k_estR;
            daxpy_(&n, &doc, diff.memptr(), &inc, dyadj, &inc);
            daxpy_(&n, &gi_, diff.memptr(), &inc, u.memptr(), &inc);
            vbtmp = as_scalar(k_estR_tmp.t() * K * k_estR_tmp);
            vbtmp += s2vara_ * dfvara_;
            vbtmp /= (dfvara_ + qk);
            vbtmp = invchisq_sample(qk + dfvara_, vbtmp);
            // vb = var(k_estR_tmp) / sumvx;
            vb = vbtmp;
            k_estR = k_estR_tmp;
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
            epsl_RHS = epsl_Z.t() * epsl_yadj;
            epsl_RHS += epsl_ZZ * epsl_estR_tmp;
            Gibbs(epsl_LHS, epsl_estR_tmp, epsl_RHS, vare_);
            // epsl_estR_tmp = spsolve(epsl_LHS, epsl_RHS, "lapack");
            gi_ = -1;
            daxpy_(&qe, &gi_, epsl_estR_tmp.memptr(), &inc, epsl_estR.memptr(), &inc);
            epsl_yadj = epsl_Z * epsl_estR;
            yadj.tail(ne) += epsl_yadj;
            u.tail(ne) -= epsl_yadj;
            vepstmp = as_scalar(epsl_estR_tmp.t() * epsl_Gi_ * epsl_estR_tmp);
            vepstmp += s2vara_ * dfvara_;
            vepstmp /= (dfvara_ + qe);
            vepstmp = invchisq_sample(qe + dfvara_, vepstmp);
            epsl_estR = epsl_estR_tmp;
            // veps = var(epsl_estR);
            veps = vepstmp;
        }
        
        switch(model_index){
            case 1:
                for(int i = 0; i < m; i++){
                    if(!vx[i])   continue;
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
                varg = (ddot_(&m, g.memptr(), &inc, g.memptr(), &inc) + s2varg_ * dfvara_) / (dfvara_ + m - nvar0);
                varg = invchisq_sample(m - nvar0 + dfvara_, varg);
                // sumvg = sum(vx * varg);
                break;
            case 2:
                for(int i = 0; i < m; i++){
                    if(!vx[i])   continue;
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
                logpi = log(Pi);
                s[0] = logpi[0];
                // sumvg = 0;
                for(int i = 0; i < m; i++){
                    if(!vx[i])   continue;
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
                fold_snp_num[0] = m - nvar0 - fold_snp_num[1];
                NnzSnp = fold_snp_num[1];
                if(!fixpi)  Pi = rdirichlet_sample(n_fold, (fold_snp_num + 1));
                break;
            case 4:
                logpi = log(Pi);
                s[0] = logpi[0];
                // sumvg = 0;
                vargi = 0;
                for(int i = 0; i < m; i++){
                    if(!vx[i])   continue;
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
                fold_snp_num[0] = m - nvar0 - fold_snp_num[1];
                NnzSnp = fold_snp_num[1];
                varg = (vargi + s2varg_ * dfvara_) / (dfvara_ + NnzSnp);
                varg = invchisq_sample(NnzSnp + dfvara_, varg);
                if(nk) va = varg;
                if(!fixpi)  Pi = rdirichlet_sample(n_fold, (fold_snp_num + 1));
                break;
            case 5:
                for(int i = 0; i < m; i++){
                    if(!vx[i])   continue;
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
                shape = shape0 + m - nvar0;
                rate = rate0 + sum(vargL) / 2;
                lambda2 = gamma_sample(shape, 1 / rate);
                lambda = sqrt(lambda2);
                break;
            case 6:
                logpi = log(Pi);
                s[0] = logpi[0];
                // sumvg = 0;
                varg = 0;
                for(int j = 1; j < n_fold; j++){
                    vare_vara_fold[j] = vare_ / vara_fold[j];
                }
                for(int i = 0; i < m; i++){
                    if(!vx[i])   continue;
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
                fold_snp_num[0] -= nvar0;
                if(!fixpi)  Pi = rdirichlet_sample(n_fold, (fold_snp_num + 1));
                break;
        }

        // genetic variance
        vara_ = var(u);
        // vara_ = invchisq_sample(n + dfvara_, vara_);
    
        // sample residual variance from inv-chisq distribution
        vare_ = (ddot_(&n, dyadj, &inc, dyadj, &inc) + s2vare_ * dfvare_) / (n + dfvare_);
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

        // collect the parameters to obtain posterior estimation
        if(iter >= nburn && (iter + 1 - nburn) % thin == 0){
            mu_store[count] = mu;
            if(!fixpi)  pi_store.col(count) = Pi;
            vara_store[count] = vara_;
            vare_store[count] = vare_;

            if(nk){
                va_store[count] = va;
                vb_store[count] = vb;
                // k_estR_store += k_estR;
                daxpy_(&qk, &doc, k_estR.memptr(), &inc, k_estR_store.memptr(), &inc);
            }

            g_store.col(count) = g;
            // daxpy_(&m, &doc, g.memptr(), &inc, g_store.memptr(), &inc);
            // g_store += g;
            if(nc){
                beta_store.col(count) = beta;
            }
            double vt = vara_ + vare_;
            if(nr){
                vt += sum(vr);
                vr_store.col(count) = vr;
                estR_store.col(count) = estR;
            }
            if(ne){
                veps_store[count] = veps;
                epsl_J_beta_store[count] = epsl_J_beta;
                epsl_estR_store.col(count) = epsl_estR;
                // daxpy_(&qe, &doc, epsl_estR.memptr(), &inc, epsl_estR_store.memptr(), &inc);
            }
            hsq_store[count] = vara_ / vt;

            count++;
        }

        if(verbose && (iter + 1) % outfreq == 0){

            // print iteration details
            if(verbose){
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
                if(model == "BSLMM")    Rcpp::Rcout << std::fixed << va << " " << std::fixed << vb << " ";
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

        if(count == n_records)  break;
    }

    List results;
    List MCMCsample;

    if(nr){
        vr = conv_to<vec>::from(mean(vr_store, 1));
        vrsd = conv_to<vec>::from(stddev(vr_store, 0, 1));
        results["Vr"] = vr;
        MCMCsample["Vr"] = vr_store;
    }
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

    double Mu = mean(mu_store);
    vec e = y - (Mu * one);
    double musd = stddev(mu_store);
    results["mu"] = Mu;
    MCMCsample["mu"] = mu_store.t();

    if(nc){
        beta = conv_to<vec>::from(mean(beta_store, 1));
        betasd = conv_to<vec>::from(stddev(beta_store, 0, 1));
        e -= (C_ * beta);
        results["beta"] = beta;
        MCMCsample["beta"] = beta_store;
    }

    if(nk){
        vec ghat = X.t() * (k_Z * (K * k_estR_store  / count)) / sumvx;
        g_store.each_col() += ghat;
        va = mean(va_store);
        vasd = stddev(va_store);
        vb = mean(vb_store);
        vbsd = stddev(vb_store);
    }
    g = conv_to<vec>::from(mean(g_store, 1));
    e -= (X * g);
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

    if(ne){
        veps = mean(veps_store);
        vepssd = stddev(veps_store);
        epsl_J_beta = mean(epsl_J_beta_store);
        epsl_J_betasd = stddev(epsl_J_beta_store);
        epsl_estR = conv_to<vec>::from(mean(epsl_estR_store, 1));
        e -= (epsl_J_beta * epsl_y_J_);
        e.tail(ne) -= (epsl_Z * epsl_estR);
        results["Veps"] = veps;
        results["J"] = epsl_J_beta;
        results["epsilon"] = epsl_estR;
        MCMCsample["Veps"] = veps_store.t();
        MCMCsample["J"] = epsl_J_beta_store.t();
        MCMCsample["epsilon"] = epsl_estR_store;
    }

    if(nr){
        List listr(2);
        estR = conv_to<vec>::from(mean(estR_store, 1));
        vector<string> r_levers;
        for(int i = 0; i < nr; i++){
            for(int j = 0; j < (Z_levels[i]).size(); j++){
                r_levers.push_back(Z_levels[i][j]);
            }
            e -= (Z[i] * estR.subvec(R_idx[i] + 1, R_idx[i + 1]));
        }
        listr[0] = wrap(r_levers.begin(), r_levers.end());
        listr[1] = wrap(estR);
        DataFrame r = listr;
        Rcpp::CharacterVector names(2);
        names[0] = "Levels";
        names[1] = "Estimation";
        r.attr("names") = names;
        results["r"] = r;
        MCMCsample["r"] = estR_store;
    }
    results["g"] = u;
    results["e"] = e;

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
        Rcpp::Rcout << "    Mu " << std::fixed << Mu << "±" << std::fixed << musd << std::endl;
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
        if(nc){
            Rcpp::Rcout << "    Effects for covariates ";
            for(int i = 0; i < nc; i++){
                Rcpp::Rcout << std::fixed << beta[i] << "±" << std::fixed << betasd[i];
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
        if(nk)  Rcpp::Rcout << "    Va: " << std::fixed << va << "±" << std::fixed << vasd << ", Vb: " << std::fixed << vb << "±" << std::fixed << vbsd << std::endl;
        // Rcpp::Rcout << "    SigmaSq " << std::fixed << sumvg << "±" << std::fixed << sumvgsd << std::endl;
        if(ne){
            Rcpp::Rcout << "    J " << std::fixed << epsl_J_beta << "±" << std::fixed << epsl_J_betasd << std::endl;
            Rcpp::Rcout << "    ε var " << std::fixed << veps << "±" << std::fixed << vepssd << std::endl;
        }
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

// [[Rcpp::export]]
Rcpp::List BayesK(
    arma::vec &y,
    arma::mat &X,
    std::string model,
    arma::vec Pi,
    SEXP &K,
    arma::uvec &K_index,
    const Nullable<arma::mat> C = R_NilValue,
    const Nullable<CharacterMatrix> R = R_NilValue,
    const Nullable<arma::vec> fold = R_NilValue,
    const int niter = 50000,
    const int nburn = 20000,
    const int thin = 5,
    const Nullable<arma::vec> epsl_y_J = R_NilValue,
    const Nullable<arma::sp_mat> epsl_Gi = R_NilValue,
    const Nullable<arma::uvec> epsl_index = R_NilValue,
    const Nullable<double> dfvr = R_NilValue,
    const Nullable<double> s2vr = R_NilValue,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
    const Nullable<arma::uvec> windindx = R_NilValue,
    const int outfreq = 100,
    const int threads = 0,
    const bool verbose = true
){
    if(Rf_inherits(K, "dgCMatrix")){
        arma::sp_mat K_ = Rcpp::as<arma::sp_mat>(K);
        return Bayes(y, X, model, Pi, K_, K_index, C, R, fold, niter, nburn, thin, epsl_y_J, epsl_Gi, epsl_index, dfvr, s2vr, vg, dfvg, s2vg, ve, dfve, s2ve, windindx, outfreq, threads, verbose);
    }else{
        arma::mat K_ = Rcpp::as<arma::mat>(K);
        return Bayes(y, X, model, Pi, K_, K_index, C, R, fold, niter, nburn, thin, epsl_y_J, epsl_Gi, epsl_index, dfvr, s2vr, vg, dfvg, s2vg, ve, dfve, s2ve, windindx, outfreq, threads, verbose);
    }
}

// [[Rcpp::export]]
Rcpp::List Bayes(
    arma::vec &y,
    arma::mat &X,
    std::string model,
    arma::vec Pi,
    const Nullable<arma::mat> C = R_NilValue,
    const Nullable<CharacterMatrix> R = R_NilValue,
    const Nullable<arma::vec> fold = R_NilValue,
    const int niter = 50000,
    const int nburn = 20000,
    const int thin = 5,
    const Nullable<arma::vec> epsl_y_J = R_NilValue,
    const Nullable<arma::sp_mat> epsl_Gi = R_NilValue,
    const Nullable<arma::uvec> epsl_index = R_NilValue,
    const Nullable<double> dfvr = R_NilValue,
    const Nullable<double> s2vr = R_NilValue,
    const Nullable<double> vg = R_NilValue,
    const Nullable<double> dfvg = R_NilValue,
    const Nullable<double> s2vg = R_NilValue,
    const Nullable<double> ve = R_NilValue,
    const Nullable<double> dfve = R_NilValue,
    const Nullable<double> s2ve = R_NilValue,
    const Nullable<arma::uvec> windindx = R_NilValue,
    const int outfreq = 100,
    const int threads = 0,
    const bool verbose = true
){
    arma::mat K;
    arma::uvec K_index;
    return Bayes(y, X, model, Pi, K, K_index, C, R, fold, niter, nburn, thin, epsl_y_J, epsl_Gi, epsl_index, dfvr, s2vr, vg, dfvg, s2vg, ve, dfve, s2ve, windindx, outfreq, threads, verbose);
}
