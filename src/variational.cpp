/* 
 * This file is part of the hibayes (https://github.com/YinLiLin/hibayes).
 * Copyright (c) 2023 Haohao Zhang.
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "hibayes.h"
#include "solver.h"
#include "magic_enum.hpp"

enum class Model { BayesA, BayesB, BayesBpi, BayesC, BayesCpi, BayesL, BayesR, BayesRR, BSLMM };

// [[Rcpp::export]]
Rcpp::List VariationalBayes(
    arma::vec &y,
    arma::mat &X,
    std::string model_str,
    arma::vec Pi,
    const Nullable<arma::mat> C = R_NilValue,
    const Nullable<CharacterMatrix> R = R_NilValue,
    const Nullable<arma::vec> fold = R_NilValue,
    const Nullable<double> dfvr = R_NilValue,
    const Nullable<double> s2vr = R_NilValue,
    const Nullable<double> vg = R_NilValue,     // additive gentic variance
    const Nullable<double> dfvg = R_NilValue,   // Nu of Inv-Chi2 distribution of additive gentic variance
    const Nullable<double> s2vg = R_NilValue,   // S2 of Inv-Chi2 distribution of additive gentic variance
    const Nullable<double> ve = R_NilValue,     // residual variance
    const Nullable<double> dfve = R_NilValue,   // Nu of Inv-Chi2 distribution of residual variance
    const Nullable<double> s2ve = R_NilValue,   // S2 of Inv-Chi2 distribution of residual variance
    const int threads = 0,
    const int max_iteration = 1000,
    const double threshold = 1e-5,
    const bool random_init = false,
    const bool verbose = true
) {
    omp_setup(threads);

    Model model_index = magic_enum::enum_cast<Model>(model_str).value();

    // Check parameters
    if (model_index != Model::BayesR && Pi.n_elem != 2)
        throw Rcpp::exception("length of Pi should be 2, the first value is the proportion of non-effect markers.");
        
    if (sum(Pi) != 1)
        throw Rcpp::exception("sum of Pi should be 1.");
    if (Pi[0] == 1)
        throw Rcpp::exception("all markers have no effect size.");
    for (int i = 0; i < Pi.n_elem; i++) {
        if(Pi[i] < 0 || Pi[i] > 1){
            throw Rcpp::exception("elements of Pi should be at the range of [0, 1]");
        }
    }

    if (dfvg.isNotNull() && as<double>(dfvg) <= 2)
        throw Rcpp::exception("dfvg should not be less than 2.");

    // Check model specific parameters
    vec fold_;
    if (model_index == Model::BayesR) {
        if(fold.isNotNull()) {
            fold_ = as<arma::vec>(fold);
        } else {
            throw Rcpp::exception("'fold' should be provided for BayesR model.");
        }

        if(fold_.n_elem != Pi.n_elem)
            throw Rcpp::exception("length of Pi and fold not equals.");
    }

    // Check input
    if (y.has_nan())
        throw Rcpp::exception("NAs are not allowed in y.");
    if (y.n_elem != X.n_rows)
        throw Rcpp::exception("Number of individuals not equals.");
    
    // Check covariates
    // int nc = 0;             // number of covariates
    // mat C_;                 // covariates matrix
    // if(C.isNotNull()){
    //     C_ = as<arma::mat>(C);
    //     if(C_.n_rows != X.n_rows)
    //         throw Rcpp::exception("Number of individuals does not match for covariates.");
    //     if(C_.has_nan())
    //         throw Rcpp::exception("Individuals with phenotypic value should not have missing covariates.");
    //     nc = C_.n_cols;
    // }

    // Check random effects
    // TODO: check random effects

    // Initialize =================================================================
    uword n = y.n_elem;     // Number of individuals
    uword m = X.n_cols;     // Number of markers
    double vary = var(y);   // Variance of y
    double h2 = 0.5;        // H2 (default 0.5)
    uword i, j, k, p;       // loop variable: i -> individual/other, p -> marker, j -> covariate, k -> random effect

    // Covariates =========================
    double mu_exp = mean(y);    // E(mu)
    double mu_var = 1.0;        // V(mu)
    double mu_exp2 = pow(mu_exp, 2.0) + mu_var;   // E(mu^2)

    // if (nc) {
    //     vec cpc = zeros(nc);    // Vector of dot(c_j * c_j)
    //     // #pragma omp parallel for
    //     for(uword j = 0; j < nc; j++)
    //         cpc[j] = accu(square(C_.col(j)));
    // }

    // Gentic marker ======================
    // Cache feature of markers
    vec xpx = zeros(m);     // Vector of (X_p * X_p)
    vec vx = zeros(m);      // Vector of variance of X_p
    // #pragma omp parallel for
    for(uword p = 0; p < m; p++) {
        auto Xp = X.col(p);
        xpx[p] = accu(Xp % Xp);
        vx[p] = var(Xp);
    }
    double sumvx = sum(vx);

    // parameters of the additive gentic effect distribution
    double dfvg_ = dfvg.isNotNull() ? as<double>(dfvg) : 4;                           // Nu of Inv-Chi2 distribution
    double vg_ = vg.isNotNull()? as<double>(vg) : ((dfvg_ - 2) / dfvg_) * vary * h2;  
    double s2vg_ = s2vg.isNotNull()? as<double>(s2vg) : vg_ * (dfvg_ - 2) / dfvg_;    // S2 of Inv-Chi2 distribution

    // parameters of the marker effect distribution
    double vp_ = vg_ / ((1 - Pi[0]) * sumvx);
    double s2vp_ = s2vg_ / ((1 - Pi[0]) * sumvx);

    vec alpha_exp;      // marker effects expectation
    vec alpha_var;      // marker effects variance
    vec alpha_exp2;     // marker effects expectation^2

    vec alpha_sigma2;   // marker effects ??
    vec alpha_S2;       // marker effects S2 of the Inverse-Chi2 distribution
    
    double sigma2;      // 
    double S2;          // S2 of the Inverse-Chi2 distribution

    vec gamma_exp;      // probability of marker have effect
    vec gamme_exp2;     // probability of marker have effect^2

    // Initialize - model
    switch (model_index) {
    case Model::BayesA:
        alpha_exp = zeros(m);
        alpha_var = zeros(m);
        alpha_exp2 = zeros(m);  // square(alpha_exp) + alpha_var;

        alpha_sigma2 = vec(m, fill::value(1.0 / m));
        alpha_S2 = vec(m, fill::value(m));
        break;
    case Model::BayesB:
    case Model::BayesBpi:
    case Model::BayesC:
        alpha_exp = zeros(m);
        alpha_var = zeros(m);
        alpha_exp2 = zeros(m);

        gamma_exp = vec(m, fill::value(0.5));
        gamme_exp2 = vec(m, fill::value(0.5));      // gamma_exp ^ 2 + gamma_exp * (1 - gamma_exp)
        break;
    case Model::BayesCpi:
    case Model::BayesL:
    case Model::BayesR:
    case Model::BayesRR:
    case Model::BSLMM:
        throw Rcpp::exception("Not implemented yet");
    }
    
    // Residual ===========================
    vec y_res = y - mu_exp;  // residual

    // parameters of the residual distribution
    double dfve_ = dfve.isNotNull()? as<double>(dfve) : -2;     // Nu of Inv-Chi2 distribution     
    double s2ve_ = s2ve.isNotNull()? as<double>(s2ve) : 0;      // S2 of Inv-Chi2 distribution
    double ve_ = ve.isNotNull()? as<double>(ve) : vary * (1 - h2); //TODO: / (nr + 1); 
    
    double tau_exp = 1.0 / ve_;
    double tau_exp_, tau_var_;


    if(verbose){
        Rcpp::Rcout.precision(4);
        Rcpp::Rcout << "Prior parameters:" << std::endl;
        Rcpp::Rcout << "    Model fitted at [" << model_str << "]" << std::endl;
        Rcpp::Rcout << "    Number of observations " << n << std::endl;
        // Rcpp::Rcout << "    Number of covariates " << (nc + 1) << std::endl;
        // Rcpp::Rcout << "    Number of envir-random effects " << nr << std::endl;
        Rcpp::Rcout << "    Number of markers " << m << std::endl;
        for(uword i = 0; i < Pi.n_elem; i++) {
            switch (i) {
                case 0:  Rcpp::Rcout << "    π for markers in zero effect size " << Pi[i] << endl; break;
                case 1:  Rcpp::Rcout << "    π for markers in non-zero effect size " << Pi[i]; break;
                default: Rcpp::Rcout << " " << Pi[i];
            }
        }
        Rcpp::Rcout << endl;
        Rcpp::Rcout << "    Phenotypic var " << std::fixed << vary << std::endl;

        Rcpp::Rcout << "MCMC started: " << std::endl;
        Rcpp::Rcout << " Iter" << "  ";
        Rcpp::Rcout << "NumNZSnp" << "  ";
        Rcpp::Rcout << "Vg" << "  ";
        Rcpp::Rcout << "Ve" << "  ";
        Rcpp::Rcout << "h2" << "  ";
        Rcpp::Rcout << "Timeleft" << std::endl;
    }

    MyTimer timer;
    timer.step("start");
    // Optimization ==============================================================
    for (int iter = 1; iter <= max_iteration; iter++) {
        double vS2;
		// To check convergence
		double Check1 = 0.0;        // sum((new_para - old_para)^2)
        double Check2 = 0.0;        // sum((new_para)^2)

        /* Probability of gamma */
        double sumGammaB2 = 0.0;
        double sumGamma   = 0.0;
	    double gamma_constant;

        double sumGammaB2_ = 0.0;
        double sumGamma_   = 0.0;
        double S2_ = 0.0;

		// For update of residual variance
        double a1, b1;
		double sumVarB = 0.0;

        uvec order = linspace<uvec>(0, m - 1, m);     // order of snp {1 ... m}, used to update SNP Effect in random order

		// Update of alpha
        switch(model_index) {
        case Model::BayesA:
        //#region BayesA
            s2vp_ = pow((dfvg_ - 2) / dfvg_, 2.0) * vary * h2 / ((1 - Pi[0]) * sumvx); 

            vS2 = dfvg_ * s2vp_;

            order = shuffle(order);     // update SNP Effect in random order
            for (uword o = 0; o < m; o++) {
                uword p = order[o];     // index of snp
                double alpha_exp_;      // new value of alpha_exp[p]
                double alpha_exp2_;     // new value of alpha_exp2[p]
                double alpha_var_;      // new value of alpha_var[p]
                double alpha_S2_;       // new value of alpha_S2[p]

                alpha_var_ = 1.0 / (xpx[p] * tau_exp + 1.0 / alpha_S2[p]);
                alpha_exp_  = alpha_var_ * (dot(X.col(p), y_res) + xpx[p] * alpha_exp[p]) * tau_exp;
                alpha_exp2_ = pow(alpha_exp_, 2.0) + alpha_var_;

                /* update residuals */
                y_res   += (X.col(p) * (alpha_exp[p] - alpha_exp_));
                sumVarB += (xpx[p] * alpha_var[p]);     // ??

                Check1 += pow(alpha_exp_ - alpha_exp[p], 2.0);
                Check2 += pow(alpha_exp_, 2.0);

                alpha_exp[p]  = alpha_exp_;
                alpha_exp2[p] = alpha_exp2_;
                alpha_var[p]  = alpha_var_;

                /* update of Sigma2 */
                alpha_sigma2[p] = (alpha_exp2[p] + vS2) / (dfvg_ - 1.0);
                alpha_S2_ = (alpha_exp2[p] + vS2) / (dfvg_ + 1.0);
                Check1 += pow(alpha_S2_ - alpha_S2[p], 2.0); 
                Check2 += pow(alpha_S2_, 2.0);
                alpha_S2[p] = alpha_S2_;
            }
        break;
        //#endregion
        case Model::BayesB:
            throw Rcpp::exception("Not implemented yet");
        break;
        case Model::BayesBpi:
            throw Rcpp::exception("Not implemented yet");
        break;
        case Model::BayesC:
            s2vp_ = pow((dfvg_ - 2) / dfvg_, 2.0) * vary * h2 / ((1 - Pi[0]) * sumvx); 
            vS2 = dfvg_ * s2vg_;

            gamma_constant = 0.5 * R::digamma(0.5 * (dfvg_ + sumGamma)) - 0.5 * log(0.5 * (sumGammaB2 + vS2)) + log(datum::pi);

            order = shuffle(order);     // update SNP Effect in random order

            for (uword o = 0; o < m; o++) {
                uword p = order[o];     // index of snp
                double alpha_exp_;      // new value of alpha_exp[p]
                double alpha_exp2_;     // new value of alpha_exp2[p]
                double alpha_var_;      // new value of alpha_var[p]
                double gamma_exp_;      // new value of gamma_exp[p]

                double temp = (dot(X.col(p), y_res) + xpx[p] * alpha_exp[p] * gamma_exp[p]) * tau_exp;

                alpha_var_ = 1.0 / (xpx[p] * tau_exp + 1.0 / S2);
                alpha_exp_  = alpha_var_ * temp;
                alpha_exp2_ = pow(alpha_exp_, 2.0) + alpha_var_;

                /* update Gamma */
                gamma_exp_ = 0.5 * alpha_var_ * temp * temp + 0.5 * log(alpha_var_);
                gamma_exp_ += gamma_constant;
                
                // to avoid overflow
                if (gamma_exp_ < 20.0) {
                    gamma_exp_ = exp(gamma_exp_);
                    gamma_exp_ = gamma_exp_ / (1.0 + gamma_exp_ - Pi[0]);
                } else {
                    gamma_exp_ = 1.0;   // exp(20) / (exp(20) + 1) ≈ 1 - 2.06e-9
                }

                /* update residuals */
                y_res += (X.col(p) * (alpha_exp[p] * gamma_exp[p] - alpha_exp_ * gamma_exp_));

                Check1 += pow(alpha_exp_ - alpha_exp[p], 2.0);
                Check2 += pow(alpha_exp_, 2.0);

                alpha_exp[p]  = alpha_exp_;
                alpha_exp2[p] = alpha_exp2_;
                alpha_var[p]  = alpha_var_;

                gamma_exp[p]  = gamma_exp_;
                gamme_exp2[p] = pow(gamma_exp_, 2.0) + gamma_exp_ * (1.0 - gamma_exp_);

                sumVarB     += (xpx[p] * gamma_exp_ * (alpha_exp2_ - gamma_exp_ * alpha_exp_ * alpha_exp_));
                sumGammaB2_ += alpha_exp2_ * gamma_exp_;
                sumGamma_   += gamma_exp_;
            }

            /* update of Sigma2 */
            // sumGammaB2[0] = sumGammaB2[1];
            // sumGamma[0] = sumGamma[1];
            // X[0].expSigma2[0] = (sumGammaB2[0] + vS2) / (H[0].v + sumGamma[0] - 2.0);
            // X[0].S2[0] = prop = (sumGammaB2[0] + vS2) / (H[0].v + sumGamma[0]);
            sigma2 = (sumGammaB2 + vS2) / (dfvg_ + sumGamma - 2.0);
            S2_ = (sumGammaB2 + vS2) / (dfvg_ + sumGamma);

            sumGammaB2 = sumGammaB2_;
            sumGamma   = sumGamma_;

            Check1 += pow(S2_ - S2, 2.0); 
            Check2 += pow(S2_, 2.0);
            S2 = S2_;
        break;
        case Model::BayesCpi:
        case Model::BayesL:
        case Model::BayesR:
        case Model::BayesRR:
        case Model::BSLMM:
        default:
            throw Rcpp::exception("Not implemented yet");
        }

		// Update of residual precision
		a1 = 0.5 * n;
		b1 = 0.5 * (accu(square(y_res)) + sumVarB);

		tau_exp_ = a1 / b1;
        tau_var_ = a1 / pow(b1, 2.0);

		Check1 += pow(tau_exp_ - tau_exp, 2.0);
		Check2 += pow(tau_exp_, 2.0);
		tau_exp = tau_exp_;

		// Lower bound
		// Rmonitor[ite - 1] = 1.0 / Tau0[0];

		if (iter == max_iteration || (Check1 / Check2) < threshold) {
			// LB[0] = - 0.5 * (double) Ny * log2pi - a1 * log(b1) + mylgamma(a1);
			// for (method = 0; method < Nm; method++)
			// {
			// 	LB[0] += LowerBound(Methodcode[method], Ny, P[method], X + method, a2[method], b2[method], H + method, Tau0[0]);
			// }
			break;
		}
        Rcpp::Rcout << "Iteration: " << iter << " finished." << endl;
	}

    List results;
        
    results["Vg"] = dot(alpha_var, vx);
    results["Ve"] = 1.0 / tau_exp;
    results["h2"] = dot(alpha_var, vx) / (dot(alpha_var, vx) + 1.0 / tau_exp);
    results["mu"] = mu_exp;
    results["Tau"] = tau_exp;
    results["vary"] = vary;
    // if(nc){
    //     beta = conv_to<vec>::from(mean(beta_store, 1));
    //     betasd = conv_to<vec>::from(stddev(beta_store, 0, 1));
    //     e -= (C_ * beta);
    //     results["beta"] = beta;
    // }

    results["alpha"] = alpha_exp;
    results["g"] = X * alpha_exp;
    results["e"] = y - X * alpha_exp - mu_exp;

    timer.step("end");
    NumericVector res(timer);
    uword tt  = floor((res[1] - res[0]) / 1e9);
    uword hor = floor(tt / 3600);
    uword min = floor(tt % 3600 / 60);
    uword sec = floor(tt % 3600 % 60);

    if (verbose)
        Rprintf("Finished within total run time: %02dh%02dm%02ds \n", hor, min, sec);
    
    return results;
}
