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
#include <algorithm>
#include <random>

#include "hibayes.h"
#include "solver.h"
#include "magic_enum.hpp"


enum class Model { BayesA, BayesB, BayesBpi, BayesC, BayesCpi, BayesL, BayesR, BayesRR, BSLMM };

vector<arma::span> createBlocks(const uword m, const uword block_size=64) {
    vector<span> blocks;
    if (block_size == 1) {
        for (uword i = 0; i < m; i ++) {
            blocks.push_back(span(i));
        }
    } else {
        for (uword i = 0; i < m; i += block_size) {
            blocks.push_back(span(i, min(i + block_size, m) - 1));
        }
    }
    return blocks;
}


void orth_augment(arma::mat &X, arma::vec &d, const arma::uword block_size) {
    uword n = X.n_rows;
    uword m = X.n_cols;
    X.resize(n + min(block_size, m), m);

    uword n_block = (m % block_size) == 0 ? m / block_size : m / block_size + 1;
    d.set_size(n_block);
    
    for (uword i = 0, b = 0; i < m; i += block_size) {
        uword p = min(block_size, m - i);
        auto X0 = X.submat(0, i, n - 1, i + p - 1);
        mat Xa = X0.t() * X0;
        d[b] = max(eig_sym(Xa));
        Xa.diag() -= (d[b] + 0.0001);
        X.submat(n, i, n + p - 1, i + p - 1) = chol(-Xa);

        b++;
    }
}

// // [[Rcpp::export]]
// Rcpp::List orth_augment_test(){
//     mat A = randu<mat>(3, 5);
//     vec d;
//     orth_augment(A, d, 2);
    
//     List results;
        
//     results["A"] = A;
//     results["d"] = d;

//     return results;
// }

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
    const Nullable<double> dfvg = R_NilValue,   // Nu of Inv-Chi2 distribution of additive gentic variance
    const Nullable<double> s2vg = R_NilValue,   // S2 of Inv-Chi2 distribution of additive gentic variance
    const Nullable<double> dfve = R_NilValue,   // Nu of Inv-Chi2 distribution of residual variance
    const Nullable<double> s2ve = R_NilValue,   // S2 of Inv-Chi2 distribution of residual variance
    arma::uword block_size = 64,
    const int threads = 0,
    const int max_iteration = 1000,
    const int seed = 42,
    const double threshold = 1e-5,
    const bool random_init = false,
    const bool verbose = true
) {
    omp_setup(threads);
    std::default_random_engine rng(seed);

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

    Rcpp::Rcout << "Initialize." << std::endl;
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
    // Orthogonal data augmentation
    vec xd;
    double yatya = 0.0;
    if (block_size > 1) {
        block_size = min(block_size, m);
        orth_augment(X, xd, block_size);
        y.resize(n + min(block_size, m));
        y.subvec(n, y.n_elem - 1) = zeros(min(block_size, m));
    }
    
    // Cache feature of markers
    vec xtx = zeros(m);     // Vector of (X_p * X_p)
    vec vx = zeros(m);      // Vector of variance of X_p
    // #pragma omp parallel for
    for(uword p = 0; p < m; p++) {
        auto Xp = X.col(p);
        xtx[p] = dot(Xp, Xp);
        vx[p] = var(Xp);
    }
    double sumvx = sum(vx);

    // parameters of the additive gentic effect distribution
    double dfvg_ = dfvg.isNotNull()? as<double>(dfvg) : 4;                           // Nu of Inv-Chi2 distribution
    // double vg_   = vg.isNotNull()? as<double>(vg) : ((dfvg_ - 2) / dfvg_) * vary * h2;  
    double s2vg_ = s2vg.isNotNull()? as<double>(s2vg) : h2 * (dfvg_ - 2) / (dfvg_ * (1 - Pi[0]) * sumvx);    // S2 of Inv-Chi2 distribution

    // parameters of the marker effect distribution
    // double vp_   = vg_   / ((1 - Pi[0]) * sumvx);
    // double s2vp_ = s2vg_ / ((1 - Pi[0]) * sumvx);

    vec beta_exp;      // E[β]
    vec beta_var;      // V[β]
    vec beta_exp2;     // E[β^2]

    vec sigma2_g;       // E[σ^2_g]: σ^2_g parameter of marker effects distribution β ~ N(0, σ^2)
    vec S2_g;           // S^2 parameter of σ^2_g ~ X^{-2}(v, S^2)
    
    vec gamma_exp;      // E[γ]: probability of marker have effect
    vec gamma_exp2;     // E[γ^2]

    // Initialize - model
    switch (model_index) {
    case Model::BayesA:
        beta_exp  = zeros(m);  // E[β] 
        beta_var  = zeros(m);  // V[β]
        beta_exp2 = zeros(m);  // E[β^2] = E[β]^2 + V[β]

        sigma2_g = vec(m, fill::value(1.0 / m));    // E[σ^2_g]
        S2_g = vec(m, fill::value(m));              // S'^2_p = (vS2 + E[β^2]) / v'
        break;
    case Model::BayesB:
    case Model::BayesBpi:
        break;
    case Model::BayesC:
    case Model::BayesCpi:
        beta_exp  = zeros(m);  // E[β] 
        beta_var  = zeros(m);  // V[β] 
        beta_exp2 = zeros(m);  // E[β^2] 

        gamma_exp  = vec(m, fill::value(0.5));      // E[γ]
        gamma_exp2 = vec(m, fill::value(0.5));      // E[γ^2]: gamma_exp ^ 2 + gamma_exp * (1 - gamma_exp)

        sigma2_g = vec(1, fill::value(1.0 / m));    // E[σ^2_g]
        S2_g = vec(1, fill::value(m));
        break;
    case Model::BayesL:
    case Model::BayesR:
    case Model::BayesRR:
    case Model::BSLMM:
        throw Rcpp::exception("Not implemented yet");
    }
    
    // Residual ===========================
    vec y_res = y - mu_exp;  // residual

    // parameters of the residual distribution
    double dfve_ = dfve.isNotNull()? as<double>(dfve) : -2;        // Nu of Inv-Chi2 distribution     
    double s2ve_ = s2ve.isNotNull()? as<double>(s2ve) : 0;         // S2 of Inv-Chi2 distribution
    // double ve_ = ve.isNotNull()? as<double>(ve) : vary * (1 - h2); //TODO: / (nr + 1); 
    
    double sigma2_e = 1.0 / (vary * (1 - h2));
    double sigma2_e_;


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

        Rcpp::Rcout << "Inference started: " << std::endl;
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
        double sumGammaB2 = 0.0;    // Σ(γ * E[β^2]) 
        double sumGamma   = 0.0;    // Σ(γ) : probability of marker have effect
	    double gamma_constant;

        double sumGammaB2_ = 0.0;
        double sumGamma_   = 0.0;
        double S2_g_ = 0.0;
        
        // double yty_a = 0.0;
        vec y_a = zeros(block_size);

		// For update of residual variance
        uvec order = linspace<uvec>(0, m - 1, m);     // order of snp {1 ... m}, used to update SNP Effect in random order
        auto blocks = createBlocks(m, block_size);
        
		// Update of beta
        if (model_index == Model::BayesA) {
            // s2vp_ = pow((dfvg_ - 2) / dfvg_, 2.0) * vary * h2 / ((1 - Pi[0]) * sumvx); 

            vS2 = dfvg_ * s2vg_;

            order = shuffle(order);     // update SNP Effect in random order
            for (uword o = 0; o < m; o++) {
                uword p = order[o];    // index of snp
                double beta_exp_;      // new value of beta_exp[p]
                // double beta_exp2_;     // new value of beta_exp2[p]
                double beta_var_;      // new value of beta_var[p]
                double S2_g_;          // new value of S2_g[p]

                // beta_var_ = 1.0 / (xtx[p] * tau_exp + 1.0 / S2_g[p]);
                beta_var_  = (sigma2_e * S2_g[p]) / (xtx[p] * S2_g[p] + sigma2_e);
                beta_exp_  = beta_var_ * (dot(X.col(p), y_res) + xtx[p] * beta_exp[p]) / sigma2_e;
                // beta_exp2_ = pow(beta_exp_, 2.0) + beta_var_;

                /* update residuals */
                y_res   += (X.col(p) * (beta_exp[p] - beta_exp_));
                // sumVarB += (xtx[p] * beta_var[p]);     // sum(V[β_p])

                Check1 += pow(beta_exp_ - beta_exp[p], 2.0);
                Check2 += pow(beta_exp_, 2.0);

                beta_exp[p]  = beta_exp_;
                // beta_exp2[p] = beta_exp2_;
                beta_var[p]  = beta_var_;

                /* update of Sigma2 */
                sigma2_g[p] = (beta_exp2[p] + vS2) / (dfvg_ - 2.0);
                S2_g_       = (beta_exp2[p] + vS2) / (dfvg_ + 1.0);

                Check1 += pow(S2_g_ - S2_g[p], 2.0); 
                Check2 += pow(S2_g_, 2.0);

                S2_g[p] = S2_g_;
            }
        } else if (model_index == Model::BayesB || model_index == Model::BayesBpi) {
            throw Rcpp::exception("Not implemented yet");
        } else if (model_index == Model::BayesC || model_index == Model::BayesCpi) {
            // s2vp_ = pow((dfvg_ - 2) / dfvg_, 2.0) * vary * h2 / ((1 - Pi[0]) * sumvx); 
            vS2 = dfvg_ * s2vg_;

            gamma_constant = 0.5 * R::digamma(0.5 * (dfvg_ + sumGamma)) - 0.5 * log(0.5 * (sumGammaB2 + vS2));

            // update SNP Effect in random order
            // order = shuffle(order);     
            std::shuffle(blocks.begin(), blocks.end(), rng);
            
            for (span b: blocks) {
                uword n_elem = b.b - b.a + 1;
                vec beta_exp_(n_elem);       // new value of beta_exp[p]
                vec beta_exp2_(n_elem);      // new value of beta_exp2[p]
                vec beta_var_(n_elem);       // H_p
                vec gamma_exp_(n_elem);      // new value of gamma_exp[p]
                
                // yty_a = 0.0;
                if (block_size > 1) {
                    // y_res.subvec(n, n + block_size - 1) = X.submat(n, b.start, n + block_size - 1, b.end - 1) * (beta_exp.subvec(b) % gamma_exp.subvec(b));
                    y_res.subvec(n, n + block_size - 1).fill(0.0);
                }

                vec temp = (X.cols(b).t() * y_res + xtx.subvec(b) % beta_exp.subvec(b) % gamma_exp.subvec(b)) / sigma2_e;
                beta_var_ = (sigma2_e * S2_g[0]) / (xtx.subvec(b) * S2_g[0] + sigma2_e);
                beta_exp_ = beta_var_ % temp;

                /* update Gamma */
                gamma_exp_ = 0.5 * beta_var_ % temp % temp + 0.5 * log(beta_var_);
                gamma_exp_ += gamma_constant;
                gamma_exp_.elem(find(gamma_exp_ > 20.0)).fill(20.0);    // to avoid overflow
                gamma_exp_ = Pi[1] * exp(gamma_exp_);
                gamma_exp_ /= (gamma_exp_ + Pi[0]);
                
                Check1 += accu(square(beta_exp_ - beta_exp.subvec(b)));
                Check2 += accu(square(beta_exp_));

                sumGammaB2_ += accu(beta_exp2_ % gamma_exp_);
                sumGamma_   += accu(gamma_exp_);

                // Update Residual

                // vec y_res_a = X.submat(n, b.start, n + block_size - 1, b.end - 1) * (beta_exp.subvec(b) % gamma_exp.subvec(b) - beta_exp_ % gamma_exp_);
                // yty_a += dot(y_res_a, y_res_a);
                y_res += X.cols(b) * (beta_exp.subvec(b) % gamma_exp.subvec(b) - beta_exp_ % gamma_exp_);
                // if (block_size > 1) {
                //     y_a += y_res.subvec(n, y_res.n_rows - 1);
                // }
                beta_exp.subvec(b) = beta_exp_;
                beta_var.subvec(b) = beta_var_;
                gamma_exp.subvec(b) = gamma_exp_;
                
                // y_res.subvec(n, n + block_size - 1) = zeros(block_size);
            }

            /* update of Sigma2 */
            //   v' = v + sumGamma
            //  S2' = (vS2 + sum(β^2)) / (v + sumGamma)
            // mean = vS2 / (v - 2)

            sigma2_g[0] = (sumGammaB2 + vS2) / (dfvg_ + sumGamma - 2.0);
            S2_g_       = (sumGammaB2 + vS2) / (dfvg_ + sumGamma);       // S~ 

            sumGammaB2 = sumGammaB2_;
            sumGamma   = sumGamma_;

            Check1 += pow(S2_g_ - S2_g[0], 2.0); 
            Check2 += pow(S2_g_, 2.0);
            S2_g[0] = S2_g_;

            // Update π
            if (model_index == Model::BayesCpi) {
                double Pi1_ =  (sumGamma + 1) / (m + 2); // E[π], π ~ Beta(m - k + 1, k + 1)
                // double Pi1_ =  sumGamma / m;
                Check1 += pow(Pi1_ - Pi[1], 2.0);
                Check2 += pow(Pi1_, 2.0);

                Pi[0] = 1 - Pi1_;
                Pi[1] = Pi1_;
            }
        } else if (model_index == Model::BayesL) {
            throw Rcpp::exception("Not implemented yet");
        } else if (model_index == Model::BayesR) {
            throw Rcpp::exception("Not implemented yet");
        } else if (model_index == Model::BayesRR) {
            throw Rcpp::exception("Not implemented yet");
        } else if (model_index == Model::BSLMM) {
        } else {
            throw Rcpp::exception("Not implemented yet");
        }

        if (block_size > 1) {
            y_res.subvec(n, y_res.n_elem - 1) = y_a;
            Rcpp::Rcout << "Debug: " << y_a.t();
        }

        // Update mu
        double mu_exp_;
        mu_exp_ = mu_exp + mean(y_res);

        Check1 += pow(mu_exp_ - mu_exp, 2.0);
        Check2 += pow(mu_exp_, 2.0);

        y_res += (mu_exp - mu_exp_);
        mu_exp = mu_exp_;

		// Update of residual precision
        sigma2_e_ = dot(y_res, y_res) / (y_res.n_elem - 2.0);

		Check1 += pow(sigma2_e_ - sigma2_e, 2.0);
		Check2 += pow(sigma2_e_, 2.0);
		sigma2_e = sigma2_e_;

		// Lower bound
		if (iter == max_iteration || (Check1 / Check2) < threshold) {
			// LB[0] = - 0.5 * (double) Ny * log2pi - a1 * log(b1) + mylgamma(a1);
			// for (method = 0; method < Nm; method++)
			// {
			// 	LB[0] += LowerBound(Methodcode[method], Ny, P[method], X + method, a2[method], b2[method], H + method, Tau0[0]);
			// }
			break;
		}
        Rcpp::Rcout << "Iteration: " << iter << " finished." << endl;
        // Rcpp::Rcout << "π: " << Pi[1] << endl;
	}

    Rcpp::Rcout << "Iteration finished." << endl;
    List results;
        
    results["Vg"] = dot(beta_var, vx);
    results["Ve"] = sigma2_e;
    results["h2"] = dot(beta_var, vx) / (dot(beta_var, vx) + 1.0 / sigma2_e);
    results["mu"] = mu_exp;
    results["vary"] = vary;
    results["Pi"] = Pi;
    // if(nc){
    //     beta = conv_to<vec>::from(mean(beta_store, 1));
    //     betasd = conv_to<vec>::from(stddev(beta_store, 0, 1));
    //     e -= (C_ * beta);
    //     results["beta"] = beta;
    // }

    results["beta"] = beta_exp % gamma_exp;
    results["gamma"] = gamma_exp;
    results["g"] = X.rows(0, n - 1) * beta_exp;
    results["e"] = y.rows(0, n - 1) - X.rows(0, n - 1) * beta_exp - mu_exp;

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
