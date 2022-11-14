#include "hibayes.h"
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <progress.hpp>
#include "progress_bar.hpp"

// [[Rcpp::depends(bigmemory, BH)]]
// [[Rcpp::depends(RcppProgress)]]

class MinimalProgressBar: public ProgressBar{
    public:
    MinimalProgressBar(const char *str = "Calculating in process")  {
        _finalized = false;
        _str = str;
    }
    ~MinimalProgressBar() {}
    void display() {}
    void update(float progress) {
        if (_finalized) return;
        int pi = (int)(progress * point_length);
        if(point[pi]){
            point[pi] = false;
            REprintf("\r");
            REprintf(_str);
            REprintf("...finished %u%", (int)(progress * 100));
        }
    }
    void end_display() {
    if (_finalized) return;
        REprintf("\r");
        REprintf(_str);
        REprintf("...[finished 100%]");
        REprintf("\n");
        _finalized = true;
    }
    private:
    bool _finalized;
    const char *_str;
    int point_length = 100;
    LogicalVector point = rep(true, point_length);
};

template <typename T>
SEXP BigStat(XPtr<BigMatrix> pMat, const int threads = 0){

    omp_setup(threads);

	MatrixAccessor<T> bigm = MatrixAccessor<T>(*pMat);

	int ind = pMat->nrow();
	int j, k, m = pMat->ncol();
	double scale_mean, p1 = 0.0;
	NumericVector mean(m);
	NumericVector sd(m);
	NumericVector sum(m);

	#pragma omp parallel for private(p1, k)
	for (j = 0; j < m; j++){
		p1 = 0.0;
		for(k = 0; k < ind; k++){
			p1 += bigm[j][k];
		}
		sum[j] = p1;
		mean[j] = p1 / ind;
	}

	#pragma omp parallel for private(p1, k, scale_mean)
	for (j = 0; j < m; j++){
		p1 = 0.0;
		for(k = 0; k < ind; k++){
			scale_mean = (bigm[j][k] - mean[j]);
			p1 += scale_mean * scale_mean;
		}
		sd[j] = sqrt(p1);
	}
	return List::create(Named("mean") = mean, Named("sum") = sum, Named("xx") = sd);
}

// [[Rcpp::export]]
SEXP BigStat(SEXP pBigMat, const int threads = 0){
	
	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return BigStat<char>(xpMat, threads);
	case 2:
		return BigStat<short>(xpMat, threads);
	case 4:
		return BigStat<int>(xpMat, threads);
	case 6:
		return BigStat<float>(xpMat, threads);
	case 8:
		return BigStat<double>(xpMat, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP tXXmat_Geno(XPtr<BigMatrix> pMat, const Nullable<double> chisq = R_NilValue, const int threads=0, const bool verbose=true){

	omp_setup(threads);

	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int m = pMat->ncol();
	int ind = pMat->nrow();
	int i, j, k;
	bool sparse = false;
	double m1, m2, sum1, sum2, p1, p2, p12, r, chisq_;
	MinimalProgressBar pb;
	List Stat = BigStat(pMat, threads);
	NumericVector mean_all = Stat[0];
	NumericVector sum_all = Stat[1];
	NumericVector xx_all = Stat[2];
	if(chisq.isNotNull()){
		chisq_ = as<double>(chisq);
		if(chisq_ > 0)	sparse = true;
	}
	
	Progress p(m, verbose, pb);
	
	if(sparse){
		if(verbose)	Rcerr << "Genome-Wide sparse matrix" << endl;
		arma::sp_mat ldmat(m, m);
		arma::vec r2vec(m); r2vec.zeros();

		#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, k, p12, p2, m2, sum2, r)
		for (j = 0; j < m; j++){
			if ( ! Progress::check_abort() ) {
				p.increment();
				p1 = xx_all[j];
				m1 = mean_all[j];
				sum1 = sum_all[j];

				for(i = j; i < m; i++){
					p12 = 0;
					p2 = xx_all[i];
					m2 = mean_all[i];
					sum2 = sum_all[i];
					for(k = 0; k < ind; k++){
						p12 += (genomat[i][k]) * (genomat[j][k]);
					}
					p12 -= sum1 * m2 + sum2 * m1 - ind * m1 * m2;
					r = p12 / (p1 * p2);
					if(r * r * ind <= chisq_){
						// nothing
					}else{

						#pragma omp critical
						ldmat(j, i) = ldmat(i, j) = p12 / ind;
					}
				}
			}
		}
		return wrap(ldmat);
	}else{
		if(verbose)	Rcerr << "Genome-Wide dense matrix" << endl;
		arma::mat ldmat(m, m);
		#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
		for (j = 0; j < m; j++){
			if ( ! Progress::check_abort() ) {
				p.increment();
				p1 = xx_all[j];
				m1 = mean_all[j];
				sum1 = sum_all[j];
				ldmat(j, j) = p1 * p1 / ind;
				for(i = j + 1; i < m; i++){
					p12 = 0;
					p2 = xx_all[i];
					m2 = mean_all[i];
					sum2 = sum_all[i];
					for(k = 0; k < ind; k++){
						p12 += (genomat[i][k]) * (genomat[j][k]);
					}
					p12 -= sum1 * m2 + sum2 * m1 - ind * m1 * m2;
					// r = p12 / (p1 * p2);
					ldmat(i, j) = ldmat(j, i) = p12 / ind;
				}
			}
		}
		return wrap(ldmat);
	}
}

// [[Rcpp::export]]
SEXP tXXmat_Geno(SEXP pBigMat, const Nullable<double> chisq = R_NilValue, const int threads=0, const bool verbose=true){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return tXXmat_Geno<char>(xpMat, chisq, threads, verbose);
	case 2:
		return tXXmat_Geno<short>(xpMat, chisq, threads, verbose);
	case 4:
		return tXXmat_Geno<int>(xpMat, chisq, threads, verbose);
	case 6:
		return tXXmat_Geno<float>(xpMat, chisq, threads, verbose);
	case 8:
		return tXXmat_Geno<double>(xpMat, chisq, threads, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

// template <typename T>
// SEXP tXXmat_Geno_1(XPtr<BigMatrix> pMat, const Nullable<double> chisq = R_NilValue, const int threads=0, const bool verbose=true){

// 	if (threads == 0) {
// 		omp_set_num_threads(omp_get_num_procs());
// 	}else if(threads > 0) {
// 		omp_set_num_threads(threads);
// 	}

// 	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

// 	int m = pMat->ncol();
// 	int ind = pMat->nrow();
// 	int i, j, k;
// 	bool sparse = false;
// 	double m1, m2, sum1, sum2, p1, p2, p12, r, chisq_;
// 	MinimalProgressBar pb;
// 	List Stat = BigStat(pMat, threads);
// 	NumericVector mean_all = Stat[0];
// 	NumericVector sum_all = Stat[1];
// 	NumericVector xx_all = Stat[2];
// 	if(chisq.isNotNull()){
// 		chisq_ = as<double>(chisq);
// 		if(chisq_ > 0)	sparse = true;
// 	}
	
// 	Progress p(m, verbose, pb);
	
// 	if(sparse){
// 		arma::sp_mat ldmat(m, m);
// 		#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
// 		for (j = 0; j < m; j++){
// 			if ( ! Progress::check_abort() ) {
// 				p.increment();
// 				p1 = xx_all[j];
// 				m1 = mean_all[j];
// 				sum1 = sum_all[j];
// 				ldmat(j, j) = p1 * p1 / ind;
// 				for(i = j + 1; i < m; i++){
// 					p12 = 0;
// 					p2 = xx_all[i];
// 					m2 = mean_all[i];
// 					sum2 = sum_all[i];
// 					for(k = 0; k < ind; k++){
// 						p12 += (genomat[i][k]) * (genomat[j][k]);
// 					}
// 					p12 -= sum1 * m2 + sum2 * m1 - ind * m1 * m2;
// 					r = p12 / (p1 * p2);
// 					if(r * r * ind <= chisq_){
// 						// nothing
// 					}else{
// 						ldmat(i, j) = ldmat(j, i) = p12 / ind;
// 					}
// 				}
// 			}
// 		}
// 		return wrap(ldmat);
// 	}else{
// 		arma::mat ldmat(m, m);
// 		#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
// 		for (j = 0; j < m; j++){
// 			if ( ! Progress::check_abort() ) {
// 				p.increment();
// 				p1 = xx_all[j];
// 				m1 = mean_all[j];
// 				sum1 = sum_all[j];
// 				ldmat(j, j) = p1 * p1 / ind;
// 				for(i = j + 1; i < m; i++){
// 					p12 = 0;
// 					p2 = xx_all[i];
// 					m2 = mean_all[i];
// 					sum2 = sum_all[i];
// 					for(k = 0; k < ind; k++){
// 						p12 += (genomat[i][k]) * (genomat[j][k]);
// 					}
// 					p12 -= sum1 * m2 + sum2 * m1 - ind * m1 * m2;
// 					// r = p12 / (p1 * p2);
// 					ldmat(i, j) = ldmat(j, i) = p12 / ind;
// 				}
// 			}
// 		}
// 		return wrap(ldmat);
// 	}
// }

// // [[Rcpp::export]]
// SEXP tXXmat_Geno_1(SEXP pBigMat, const Nullable<double> chisq = R_NilValue, const int threads=0, const bool verbose=true){

// 	XPtr<BigMatrix> xpMat(pBigMat);

// 	switch(xpMat->matrix_type()){
// 	case 1:
// 		return tXXmat_Geno_1<char>(xpMat, chisq, threads, verbose);
// 	case 2:
// 		return tXXmat_Geno_1<short>(xpMat, chisq, threads, verbose);
// 	case 4:
// 		return tXXmat_Geno_1<int>(xpMat, chisq, threads, verbose);
// 	case 6:
// 		return BigStat<float>(xpMat, chisq, threads, verbose);
// 	case 8:
// 		return tXXmat_Geno_1<double>(xpMat, chisq, threads, verbose);
// 	default:
// 		throw Rcpp::exception("unknown type detected for big.matrix object!");
// 	}
// }

template <typename T>
SEXP tXXmat_Geno_gwas(XPtr<BigMatrix> pMat, SEXP gwasgeno, const LogicalVector refindx, const NumericVector gwasindx, const Nullable<double> chisq = R_NilValue, const int threads=0, const bool verbose=true){

	omp_setup(threads);

	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);
	XPtr<BigMatrix> gwasMat(gwasgeno);
	MatrixAccessor<T> gwasgenomat = MatrixAccessor<T>(*gwasMat);

	int m = pMat->ncol();
	int ind = pMat->nrow();
	int mgwas = gwasMat->ncol();
	int indgwas = gwasMat->nrow();
	int i, j, k;
	bool sparse = false;
	double m1, m2, sum1, sum2, p1, p2, p12, r, chisq_;
	List Stat = BigStat(pMat, threads);
	NumericVector mean_all = Stat[0];
	NumericVector sum_all = Stat[1];
	NumericVector xx_all = Stat[2];
	List Stat_gwas = BigStat(gwasMat, threads);
	NumericVector mean_all_gwas = Stat_gwas[0];
	NumericVector sum_all_gwas = Stat_gwas[1];
	NumericVector xx_all_gwas = Stat_gwas[2];
	if(chisq.isNotNull()){
		chisq_ = as<double>(chisq);
		if(chisq_ > 0)	sparse = true;
	}

	if(sparse){
		if(verbose)	Rcerr << "Genome-Wide sparse matrix" << endl;
		arma::sp_mat ldmat(m, m);
		for(int xx = 0; xx < 2; xx++){
			if(xx == 0){
				MinimalProgressBar pb;
				Progress p(m, verbose, pb);

				#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
				for (j = 0; j < m; j++){
					if ( ! Progress::check_abort() ) {
						p.increment();
						p1 = xx_all[j];
						m1 = mean_all[j];
						sum1 = sum_all[j];
						for(i = j; i < m; i++){
							if(refindx[j] && refindx[i]){
								// nothing
							}else{
								p12 = 0;
								p2 = xx_all[i];
								m2 = mean_all[i];
								sum2 = sum_all[i];
								for(k = 0; k < ind; k++){
									p12 += (genomat[i][k]) * (genomat[j][k]);
								}
								p12 -= sum1 * m2 + sum2 * m1 - ind * m1 * m2;
								r = p12 / (p1 * p2);
								if(r * r * ind <= chisq_){
									// nothing
								}else{
									
									#pragma omp critical
									ldmat(j, i) = ldmat(i, j) = p12 / ind;
								}
							}
						}
					}
				}
			}else{
				if(verbose)	Rcerr << "Update LD for SNPs in GWAS panel" << endl;
				MinimalProgressBar pb;
				Progress p(mgwas, verbose, pb);

				#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
				for (j = 0; j < mgwas; j++){
					if ( ! Progress::check_abort() ) {
						p.increment();
						p1 = xx_all_gwas[j];
						m1 = mean_all_gwas[j];
						sum1 = sum_all_gwas[j];
						for(i = j; i < mgwas; i++){
							p12 = 0;
							p2 = xx_all_gwas[i];
							m2 = mean_all_gwas[i];
							sum2 = sum_all_gwas[i];
							for(k = 0; k < indgwas; k++){
								p12 += (gwasgenomat[i][k]) * (gwasgenomat[j][k]);
							}
							p12 -= sum1 * m2 + sum2 * m1 - indgwas * m1 * m2;
							r = p12 / (p1 * p2);
							if(r * r * ind <= chisq_){
								// nothing
							}else{
								
								#pragma omp critical
								ldmat(gwasindx[j], gwasindx[i]) = ldmat(gwasindx[i], gwasindx[j]) = p12 / indgwas;
							}
						}
					}
				}
			}
		}
		return wrap(ldmat);
	}else{
		if(verbose)	Rcerr << "Genome-Wide dense matrix" << endl;
		arma::mat ldmat(m, m);
		for(int xx = 0; xx < 2; xx++){
			if(xx == 0){
				MinimalProgressBar pb;
				Progress p(m, verbose, pb);

				#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
				for (j = 0; j < m; j++){
					if ( ! Progress::check_abort() ) {
						p.increment();
						p1 = xx_all[j];
						m1 = mean_all[j];
						sum1 = sum_all[j];
						ldmat(j, j) = p1 * p1 / ind;
						for(i = j + 1; i < m; i++){
							if(refindx[j] && refindx[i]){
								// nothing
							}else{
								p12 = 0;
								p2 = xx_all[i];
								m2 = mean_all[i];
								sum2 = sum_all[i];
								for(k = 0; k < ind; k++){
									p12 += (genomat[i][k]) * (genomat[j][k]);
								}
								p12 -= sum1 * m2 + sum2 * m1 - ind * m1 * m2;
								// r = p12 / (p1 * p2);
								ldmat(i, j) = ldmat(j, i) = p12 / ind;
							}
						}
					}
				}
			}else{
				if(verbose)	Rcerr << "Update LD for SNPs in GWAS panel" << endl;
				MinimalProgressBar pb;
				Progress p(mgwas, verbose, pb);

				#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
				for (j = 0; j < mgwas; j++){
					if ( ! Progress::check_abort() ) {
						p.increment();
						p1 = xx_all_gwas[j];
						m1 = mean_all_gwas[j];
						sum1 = sum_all_gwas[j];
						for(i = j; i < mgwas; i++){
							p12 = 0;
							p2 = xx_all_gwas[i];
							m2 = mean_all_gwas[i];
							sum2 = sum_all_gwas[i];
							for(k = 0; k < indgwas; k++){
								p12 += (gwasgenomat[i][k]) * (gwasgenomat[j][k]);
							}
							p12 -= sum1 * m2 + sum2 * m1 - indgwas * m1 * m2;
							// r = p12 / (p1 * p2);
							ldmat(gwasindx[i], gwasindx[j]) = ldmat(gwasindx[j], gwasindx[i]) = p12 / indgwas;
						}
					}
				}
			}
		}
		return wrap(ldmat);
	}
}

// [[Rcpp::export]]
SEXP tXXmat_Geno_gwas(SEXP pBigMat, SEXP gwasgeno, const LogicalVector refindx, const NumericVector gwasindx, const Nullable<double> chisq = R_NilValue, const int threads=0, const bool verbose=true){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return tXXmat_Geno_gwas<char>(xpMat, gwasgeno, refindx, gwasindx, chisq, threads, verbose);
	case 2:
		return tXXmat_Geno_gwas<short>(xpMat, gwasgeno, refindx, gwasindx, chisq, threads, verbose);
	case 4:
		return tXXmat_Geno_gwas<int>(xpMat, gwasgeno, refindx, gwasindx, chisq, threads, verbose);
	case 6:
		return tXXmat_Geno_gwas<float>(xpMat, gwasgeno, refindx, gwasindx, chisq, threads, verbose);
	case 8:
		return tXXmat_Geno_gwas<double>(xpMat, gwasgeno, refindx, gwasindx, chisq, threads, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP tXXmat_Chr(XPtr<BigMatrix> pMat, const NumericVector chr, const Nullable<double> chisq = R_NilValue, const int threads=0, const bool verbose=true){

	omp_setup(threads);

	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);

	int m = pMat->ncol();
	int ind = pMat->nrow();
	int i, j, k;
	bool sparse = false;
	double m1, m2, sum1, sum2, p1, p2, p12, r, chisq_;
	List Stat = BigStat(pMat, threads);
	NumericVector mean_all = Stat[0];
	NumericVector sum_all = Stat[1];
	NumericVector xx_all = Stat[2];
	if(chisq.isNotNull()){
		chisq_ = as<double>(chisq);
		sparse = true;
	}
	arma::vec vecchr = Rcpp::as<arma::vec>(chr);
	arma::vec unichr = unique(vecchr);
	
	if(sparse){
		if(verbose)	Rcerr << "Chromosome-Wide sparse matrix" << endl;
		arma::sp_mat ldmat(m, m);
		for(int cc = 0; cc < unichr.n_elem; cc++){
			uvec chrindx = find(vecchr == unichr[cc]);
			if(verbose)	Rcpp::Rcout << "Loop on chromosome No." << cc + 1 << " with total number of SNPs " << chrindx.n_elem << std::endl;
			MinimalProgressBar pb;
			Progress p(chrindx.n_elem, verbose, pb);

			#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
			for (j = 0; j < chrindx.n_elem; j++){
				if ( ! Progress::check_abort() ) {
					p.increment();
					p1 = xx_all[chrindx[j]];
					m1 = mean_all[chrindx[j]];
					sum1 = sum_all[chrindx[j]];
					for(i = j; i < chrindx.n_elem; i++){
						p12 = 0;
						p2 = xx_all[chrindx[i]];
						m2 = mean_all[chrindx[i]];
						sum2 = sum_all[chrindx[i]];
						for(k = 0; k < ind; k++){
							p12 += (genomat[chrindx[i]][k]) * (genomat[chrindx[j]][k]);
						}
						p12 -= sum1 * m2 + sum2 * m1 - ind * m1 * m2;
						r = p12 / (p1 * p2);
						if(r * r * ind <= chisq_){
							// nothing
						}else{
							
							#pragma omp critical
							ldmat(chrindx[i], chrindx[j]) = ldmat(chrindx[j], chrindx[i]) = p12 / ind;
						}
					}
				}
			}
		}
		return wrap(ldmat);
	}else{
		if(verbose)	Rcerr << "Chromosome-Wide dense matrix" << endl;
		arma::sp_mat ldmat(m, m);
		for(int cc = 0; cc < unichr.n_elem; cc++){
			uvec chrindx = find(vecchr == unichr[cc]);
			if(verbose)	Rcpp::Rcout << "Loop on chromosome No." << cc + 1 << " with total number of SNPs " << chrindx.n_elem << std::endl;
			MinimalProgressBar pb;
			arma::mat ldmat_chr(chrindx.n_elem, chrindx.n_elem);
			Progress p(chrindx.n_elem, verbose, pb);

			#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
			for (j = 0; j < chrindx.n_elem; j++){
				if ( ! Progress::check_abort() ) {
					p.increment();
					p1 = xx_all[chrindx[j]];
					m1 = mean_all[chrindx[j]];
					sum1 = sum_all[chrindx[j]];
					ldmat_chr(j, j) = p1 * p1 / ind;
					for(i = j + 1; i < chrindx.n_elem; i++){
						p12 = 0;
						p2 = xx_all[chrindx[i]];
						m2 = mean_all[chrindx[i]];
						sum2 = sum_all[chrindx[i]];
						for(k = 0; k < ind; k++){
							p12 += (genomat[chrindx[i]][k]) * (genomat[chrindx[j]][k]);
						}
						p12 -= sum1 * m2 + sum2 * m1 - ind * m1 * m2;
						// r = p12 / (p1 * p2);
						ldmat_chr(i, j) = ldmat_chr(j, i) = p12 / ind;
					}
				}
			}
			for (j = 0; j < chrindx.n_elem; j++){
				for(i = j; i < chrindx.n_elem; i++){
					ldmat(chrindx[i], chrindx[j]) = ldmat(chrindx[j], chrindx[i]) = ldmat_chr(i, j);
				}
			}
		}
		return wrap(ldmat);
	}
}

// [[Rcpp::export]]
SEXP tXXmat_Chr(SEXP pBigMat, const NumericVector chr, const Nullable<double> chisq = R_NilValue, const int threads=0, const bool verbose=true){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return tXXmat_Chr<char>(xpMat, chr, chisq, threads, verbose);
	case 2:
		return tXXmat_Chr<short>(xpMat, chr, chisq, threads, verbose);
	case 4:
		return tXXmat_Chr<int>(xpMat, chr, chisq, threads, verbose);
	case 6:
		return tXXmat_Chr<float>(xpMat, chr, chisq, threads, verbose);
	case 8:
		return tXXmat_Chr<double>(xpMat, chr, chisq, threads, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}

template <typename T>
SEXP tXXmat_Chr_gwas(XPtr<BigMatrix> pMat, const NumericVector chr, SEXP gwasgeno, const NumericVector gwaschr, const LogicalVector refindx, const NumericVector gwasindx, const Nullable<double> chisq = R_NilValue, const int threads=0, const bool verbose=true){

	omp_setup(threads);

	MatrixAccessor<T> genomat = MatrixAccessor<T>(*pMat);
	XPtr<BigMatrix> gwasMat(gwasgeno);
	MatrixAccessor<T> gwasgenomat = MatrixAccessor<T>(*gwasMat);

	int m = pMat->ncol();
	int ind = pMat->nrow();
	// int mgwas = gwasMat->ncol();
	int indgwas = gwasMat->nrow();
	int i, j, k;
	bool sparse = false;
	double m1, m2, sum1, sum2, p1, p2, p12, r, chisq_;
	List Stat = BigStat(pMat, threads);
	NumericVector mean_all = Stat[0];
	NumericVector sum_all = Stat[1];
	NumericVector xx_all = Stat[2];
	List Stat_gwas = BigStat(gwasMat, threads);
	NumericVector mean_all_gwas = Stat_gwas[0];
	NumericVector sum_all_gwas = Stat_gwas[1];
	NumericVector xx_all_gwas = Stat_gwas[2];
	if(chisq.isNotNull()){
		chisq_ = as<double>(chisq);
		sparse = true;
	}
	arma::vec vecchr = Rcpp::as<arma::vec>(chr);
	arma::vec unichr = unique(vecchr);
	arma::vec vecgwaschr = Rcpp::as<arma::vec>(gwaschr);
	arma::vec unigwaschr = unique(vecgwaschr);
	
	if(sparse){
		if(verbose)	Rcerr << "Chromosome-Wide sparse matrix" << endl;
		arma::sp_mat ldmat(m, m);
		for(int cc = 0; cc < unichr.n_elem; cc++){
			uvec chrindx = find(vecchr == unichr[cc]);
			if(verbose)	Rcpp::Rcout << "Loop on chromosome No." << cc + 1 << " with total number of SNPs in reference panel" << chrindx.n_elem << std::endl;
			MinimalProgressBar pb;
			Progress p(chrindx.n_elem, verbose, pb);

			#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
			for (j = 0; j < chrindx.n_elem; j++){
				if ( ! Progress::check_abort() ) {
					p.increment();
					p1 = xx_all[chrindx[j]];
					m1 = mean_all[chrindx[j]];
					sum1 = sum_all[chrindx[j]];
					for(i = j; i < chrindx.n_elem; i++){
						if(refindx[j] && refindx[i]){
						// nothing
						}else{
							p12 = 0;
							p2 = xx_all[chrindx[i]];
							m2 = mean_all[chrindx[i]];
							sum2 = sum_all[chrindx[i]];
							for(k = 0; k < ind; k++){
								p12 += (genomat[chrindx[i]][k]) * (genomat[chrindx[j]][k]);
							}
							p12 -= sum1 * m2 + sum2 * m1 - ind * m1 * m2;
							// r = p12 / (p1 * p2);
							r = p12 / (p1 * p2);
							if(r * r * ind <= chisq_){
								// nothing
							}else{
								
								#pragma omp critical
								ldmat(chrindx[i], chrindx[j]) = ldmat(chrindx[j], chrindx[i]) = p12 / ind;
							}
							
						}
					}
				}
			}
			chrindx = find(vecgwaschr == unichr[cc]);
			if(chrindx.n_elem > 0){
				if(verbose)	Rcpp::Rcout << "Loop on chromosome No." << cc + 1 << " with total number of SNPs in GWAS panel" << chrindx.n_elem << std::endl;
				MinimalProgressBar pb;
				Progress p(chrindx.n_elem, verbose, pb);
				#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
				for (j = 0; j < chrindx.n_elem; j++){
					if ( ! Progress::check_abort() ) {
						p.increment();
						p1 = xx_all_gwas[chrindx[j]];
						m1 = mean_all_gwas[chrindx[j]];
						sum1 = sum_all_gwas[chrindx[j]];
						for(i = j; i < chrindx.n_elem; i++){
							p12 = 0;
							p2 = xx_all_gwas[chrindx[i]];
							m2 = mean_all_gwas[chrindx[i]];
							sum2 = sum_all_gwas[chrindx[i]];
							for(k = 0; k < indgwas; k++){
								p12 += (gwasgenomat[chrindx[i]][k]) * (gwasgenomat[chrindx[j]][k]);
							}
							p12 -= sum1 * m2 + sum2 * m1 - ind * m1 * m2;
							r = p12 / (p1 * p2);
							if(r * r * ind <= chisq_){
								// nothing
							}else{
								
								#pragma omp critical
								ldmat(gwasindx[chrindx[i]], gwasindx[chrindx[j]]) = ldmat(gwasindx[chrindx[j]], gwasindx[chrindx[i]]) = p12 / indgwas;
							}
						}
					}
				}
			}
		}
		return wrap(ldmat);
	}else{
		if(verbose)	Rcerr << "Chromosome-Wide dense matrix" << endl;
		arma::sp_mat ldmat(m, m);
		for(int cc = 0; cc < unichr.n_elem; cc++){
			uvec chrindx = find(vecchr == unichr[cc]);
			if(verbose)	Rcpp::Rcout << "Loop on chromosome No." << cc + 1 << " with total number of SNPs in reference panel" << chrindx.n_elem << std::endl;
			MinimalProgressBar pb;
			arma::mat ldmat_chr(chrindx.n_elem, chrindx.n_elem);
			Progress p(chrindx.n_elem, verbose, pb);

			#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
			for (j = 0; j < chrindx.n_elem; j++){
				if ( ! Progress::check_abort() ) {
					p.increment();
					p1 = xx_all[chrindx[j]];
					m1 = mean_all[chrindx[j]];
					sum1 = sum_all[chrindx[j]];
					ldmat_chr(j, j) = p1 * p1 / ind;
					for(i = j + 1; i < chrindx.n_elem; i++){
						if(refindx[j] && refindx[i]){
						// nothing
						}else{
							p12 = 0;
							p2 = xx_all[chrindx[i]];
							m2 = mean_all[chrindx[i]];
							sum2 = sum_all[chrindx[i]];
							for(k = 0; k < ind; k++){
								p12 += (genomat[chrindx[i]][k]) * (genomat[chrindx[j]][k]);
							}
							p12 -= sum1 * m2 + sum2 * m1 - ind * m1 * m2;
							// r = p12 / (p1 * p2);
							ldmat_chr(i, j) = ldmat_chr(j, i) = p12 / ind;
						}
					}
				}
			}
			for (j = 0; j < chrindx.n_elem; j++){
				for(i = j; i < chrindx.n_elem; i++){
					ldmat(chrindx[i], chrindx[j]) = ldmat(chrindx[j], chrindx[i]) = ldmat_chr(i, j);
				}
			}

			chrindx = find(vecgwaschr == unichr[cc]);
			if(chrindx.n_elem > 0){
				if(verbose)	Rcpp::Rcout << "Loop on chromosome No." << cc + 1 << " with total number of SNPs in GWAS panel" << chrindx.n_elem << std::endl;
				MinimalProgressBar pb;
				Progress p(chrindx.n_elem, verbose, pb);
				arma::mat ldmat_chr(chrindx.n_elem, chrindx.n_elem);
				#pragma omp parallel for schedule(dynamic) private(j, p1, m1, sum1, i, p12, p2, m2, sum2, k, r)
				for (j = 0; j < chrindx.n_elem; j++){
					if ( ! Progress::check_abort() ) {
						p.increment();
						p1 = xx_all_gwas[chrindx[j]];
						m1 = mean_all_gwas[chrindx[j]];
						sum1 = sum_all_gwas[chrindx[j]];
						for(i = j; i < chrindx.n_elem; i++){
							p12 = 0;
							p2 = xx_all_gwas[chrindx[i]];
							m2 = mean_all_gwas[chrindx[i]];
							sum2 = sum_all_gwas[chrindx[i]];
							for(k = 0; k < indgwas; k++){
								p12 += (gwasgenomat[chrindx[i]][k]) * (gwasgenomat[chrindx[j]][k]);
							}
							p12 -= sum1 * m2 + sum2 * m1 - ind * m1 * m2;

							ldmat_chr(i, j) = ldmat_chr(j, i) = p12 / indgwas;
							// r = p12 / (p1 * p2);
							ldmat(gwasindx[chrindx[i]], gwasindx[chrindx[j]]) = ldmat(gwasindx[chrindx[j]], gwasindx[chrindx[i]]) = p12 / indgwas;
						}
					}
				}
			}
			for (j = 0; j < chrindx.n_elem; j++){
				for(i = j; i < chrindx.n_elem; i++){
					ldmat(gwasindx[chrindx[i]], gwasindx[chrindx[j]]) = ldmat(gwasindx[chrindx[j]], gwasindx[chrindx[i]]) = ldmat_chr(i, j);
				}
			}

		}
		return wrap(ldmat);
	}
}

// [[Rcpp::export]]
SEXP tXXmat_Chr_gwas(SEXP pBigMat, const NumericVector chr, SEXP gwasgeno, const NumericVector gwaschr, const LogicalVector refindx, const NumericVector gwasindx, const Nullable<double> chisq = R_NilValue, const int threads=0, const bool verbose=true){

	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()){
	case 1:
		return tXXmat_Chr_gwas<char>(xpMat, chr, gwasgeno, gwaschr, refindx, gwasindx, chisq, threads, verbose);
	case 2:
		return tXXmat_Chr_gwas<short>(xpMat, chr, gwasgeno, gwaschr, refindx, gwasindx, chisq, threads, verbose);
	case 4:
		return tXXmat_Chr_gwas<int>(xpMat, chr, gwasgeno, gwaschr, refindx, gwasindx, chisq, threads, verbose);
	case 6:
		return tXXmat_Chr_gwas<float>(xpMat, chr, gwasgeno, gwaschr, refindx, gwasindx, chisq, threads, verbose);
	case 8:
		return tXXmat_Chr_gwas<double>(xpMat, chr, gwasgeno, gwaschr, refindx, gwasindx, chisq, threads, verbose);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
