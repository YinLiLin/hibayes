#if !defined(ARMA_64BIT_WORD)
#define ARMA_64BIT_WORD 1
#endif

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec cutwind_by_bp(const arma::vec & chr, const arma::vec & pos, double bp){

	vec unichr = unique(chr);
	vec windindx(chr.n_elem);
	uvec indx1, indx2;
	int count = 1;

	for(int i = 0; i < unichr.n_elem; i++){
		double bp0 = 1;
		indx1 = find(chr == unichr[i]);
		double maxbp = max(pos(indx1));
		while(bp0 <= maxbp){
			indx2 = find(pos(indx1) >= bp0 && pos(indx1) < (bp0 + bp));
			if(indx2.n_elem != 0){
				windindx(indx1(indx2)).fill(count);
				count++;
			}
			bp0 = bp0 + bp;
		}
	}
	return windindx;
}

// [[Rcpp::export]]
arma::vec cutwind_by_num(const arma::vec & chr, const arma::vec & pos, int fixN){

	vec unichr = unique(chr);
	vec windindx(chr.n_elem);
	uvec indx1;
	int count = 1;

	for(int i = 0; i < unichr.n_elem; i++){
		indx1 = find(chr == unichr[i]);
		int chrlen = indx1.n_elem;
		if(chrlen <= fixN){
			windindx(indx1).fill(count);
			count++;
		}else{
			uvec indices = sort_index(pos(indx1));
			int st = 0;
			int end = 0;
			while(end < (chrlen - 1)){
				end = st + fixN - 1;
				if(end > (chrlen - 1))	end = (chrlen - 1);
				windindx(indx1(indices.subvec(st, end))).fill(count);
				st += fixN;
				count++;
			}
		}
	}
	return windindx;
}
