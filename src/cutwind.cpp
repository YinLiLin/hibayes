// #if !defined(ARMA_64BIT_WORD)
// #define ARMA_64BIT_WORD 1
// #endif

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec cutwind(const arma::vec & chr, const arma::vec & pos, double bp){

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
