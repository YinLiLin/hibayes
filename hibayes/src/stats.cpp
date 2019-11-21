#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double uniform_sample(const double start = 0, const double end = 1){
	return R::runif(start, end);
}

// [[Rcpp::export]]
double norm_sample(const double mean = 0, const double sd = 1){
	return R::rnorm(mean, sd);
}

// [[Rcpp::export]]
double gamma_sample(const double shape, const double scale){
	return R::rgamma(shape, scale);
}

// [[Rcpp::export]]
double invgamma_sample(const double shape, const double scale){
	return (1 / gamma_sample(shape, 1 / scale));
}

// [[Rcpp::export]]
double chisq_sample(const double shape){
	return R::rchisq(shape);
}

// [[Rcpp::export]]
double invchisq_sample(const double shape, const double scale){
	return ((shape * scale) / chisq_sample(shape));
}

// [[Rcpp::export]]
double beta_sample(const double a, const double b){
	return R::rbeta(a, b);
}

// [[Rcpp::export]]
double t_sample(const double shape){
	return R::rt(shape);
}

// [[Rcpp::export]]
double cauchy_sample(const double location, const double scale){
	return R::rcauchy(location, scale);
}

// [[Rcpp::export]]
double exponential_sample(const double scale){
	return R::rexp(scale);
}

// [[Rcpp::export]]
double laplace_sample(const double mean, const double scale){
	double u = uniform_sample();
	if(u < 0.5){
		return (mean + scale * log(2 * u));
	}else{
		return (mean - scale * log(2 * (1 - u)));
	}
}

// [[Rcpp::export]]
double rinvgaussian_sample(double mu, double lambda){
    double v,z,y,x,u;
    z = R::rnorm(0,1);
    y = z * z;
    x = mu + 0.5 * mu * mu * y / lambda - 0.5 * (mu / lambda) * sqrt(4 * mu * lambda * y + mu * mu * y * y);
    u = R::runif(0,1);
    if(u <= mu / (mu + x)){
        v = x;
    }else{
        v = mu * mu / x;
    }
    return v;
}

// [[Rcpp::export]]
NumericVector rdirichlet_sample(const double n, const NumericVector x){
	NumericVector xn(x.size());
	for(int i = 0; i < x.size(); i++){
		xn[i] = gamma_sample(x[i], 1);
	}
	return (xn / sum(xn));
}
