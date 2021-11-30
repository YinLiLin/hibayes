#include "stats.h"

double uniform_sample(double start, double end){
	// return R::runif(start, end);
    return unif_rand();
}

double norm_sample(double mean, double sd){
	// return R::rnorm(mean, sd);
    return mean + sd * norm_rand();
}

double gamma_sample(double shape, double scale){
	return R::rgamma(shape, scale);
}

double invgamma_sample(double shape, double scale){
    double d = 1 / scale;
	return (1 / gamma_sample(shape, d));
}

double chisq_sample(double shape){
	return R::rchisq(shape);
}

double invchisq_sample(double shape, double scale){
	return ((shape * scale) / chisq_sample(shape));
}

double beta_sample(double a, double b){
	return R::rbeta(a, b);
}

double t_sample(double shape){
	return R::rt(shape);
}

double cauchy_sample(double location, double scale){
	return R::rcauchy(location, scale);
}

double exponential_sample(double scale){
	return R::rexp(scale);
}

double laplace_sample(double mean, double scale){
	double u = uniform_sample();
	if(u < 0.5){
		return (mean + scale * log(2 * u));
	}else{
		return (mean - scale * log(2 * (1 - u)));
	}
}

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

vec rdirichlet_sampler(int n, vec x){
	vec xn(n);
    double c = 1.0;
	for(int i = 0; i < n; i++){
		xn[i] = gamma_sample(x[i], c);
	}
	return (xn / sum(xn));
}

IntegerVector which_c(NumericVector x, double value, int c){
    
    IntegerVector eff_index_(x.size());
    int type_len = 0;
    bool logi;
    for(int i = 0; i < x.size(); i++){
        if(c == 1){
            logi = x[i] > value;
        }else if(c == 2){
            logi = x[i] >= value;
        }else if(c == 3){
            logi = x[i] < value;
        }else if(c == 4){
            logi = x[i] <= value;
        }else if(c == 5){
            logi = x[i] == value;
        }else if(c == 6){
            logi = (x[i] >= value) && (x[i] <= 1 - value);
        }else if(c == 7){
            logi = (x[i] < value) || (x[i] > 1 - value);
        }
        if(logi){
            eff_index_[type_len] = i;
            type_len++;
        }
    }
    IntegerVector eff_index(type_len);
    for(int i = 0; i < type_len; i++){
        eff_index[i] = eff_index_[i];
    }
    return eff_index;
}

