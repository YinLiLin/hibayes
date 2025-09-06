#include "hibayes.h"
#include "solver.h"

// [[Rcpp::export]]
SEXP make_grm(
    arma::mat &Z,
    double lambda = 0.0,
    bool inverse = false,
    bool eigen = false,
    bool verbose = true
){

    blas_int n = Z.n_rows;
    blas_int m = Z.n_cols;
    if(verbose) Rcpp::Rcout << "Start construct G matrix for " <<  Z.n_rows << " individuals using " << m << " markers" << std::endl;
    
    if(verbose) Rcpp::Rcout << "Calculate mean for all markers" << std::endl;
    rowvec means = mean(Z);

    if(verbose) Rcpp::Rcout << "Center genotype matrix" << std::endl;
    Z.each_row() -= means;

    if(verbose) Rcpp::Rcout << "Compute Z * Z'" << std::endl;
	double alp = 1.0;
	double beta = 1.0;
	char uplo = 'L';
    arma::mat G = zeros<arma::mat>(n, n);
	double* C = G.memptr();
	double* A = Z.memptr();
    // G = Z * Z.t();
	dsyrk_(&uplo, "N", &n, &m, &alp, A, &n, &beta, C, &n);
    for (size_t j = 0; j < n; j++) {
		for (size_t i = (j + 1); i < n; i++) {
			G(j, i) = G(i, j);
		}
	}
    G /= mean(G.diag());

    if(inverse){
        if(verbose) Rcpp::Rcout << "Invert the G matrix" << std::endl;
        solve(G, lambda);
    }

    if(eigen){
        if(verbose) Rcpp::Rcout << "Eigen decomposition on G matrix" << std::endl;
        if(lambda)  G.diag() += lambda;
        arma::vec ev;
        eigen_sym_dc(G, ev);
        return List::create(wrap(ev), wrap(G));
    }else{
        return wrap(G);
    }
}

// [[Rcpp::export]]
List make_ped(
    std::vector<std::string> &pvec,
    std::vector<std::string> &svec,
    std::vector<std::string> &dvec,
    bool verbose = true
)
{
    if(verbose) Rcpp::Rcout << "Pre-adjust pedigree order" << std::endl;
    set<string> na_pool = {"NA", "Na", ".", "-", "NaN", "NAN", "nan", "na", "N/A", "n/a", "<NA>"};
    vector< vector<string> > ped_mat(3);
    string p, s, d;

    for(int i = 0; i < pvec.size(); i++){
        p = pvec[i];
        s = svec[i];
        d = dvec[i];
        if(na_pool.find(p) == na_pool.end()){
            ped_mat[0].push_back(p);
            if(na_pool.find(s) != na_pool.end())    s = "0";
                ped_mat[1].push_back(s);
            if(na_pool.find(d) != na_pool.end())    d = "0";
                ped_mat[2].push_back(d);
        }
    }
    
    set<std::string> ped_1(ped_mat[0].begin(), ped_mat[0].end());
    if(ped_mat[0].size() != ped_1.size())   throw Rcpp::exception("repeated records are not allowed in the first column of pedigree file.");

    int n = (ped_mat[0]).size();
    uvec rowlogi(n); rowlogi.fill(1);
    vector<string> ss, dd;
    set<std::string> idset;
    idset.insert("0");
    vector<string> id;
    for(int i = 0; i < n; i++){
        if(ped_mat[1][i] == "0" && ped_mat[2][i] == "0"){
            idset.insert(ped_mat[0][i]);
            id.push_back(ped_mat[0][i]);
            ss.push_back("0");
            dd.push_back("0");
            rowlogi[i] = 0;
        }else{
            if(ped_mat[1][i] != "0" && ped_1.find(ped_mat[1][i]) == ped_1.end() && idset.find(ped_mat[1][i]) == idset.end()){
                idset.insert(ped_mat[1][i]);
                id.push_back(ped_mat[1][i]);
                ss.push_back("0");
                dd.push_back("0"); 
            }
            if(ped_mat[2][i] != "0" && ped_1.find(ped_mat[2][i]) == ped_1.end() && idset.find(ped_mat[2][i]) == idset.end()){
                idset.insert(ped_mat[2][i]);
                id.push_back(ped_mat[2][i]);
                ss.push_back("0");
                dd.push_back("0"); 
            }
        }
    }
    
    int loopin = 0;
    while(sum(rowlogi)){
        for(int i = 0; i < n; i++){
            if(!rowlogi[i]) continue;
            if(idset.find(ped_mat[1][i]) != idset.end() && idset.find(ped_mat[2][i]) != idset.end()){
                idset.insert(ped_mat[0][i]);
                id.push_back(ped_mat[0][i]);
                ss.push_back(ped_mat[1][i]);
                dd.push_back(ped_mat[2][i]);
                rowlogi[i] = 0;
                loopin++;
            }
        }
        if(!loopin){
            for(int i = 0; i < n; i++){
                if(!rowlogi[i]) continue;
                if(idset.find(ped_mat[1][i]) != idset.end() || idset.find(ped_mat[2][i]) != idset.end()){
                    idset.insert(ped_mat[0][i]);
                    id.push_back(ped_mat[0][i]);
                    ss.push_back(ped_mat[1][i]);
                    dd.push_back(ped_mat[2][i]);
                    rowlogi[i] = 0;
                    loopin++;
                }
            }
        }
        if(!loopin){
            for(int i = 0; i < n; i++){
                if(!rowlogi[i]) continue;
                idset.insert(ped_mat[0][i]);
                id.push_back(ped_mat[0][i]);
                ss.push_back(ped_mat[1][i]);
                dd.push_back(ped_mat[2][i]);
                rowlogi[i] = 0;
            }
        }
        loopin = 0;
    }
    if(id.size() == 0)  throw Rcpp::exception("no individuals detected;");
    if(verbose) Rcpp::Rcout << id.size() << " unique individuals have been detected in pedigree" << std::endl;

    // position + 1
    vector<int> ints;
    vector<int> intd;
    map<string, int> ped_map;
    ped_map.insert(pair<string, int>("0", 0));
    for (int j = 0; j < id.size(); j++){
        ped_map.insert(pair<string, int>(id[j], j + 1));
    }
    map<string, int>::iterator iter;
    for(int i = 0; i < id.size(); i++){
        iter = ped_map.find(ss[i]);
        ints.push_back(iter->second);
        iter = ped_map.find(dd[i]);
        intd.push_back(iter->second);
    }
    return List::create(id, ints, intd);
}

// [[Rcpp::export]]
SEXP make_Ainv(
    std::vector<int> s,
    std::vector<int> d,
    bool verbose = true
) {
    
    if(verbose) Rcpp::Rcout << "Derive inverse of A matrix from pedigree" << std::endl;
    int n = s.size();
    arma::sp_mat Ainv(n, n);

    for (int x = 0; x < n; x++) {
        int sx = s[x] - 1;
        int dx = d[x] - 1;
        if (s[x] == 0 && d[x] == 0) {
            Ainv(x, x) = 1;
        } else if (s[x] > 0 && d[x] > 0) {
            Ainv( x,  x) = Ainv( x,  x) + 2;
            Ainv( x, sx) = Ainv(sx,  x) = (Ainv( x, sx) - 1);
            Ainv(dx,  x) = Ainv( x, dx) = (Ainv(dx,  x) - 1);
            Ainv(sx, sx) = Ainv(sx, sx) + 0.5;
            Ainv(sx, dx) = Ainv(dx, sx) = (Ainv(sx, dx) + 0.5);
            Ainv(dx, dx) = Ainv(dx, dx) + 0.5;
        } else if (s[x] > 0 && d[x] == 0) {
            Ainv( x,  x) = Ainv( x,  x) + (4/3);
            Ainv( x, sx) = Ainv(sx,  x) = (Ainv( x, sx) - 2/3);
            Ainv(sx, sx) = Ainv(sx, sx) + (1/3);
        } else {    // s[x] == 0 && d[x] > 0
            Ainv( x,  x) = Ainv( x,  x) + (4/3);
            Ainv( x, dx) = Ainv(dx,  x) = (Ainv( x, dx) - 2/3);
            Ainv(dx, dx) = Ainv(dx, dx) + (1/3);
        }
    }
    return wrap(Ainv);
}

// [[Rcpp::export]]
SEXP geno_impute(
    arma::sp_mat &Ang,
    arma::mat &geno,
    int threads
){
    omp_setup(threads);
    int m = geno.n_cols;
    int n = Ang.n_rows;
    arma::mat imp_geno(n, m);

    if(threads == 1){
        imp_geno = Ang * geno;
    }else{

        #pragma omp parallel for
        for(int i = 0; i < m; i++){
            imp_geno.col(i) = Ang * geno.col(i);
        }
    }
    return wrap(imp_geno);
}
