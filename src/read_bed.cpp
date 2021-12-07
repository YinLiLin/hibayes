#include <Rcpp.h>
#include "omp_set.h"
#include <fstream>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(bigmemory, BH)]]

using namespace std;
using namespace Rcpp;

int FileNrow(std::string filename) {
    // Define
    string line;
    int m = 0; 
    ifstream file(filename);
    
    if (!file) throw Rcpp::exception(("Error: can not open the file [" + filename + "].").c_str());

    while (getline(file, line))
        m++;
    
    file.close();
    return m;
}

// [[Rcpp::export]]
List rMap_c(std::string map_file, const Nullable<std::string> out = R_NilValue){

	bool fileout;
	std::string filename;
	if(out.isNotNull()){
		fileout = true;
		filename = as<std::string>(out);
	}else{
		fileout = false;
	}

	int n = FileNrow(map_file);
	string line;
	vector<string> l;
	vector<string> snp(n);
	vector<string> chr(n);
	vector<string> pos(n);
	vector<string> a1(n);
	vector<string> a2(n);

	ifstream file(map_file);
	if (!file) throw Rcpp::exception(("Error: can not open the file [" + map_file + "].").c_str());

	int idx = 0;
	string s, r, a, c, g, p;
	// float g, p;
	while (file >> c >> s >> g >> p >> r >> a) {

		snp[idx] = s;
		chr[idx] = c;
		pos[idx] = p;
		a1[idx] = r;
		a2[idx] = a;

		idx++;
	}
	// while (getline(file, line)) {

	// 	boost::split(l, line, boost::is_any_of("\t ,"));
	// 	// boost::split(l, line, "\t");

	// 	snp[idx] = l[1];
	// 	chr[idx] = l[0];
	// 	pos[idx] = l[3];
	// 	a1[idx] = l[4];
	// 	a2[idx] = l[5];

	// 	idx++;
	// }
	file.close();

	if(fileout){
		ofstream map(filename + ".map");
		map << "SNP\tCHROM\tPOS\tA1\tA2" << endl;
		for(int i = 0; i < n; i++){
			map << snp[i] << "\t" << chr[i] << "\t" << pos[i] << "\t" << a1[i] << "\t" << a2[i] << endl;
		}
		map.close();
	}

	return List::create(Named("SNP") = snp,
                       Named("Chr") = chr,
                       Named("Pos") = pos,
					   Named("A1") = a1,
					   Named("A2") = a2
	);
}

template <typename T>
void read_bed(std::string bed_file, XPtr<BigMatrix> pMat, const long maxLine, const double NA_C, const bool impt = true, const bool d=false, const int threads=0){
	string ending = ".bed";
	if (bed_file.length() <= ending.length() ||
		0 != bed_file.compare(bed_file.length() - ending.length(), ending.length(), ending))
		bed_file += ending;

	omp_setup(threads);

	long n = pMat->nrow() / 4;  // 4 individual = 1 bit
	int m = pMat->ncol();
	int nid = pMat->nrow();
	if (pMat->nrow() % 4 != 0) 
		n++; 
	char *buffer;
	long buffer_size;
	MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
	
	// map
	std::map<int, T> code;
	code[3] = d ? 0 : 0;
	code[2] = d ? 1 : 1;
	code[1] = static_cast<T>(NA_C);
	code[0] = d ? 0 : 2;

    std::vector<size_t> ggvec(3);
    if(d){
    	ggvec = {0, 1, 0};
    }else{
    	ggvec = {0, 1, 2};
    }

	// open file
	FILE *fin;
	fin = fopen(bed_file.c_str(), "rb");
	fseek(fin, 0, SEEK_END);
	long length = ftell(fin);
	rewind(fin);
	
	// get buffer_size
	buffer_size = maxLine > 0 ? (maxLine * n) : (length - 3);
	
	int n_block = (length - 3) / buffer_size;
	if ((length - 3) % buffer_size != 0) { n_block++; }
	buffer = new char [3];
	size_t n_bytes_read = static_cast<size_t>(fread(buffer, 1, 3, fin));
	
	size_t r, c, x;
	uint8_t p;
	NumericVector miss(m); miss.fill(0);

	for (int i = 0; i < n_block; i++) {
		buffer = new char [buffer_size];
		n_bytes_read = static_cast<size_t>(fread(buffer, 1, buffer_size, fin));
		
		int cond = min(buffer_size, (length - 3 - i * buffer_size));

		#pragma omp parallel for schedule(dynamic) private(r, c, p, x)
		for (size_t j = 0; j < cond; j++) {
			// bit -> item in matrix
			r = (i * buffer_size + j) / n;
			c = (i * buffer_size + j) % n * 4;
			p = buffer[j];
			for (x = 0; x < 4 && (c + x) < pMat->nrow(); x++) {
				T gg = code[(p >> (2*x)) & 0x03];
				mat[r][c + x] = gg;
				if(gg == NA_C){
					miss[r] = 1;
				}
			}
		}
	}
	fclose(fin);


	int NMISS = 0;
	for (int i = 0; i < m; i++) {
		if(miss[i]){
			Rcout << "Warning: Missing values in genotype exist!" << endl;
			break;
		}else{
			NMISS++;
		}
	}

	if(impt && (NMISS != m)){

		Rcout << "Imputing missing values by major genotype..." << endl;

	 	// impute
	 	#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < m; i++) {
			if(miss[i]){
		        std::vector<size_t> na_index = {};
		        std::vector<size_t> counts(3);
		        // size_t counts[3] = { 0 };

		        // count allele, record missing index
		       	int max = 0;
		        T major = 0;

		        if(d){
		        	for (int j = 0; j < nid; j++) {
			            switch(int(mat[i][j])) {
				            case 0: counts[0]++; break;
				            case 1: counts[1]++; break;
				            default: na_index.push_back(j);
		            	}
		        	}	
		        }else{
		        	for (int j = 0; j < nid; j++) {
			            switch(int(mat[i][j])) {
				            case 0: counts[0]++; break;
				            case 1: counts[1]++; break;
				            case 2: counts[2]++; break;
				            default: na_index.push_back(j);
		            	}
		        	}
		        }

	        	for(size_t j = 0; j < counts.size(); j++){
	        		if(counts[j] > max){
	        			max = counts[j];
	        			major = ggvec[j];
	        		}
	        	}

	        	// impute
		        for (auto&& x: na_index) {
		            mat[i][x] = major;   
		        }
		    }
		}
	}
	return;
}

// [[Rcpp::export]]
void read_bed(std::string bfile, const SEXP pBigMat, const long maxLine, const bool impt = true, const bool d=false, const int threads=0){
	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return read_bed<char>(bfile, xpMat, maxLine, NA_CHAR, impt, d, threads);
	case 2:
		return read_bed<short>(bfile, xpMat, maxLine, NA_SHORT, impt, d, threads);
	case 4:
		return read_bed<int>(bfile, xpMat, maxLine, NA_INTEGER, impt, d, threads);
	case 8:
		return read_bed<double>(bfile, xpMat, maxLine, NA_REAL, impt, d, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
