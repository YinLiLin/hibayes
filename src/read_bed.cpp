#include <Rcpp.h>
#include <omp.h>
#include <iostream>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(bigmemory, BH)]]

using namespace std;
using namespace Rcpp;

template <typename T>
void read_bed(std::string bed_file, XPtr<BigMatrix> pMat, const long maxLine, const double NA_C, const bool d=false, const int threads=0){
	string ending = ".bed";
	if (bed_file.length() <= ending.length() ||
		0 != bed_file.compare(bed_file.length() - ending.length(), ending.length(), ending))
		bed_file += ending;

	if (threads == 0) {
		omp_set_num_threads(omp_get_num_procs());
	}else if(threads > 0) {
		omp_set_num_threads(threads);
	}

	long n = pMat->nrow() / 4;  // 4 individual = 1 bit
	if (pMat->nrow() % 4 != 0) 
		n++; 
	char *buffer;
	long buffer_size;
	MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
	
	// map
	std::map<int, T> code;
	if(d){
		code[3] = 0;
		code[2] = 1;
		code[1] = static_cast<T>(NA_C);
		code[0] = 0;
	}else{
		code[3] = 0;
		code[2] = 1;
		code[1] = static_cast<T>(NA_C);
		code[0] = 2;	
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
	fread(buffer, 1, 3, fin);
	
	size_t r, c, x;
	uint8_t p;

	for (int i = 0; i < n_block; i++) {
		buffer = new char [buffer_size];
		fread(buffer, 1, buffer_size, fin);

		int cond = min(buffer_size, (length - 3 - i * buffer_size));

		#pragma omp parallel for schedule(dynamic) private(r, c, p, x)
		for (size_t j = 0; j < cond; j++) {
			// bit -> item in matrix
			r = (i * buffer_size + j) / n;
			c = (i * buffer_size + j) % n * 4;
			p = buffer[j];
			for (x = 0; x < 4 && (c + x) < pMat->nrow(); x++) {
				mat[r][c + x] = code[(p >> (2*x)) & 0x03];
			}
		}
	}
	fclose(fin);
	return;
}

// [[Rcpp::export]]
void read_bed(std::string bfile, const SEXP pBigMat, const long maxLine, const bool d=false, const int threads=0){
	XPtr<BigMatrix> xpMat(pBigMat);

	switch(xpMat->matrix_type()) {
	case 1:
		return read_bed<char>(bfile, xpMat, maxLine, NA_CHAR, d, threads);
	case 2:
		return read_bed<short>(bfile, xpMat, maxLine, NA_SHORT, d, threads);
	case 4:
		return read_bed<int>(bfile, xpMat, maxLine, NA_INTEGER, d, threads);
	case 8:
		return read_bed<double>(bfile, xpMat, maxLine, NA_REAL, d, threads);
	default:
		throw Rcpp::exception("unknown type detected for big.matrix object!");
	}
}
