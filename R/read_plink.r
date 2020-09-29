#' data load
#'
#' To load plink binary data
#'
#' @param bfile character, prefix of Plink binary format data.
#' @param maxLine number, set the number of lines to read at a time.
#' @param impute logical, whether to impute missing values in genotype by major alleles.
#' @param mode "A" or "D", additive effect or dominant effect.
#' @param out character, path and prefix of output file
#' @param threads number, the number of used threads for parallel process

#' @examples
#' bfile_path = system.file("extdata", "example", package = "hibayes")
#' data = read_plink(bfile_path, out=tempfile(), mode="A")
#' pheno = data$pheno
#' geno = data$geno
#' map = data$map

#' @export

read_plink <- 
function(
	bfile = "", 
	maxLine = 10000,
	impute = TRUE,
	mode = c("A", "D"),
	out = NULL, 
	threads = 0
){
	if(is.null(out)){
		backingfile <- paste0(basename(bfile),".bin")
		descriptorfile <- paste0(basename(bfile),".desc")
		backingpath <- "."
	}else{
		backingfile <- paste0(basename(out),".bin")
		descriptorfile <- paste0(basename(out),".desc")
		backingpath <- dirname(out)
	}
	map_file <- unlist(strsplit(descriptorfile, "", fixed = TRUE))
	sep_index <- which(map_file == ".")
	if(length(sep_index)){
		map_file <- paste0(map_file[1 : (sep_index[length(sep_index)] - 1)], collapse="")
	}else{
		map_file <- paste0(map_file, collapse="")
	}
	map_file <- paste0(backingpath, "/", map_file)
	map <- as.data.frame(rMap_c(paste0(bfile, ".bim"), out = map_file), stringsAsFactors=FALSE)
	pheno <- read.table(paste(bfile, ".fam", sep=""), header=FALSE, stringsAsFactors=FALSE)[,-c(1:5), drop=FALSE]
	colnames(pheno) <- NULL
	m <- nrow(map)
	n <- nrow(pheno)
	geno <- bigmemory::big.matrix(
		nrow = n,
		ncol = m,
		type = "char",
		dimnames = c(NULL, NULL),
		backingfile = backingfile,
		backingpath = backingpath,
		descriptorfile = descriptorfile
	)
	switch(
		match.arg(mode), 
		"A"={
			read_bed(bfile = bfile, pBigMat = geno@address, maxLine = maxLine, impt = impute, d = FALSE, threads = threads)
		},
		"D"={
			read_bed(bfile = bfile, pBigMat = geno@address, maxLine = maxLine, impt = impute, d = TRUE, threads = threads)
		}
	)
	return(list(pheno=pheno, geno=geno, map=map))	
}
