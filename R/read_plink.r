#' data load
#'
#' To load plink binary data
#'
#' @param bfile character, prefix of Plink binary format data.
#' @param maxLine number, set the number of lines to read at a time.
#' @param mode "A" or "D", additive effect or dominant effect.
#' @param backingpath the path to the directory containing the file backing cache.
#' @param descriptorfile the name of the file to hold the backingfile description, for subsequent use with ‘attach.big.matrix’; if ‘NULL’, the ‘backingfile’ is used as the root part of the descriptor file name.  The descriptor file is placed in the same directory as the backing files.
#' @param backingfile the root name for the file(s) for the cache.
#' @param threads number, the number of used threads for parallel process

#' @examples
#' bfile_path = system.file("extdata", "example", package = "hibayes")
#' data = read_plink(bfile_path, mode="A")
#' pheno = data$pheno
#' geno = data$geno
#' map = data$map

read_plink <- 
function(
	bfile = "", 
	maxLine = 10000,
	mode = c("A", "D"),
	backingpath = NULL,
	descriptorfile = NULL,
	backingfile = NULL,
	threads = 0
){
	map <- read.table(paste(bfile, ".bim", sep=""), head=FALSE)[,c(2, 1, 4)]
	pheno <- read.table(paste(bfile, ".fam", sep=""), head=FALSE)[,-c(1:5), drop=FALSE]
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
			read_bed(bfile = bfile, pBigMat = geno@address, maxLine = maxLine, d = FALSE, threads = threads)
		},
		"D"={
			read_bed(bfile = bfile, pBigMat = geno@address, maxLine = maxLine, d = TRUE, threads = threads)
		}
	)
	colnames(pheno) <- NULL
	colnames(map) <- c("SNPid", "Chr", "Pos")
	return(list(pheno=pheno, geno=geno, map=map))	
}
