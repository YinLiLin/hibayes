#' data load
#'
#' To load plink binary data
#'
#' @param bfile character, prefix of Plink binary format data.
#' @param maxLine number, set the number of lines to read at a time.
#' @param mode "A" or "D", additive effect or dominant effect.
#' @param threads number, the number of used threads for parallel process

#' @examples
#' need to update

read_plink <- 
function(
	bfile = "", 
	maxLine = 10000,
	mode = c("A", "D"),
	threads = 1
){
	map <- read.table(paste(bfile, ".bim", sep=""), head=FALSE)[,c(2, 1, 4)]
	pheno <- read.table(paste(bfile, ".fam", sep=""), head=FALSE)[,-c(1:5), drop=FALSE]
	m <- nrow(map)
	n <- nrow(fam)
	geno <- bigmemory::filebacked.big.matrix(
		nrow = n,
		ncol = m,
		type = "char",
		dimnames = c(NULL, NULL)
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
	return(list(pheno, geno, map))	
}
