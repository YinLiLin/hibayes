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
#' bfile_path = system.file("extdata", "demo", package = "hibayes")
#' data = read_plink(bfile_path, out=tempfile(), mode="A")
#' fam = data$fam
#' geno = data$geno
#' map = data$map

#' @return
#' four files will be generated in the directed folder: *.desc, *.bin, *.id, *.map, where '*' is the prefix of the argument 'out', the memory-mapping files can be fast loaded into memory by 'geno = attach.big.matrix("*.desc")'. Note that hibayes will code the genotype A1A1 as 2, A1A2 as 1, and A2A2 as 0, where A1 is the first allele of each marker in *.bim file, therefore the estimated effect size is on A1 allele, users should pay attention to it when a process involves marker effect.

#' @export

read_plink <- 
function(
	bfile = "", 
	maxLine = 10000,
	impute = TRUE,
	mode = c("A", "D"),
	out = NULL, 
	threads = 4
){
	if(is.null(out))	out <- tempfile()
	# {
	# 	backingfile <- paste0(basename(bfile),".bin")
	# 	descriptorfile <- paste0(basename(bfile),".desc")
	# 	backingpath <- "."
	# }else{
		backingfile <- paste0(basename(out),".bin")
		descriptorfile <- paste0(basename(out),".desc")
		backingpath <- dirname(out)
	# }
	# map_file <- unlist(strsplit(descriptorfile, "", fixed = TRUE))
	# sep_index <- which(map_file == ".")
	# if(length(sep_index)){
	# 	map_file <- paste0(map_file[1 : (sep_index[length(sep_index)] - 1)], collapse="")
	# }else{
	# 	map_file <- paste0(map_file, collapse="")
	# }
	# map_file <- paste0(backingpath, "/", map_file)
	map_file <- out
	map <- as.data.frame(rMap_c(paste0(bfile, ".bim"), out = map_file), stringsAsFactors=FALSE)
	fam <- read.table(paste(bfile, ".fam", sep=""), header=FALSE, stringsAsFactors=FALSE)
	colnames(fam) <- NULL
	m <- nrow(map)
	n <- nrow(fam)
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
	write.table(fam[, 2], paste0(out, ".id"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
	return(list(fam=fam, geno=geno, map=map))	
}
