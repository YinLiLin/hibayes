#' LD variance-covariance matrix calculation
#'
#' To calculate density or sparse LD variance-covariance matrix with genotype in bigmemory format.
#'
#' @param geno the reference genotype panel in bigmemory format.
#' @param map the map information of reference genotype panel, columns are: SNPs, chromosome, physical position. 
#' @param gwas.geno (optional) the genotype of gwas samples which were used to generate the summary data.
#' @param gwas.map (optional) the map information of the genotype of gwas samples, columns are: SNPs, chromosome, physical position. 
#' @param chisq chi-squre value for generating sparse matrix, if n*r2 < chisq, it would be set to zero.
#' @param ldchr lpgical, whether to calulate the LD between chromosomes.
#' @param threads the number of threads used in computation.
#' @param verbose whether to print the information.
#'
#' @examples
#' bfile_path = system.file("extdata", "geno", package = "hibayes")
#' data = read_plink(bfile_path, out=tempfile())
#' geno = data$geno
#' map = data$map
#' \donttest{
#' xx = ldmat(geno, threads=4)   #chromosome wide full ld matrix
#' xx = ldmat(geno, chisq=5, threads=4)   #chromosome wide sparse ld matrix
#' xx = ldmat(geno, map, ldchr=FALSE, threads=4)   #chromosome block ld matrix
#' xx = ldmat(geno, map, ldchr=FALSE, chisq=5, threads=4)   #chromosome block + sparse ld matrix
#' }
#' 
#' @return
#' For full ld matrix, it returns a standard R matrix, for sparse matrix, it returns a 'dgCMatrix'.
#' 
#' @export

ldmat <- function(
	geno,
	map = NULL,
	gwas.geno = NULL,
	gwas.map = NULL,
	chisq = NULL,
	ldchr = FALSE,
	threads = 4,
	verbose = TRUE
){
	if(!is.big.matrix(geno)){
		geno <- as.big.matrix(geno)
	}
	if(!is.null(chisq)){
		if(chisq < 0)	chisq = NULL
	}
	# if(!ldchr){
	# 	if(is.null(chisq))	chisq = 0
	# }
	if(!is.null(map)){
		if(length(unique(map[,2])) == 1)	ldchr <- TRUE;
		if(!is.null(chisq)){
			if(chisq == 0 & length(unique(map[,2])) == 1)	chisq = NULL
		}
		if(length(unique(map[,1])) != nrow(map))	stop("Same SNPs names detected.")
		if(all(!is.na(map[,2]))){
			# ok
		}else{
			stop("NAs are not allowed in chromosome.")
		}
		if(sum((map[,2]) == 0) != 0){
			stop("0 is not allowed in chromosome.")
		}
		snpid <- map[, 1]
		map <- as.matrix(map[,c(2, 3)])
		suppressWarnings(max.chr <- max(as.numeric(map[, 1]), na.rm=TRUE))
		if(is.infinite(max.chr))	max.chr <- 0
		suppressWarnings(map.xy.index <- which(!as.numeric(map[, 1]) %in% c(0 : max.chr)))
		if(length(map.xy.index) != 0){
			chr.xy <- unique(map[map.xy.index, 1])
			for(i in 1:length(chr.xy)){
				map[map[, 1] == chr.xy[i], 1] <- max.chr + i
			}
		}
		map <- matrix(as.numeric(map), nrow(map))
	}else{
		if(!is.null(chisq)){
			if(chisq == 0)	chisq = NULL
		}
		ldchr = TRUE
	}
	if(!is.null(gwas.map)){
		if(length(unique(gwas.map[,1])) != nrow(gwas.map))	stop("Same SNPs names detected.")
	}
	if(is.null(gwas.geno)){
		if(ldchr){
			ldmat <- tXXmat_Geno(geno@address, chisq = chisq, threads = threads, verbose = verbose)
		}else{
			if(is.null(map))	stop("map information should be provided.")
			ldmat <- tXXmat_Chr(geno@address, chr = map[, 1], chisq = chisq, threads = threads, verbose = verbose)
		}
	}else{
		if(is.null(map))	stop("map information for reference should be provided.")
		if(is.null(gwas.map))	stop("map information for gwas sample should be provided.")
		gwassnpid <- gwas.map[, 1]
		refindx <- snpid %in% gwassnpid
		if(sum(refindx) == 0)	stop("No shared SNPs between 'geno' and 'gwas.geno'. ")
		gwasindx <- gwassnpid %in% snpid
		if(sum(gwasindx) != nrow(gwas.map)){
			gwas.geno <- deepcopy(gwas.geno, cols=which(gwasindx))
		}
		matchgwasindx <- match(gwassnpid, snpid)
		if(ldchr){
			ldmat <- tXXmat_Geno_gwas(geno@address, gwas.geno@address, refindx = refindx, gwasindx = matchgwasindx, chisq = chisq, threads = threads, verbose = verbose)
		}else{
			ldmat <- tXXmat_Chr_gwas(geno@address, chr = map[, 1], gwas.geno@address, gwaschr = map[matchgwasindx, 1], refindx = refindx, gwasindx = matchgwasindx, chisq = chisq, threads = threads, verbose = verbose)
		}

	}
	return(ldmat)
}
