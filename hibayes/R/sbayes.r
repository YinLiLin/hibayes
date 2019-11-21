#' SBayes model
#'
#' Bayes linear regression model using summary level data
#'
#' @param sumstat matrix of summary data, details refer to https://cnsgenomics.com/software/gcta/#COJO.
#' @param ldm dense or sparse matrix, ld for reference panel (m * m, m is the number of SNPs), note that the order of SNPs should be consistent with summary data.
#' @param model sbayes model including: "SBayesRR", "SBayesA", "SBayesLASSO", "SBayesB", "SBayesBpi", "SBayesC", "SBayesCpi", "SBayesR".
#' @param map (optional, only for GWAS) the map information of genotype, columns are: SNPs, chromosome, physical position. 
#' @param pi percentage of zero effect SNPs. For bayesR, it is a vector for groups of SNPs, the default is c(0.95, 0.02, 0.02, 0.01).
#' @param lambda value or vector, the ridge regression value for each SNPs.
#' @param fold percentage of variance explained for groups of SNPs, the default is c(0, 0.0001, 0.001, 0.01).
#' @param niter the number of MCMC iteration.
#' @param nburn the number of iterations to be discarded.
#' @param windsize window size in bp for GWAS, the default is 1e6.
#' @param wppa the threshold of genetic variance explained by single window, the default is 0.01.
#' @param vara prior value of genetic variance.
#' @param dfvara the number of degrees of freedom for the distribution of genetic variance. 
#' @param s2vara scale parameter for the distribution of genetic variance.
#' @param vare prior value of residual variance.
#' @param dfvare the number of degrees of freedom for the distribution of residual variance.
#' @param s2vare scale parameter for the distribution of residual variance.
#' @param outfreq frequency of information output on console, the default is 100.
#' @param seed seed for random sample.
#' @param verbose whether to print the iteration information.

#' @examples
#' need to update

sbayes <- 
function(
	sumstat,
    ldm,
    model = c("SBayesB", "SBayesA", "SBayesLASSO", "SBayesRR", "SBayesBpi", "SBayesC", "SBayesCpi", "SBayesR", "CG"),
    map = NULL,
    pi = NULL,
    lambda = NULL,
    fold = NULL,
    niter = 50000,
    nburn = 25000,
    windsize = NULL,
    wppa = 0.01,
    vara = NULL,
    dfvara = NULL,
    s2vara = NULL,
    vare = NULL,
    dfvare = NULL,
    s2vare = NULL,
    outfreq = 100,
    seed = 666666,
    verbose = TRUE
){
	if(!is.null(windsize)){
		if(is.null(map)){
			stop("map information must be provided.")
		}else{
			if(ncol(map) < 3)	stop("At least 3 columns in map.")
		}
		if(all(!is.na(map[,2]))){
			# ok
		}else{
			stop("NAs are not allowed in chromosome.")
		}
		if(sum((map[,2]) == 0) != 0){
			stop("0 is not allowed in chromosome.")
		}
		if(all(!is.na(map[,3]))){
			# ok
		}else{
			stop("NAs are not allowed in physical position.")
		}
		if(sum((map[,3]) == 0) != 0){
			stop("0 is not allowed in physical position.")
		}
		if(!is.numeric(map[,3]))	stop("Factors or characters are not allowed in physical position.")
		map <- as.matrix(map[,-1])
		chr <- map[, 1]
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
		chr <- chr[order(map[,1])]
		if(max(map[,2]) < windsize)	stop("Maximum of physical position is less than wind size.")
		windindx <- cutwind(map[,1], map[,2], windsize)
		windrange <- do.call(rbind, tapply(map[, 2], windindx, range))
		windsnpN <- tapply(map[, 2], windindx, length)
		windchr <- unique(chr)[tapply(map[, 1], windindx, unique)]
		windinfo <- data.frame(paste0("wind", 1:max(windindx)), windchr, windsnpN, windrange)
		colnames(windinfo) <- c("WIND", "CHR", "NUM", "START", "END")
	}else{
		windindx <- NULL
	}
	set.seed(seed)
	if(is.null(pi)){
		if(match.arg(model) == "SBayesR"){
			pi <- c(0.95, 0.02, 0.02, 0.01)
		}else{
			pi <- 0.95
		}
	}
	if(ncol(sumstat) != 8)	stop("Inappropriate summary data format.")
	sumstat <- sumstat[, c(4, 5, 6, 8)]
	sumstat <- data.matrix(sumstat)
	switch(
		match.arg(model), 
		"SBayesRR"={
			res = SBayesRR(sumstat=sumstat, ldm=ldm, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"SBayesA"={
			res = SBayesA(sumstat=sumstat, ldm=ldm, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"SBayesLASSO"={
			res = SBayesLASSO(sumstat=sumstat, ldm=ldm, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"SBayesB"={
			res = SBayesB(sumstat=sumstat, ldm=ldm, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"SBayesBpi"={
			res = SBayesBpi(sumstat=sumstat, ldm=ldm, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"SBayesC"={
			res = SBayesC(sumstat=sumstat, ldm=ldm, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"SBayesCpi"={
			res = SBayesCpi(sumstat=sumstat, ldm=ldm, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"SBayesR"={
			res = SBayesR(sumstat=sumstat, ldm=ldm, pi=pi, fold=fold, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"CG"={
			if(!is.null(lambda)){
				if(length(lambda) == 1){
					lambda = rep(lambda, nrow(sumstat))
				}else if(length(lambda) != nrow(sumstat)){
					stop("length of lambda should equal to the length of SNPs.")
				}
			}
			res = conjgt(sumstat=sumstat, ldm=ldm, lambda = lambda, outfreq=outfreq, verbose=verbose)
		}
	)

	if(!is.null(windsize) & match.arg(model) != "CG"){
		wppa <- res$wppa
		wgve <- res$wgve
		res <- head(res, length(res) - 2)
		res$gwas <- data.frame(windinfo, wppa, wgve)
	}
	return(res)
}

#' LD calculation
#'
#' To calculate density or sparse LD matrix with genotype in bigmemory format.
#'
#' @param geno the reference genotype panel in bigmemory format.
#' @param map the map information of reference genotype panel, columns are: SNPs, chromosome, physical position. 
#' @param gwas.geno (optional) the genotype of gwas samples which were used to generate the summary data.
#' @param gwas.map (optional) the map information of the genotype of gwas samples, columns are: SNPs, chromosome, physical position. 
#' @param chisq chi-squre value for generating sparse matrix, if n*r2 < chsiq, it would be set to zero.
#' @param ldchr lpgical, whether to calulate the LD between chromosomes.
#' @param threads the number of threads used in computation.
#' @param verbose whether to print the information.

#' @examples
#' #' need to update

ldmat <- function(
	geno,
	map,
	gwas.geno = NULL,
	gwas.map = NULL,
	chisq = 0,
	ldchr = FALSE,
	threads = NULL,
	verbose = TRUE
){
	if(!is.null(map)){
		if(length(unique(map[,2])) == 1)	ldchr <- TRUE;
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
		map <- as.matrix(map[,-1])
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
