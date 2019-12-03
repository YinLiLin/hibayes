#' SBayes model
#'
#' Bayes linear regression model using summary level data
#'
#' @param sumstat matrix of summary data, details refer to https://cnsgenomics.com/software/gcta/#COJO.
#' @param ldm dense or sparse matrix, ld for reference panel (m * m, m is the number of SNPs), note that the order of SNPs should be consistent with summary data.
#' @param model sbayes model including: "SBayesRR", "SBayesA", "SBayesLASSO", "SBayesB", "SBayesBpi", "SBayesC", "SBayesCpi", "SBayesR", "CG".
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
#' bfile_path = system.file("extdata", "example", package = "hibayes")
#' data = read_plink(bfile_path)
#' geno = data$geno
#' map = data$map
#' head(map)
#' sumstat_path = system.file("extdata", "example.ma", package = "hibayes")
#' sumstat = read.table(sumstat_path, header=TRUE)
#' head(sumstat)
#' 
#' ## computate ld variance covariance matrix
#' # ldm1 = ldmat(geno, threads=4)   #chromosome wide full ld matrix
#' # ldm2 = ldmat(geno, map, ldchr=FALSE, chisq=5, threads=4)   #chromosome block + sparse ld matrix
#' 
#' ## if the order of SNPs in genotype is not consistent with the order in sumstat file, prior adjusting is necessary.
#' # indx = match(sumstat[, 1], map[, 1])
#' # ldm1 = ldm1[indx, indx]
#' 
#' ## fit model
#' # fit = sbayes(sumstat=sumstat, ldm=ldm1, model="SBayesR")

sbayes <- 
function(
	sumstat,
    ldm,
    model = c("SBayesB", "SBayesA", "SBayesLASSO", "SBayesRR", "SBayesBpi", "SBayesC", "SBayesCpi", "SBayesR", "CG"),
    map = NULL,
    pi = NULL,
    lambda = NULL,
    fold = NULL,
    niter = 20000,
    nburn = 12000,
    windsize = NULL,
    wppa = 0.01,
    vara = NULL,
    dfvara = NULL,
    s2vara = NULL,
    vare = NULL,
    dfvare = NULL,
    s2vare = NULL,
    outfreq = 10,
    seed = 666666,
    verbose = TRUE
){
	if(is.matrix(ldm)){
		sparse = FALSE
	}else if(is(ldm, "dgCMatrix")){
		sparse = TRUE
	}else{
		stop("Unrecognized type of ldm.")
	}
	if(!is.null(windsize)){
		if(is.null(map)){
			stop("map information must be provided.")
		}else{
			if(ncol(map) < 3)	stop("At least 3 columns in map.")
			map <- map[, c(1:3)]
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
			if(sparse){
				res = SBayesRR_spa(sumstat=sumstat, ldm=ldm, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}else{
				res = SBayesRR_den(sumstat=sumstat, ldm=ldm, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}
		},
		"SBayesA"={
			if(sparse){
				res = SBayesA_spa(sumstat=sumstat, ldm=ldm, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}else{
				res = SBayesA_den(sumstat=sumstat, ldm=ldm, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}
		},
		"SBayesLASSO"={
			if(sparse){
				res = SBayesLASSO_spa(sumstat=sumstat, ldm=ldm, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)	
			}else{
				res = SBayesLASSO_den(sumstat=sumstat, ldm=ldm, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}
		},
		"SBayesB"={
			if(sparse){
				res = SBayesB_spa(sumstat=sumstat, ldm=ldm, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}else{
				res = SBayesB_den(sumstat=sumstat, ldm=ldm, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}
		},
		"SBayesBpi"={
			if(sparse){
				res = SBayesBpi_spa(sumstat=sumstat, ldm=ldm, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}else{
				res = SBayesBpi_den(sumstat=sumstat, ldm=ldm, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}
		},
		"SBayesC"={
			if(sparse){
				res = SBayesC_spa(sumstat=sumstat, ldm=ldm, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}else{
				res = SBayesC_den(sumstat=sumstat, ldm=ldm, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}
		},
		"SBayesCpi"={
			if(sparse){
				res = SBayesCpi_spa(sumstat=sumstat, ldm=ldm, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}else{
				res = SBayesCpi_den(sumstat=sumstat, ldm=ldm, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}
		},
		"SBayesR"={
			if(sparse){
				res = SBayesR_spa(sumstat=sumstat, ldm=ldm, pi=pi, fold=fold, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}else{
				res = SBayesR_den(sumstat=sumstat, ldm=ldm, pi=pi, fold=fold, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
			}
		},
		"CG"={
			if(!is.null(lambda)){
				if(length(lambda) == 1){
					lambda = rep(lambda, nrow(sumstat))
				}else if(length(lambda) != nrow(sumstat)){
					stop("length of lambda should equal to the length of SNPs.")
				}
			}
			if(sparse){
				res = conjgt_spa(sumstat=sumstat, ldm=ldm, lambda = lambda, outfreq=outfreq, verbose=verbose)
			}else{
				res = conjgt_den(sumstat=sumstat, ldm=ldm, lambda = lambda, outfreq=outfreq, verbose=verbose)
			}
		}
	)

	if(!is.null(windsize) & match.arg(model) != "CG"){
		WPPA <- res$wppa
		WGVE <- res$wgve
		res <- head(res, length(res) - 2)
		res$gwas <- data.frame(windinfo, WPPA, WGVE)
	}
	return(res)
}
