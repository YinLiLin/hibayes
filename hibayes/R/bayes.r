#' Bayes model
#'
#' Bayes linear regression model using individual level data
#'
#' @param y vector of phenotype, NAs are not allowed.
#' @param X numeric matrix of genotype with individuals in rows and markers in columns, NAs are not allowed.
#' @param model bayes model including: "BayesRR", "BayesA", "BayesLASSO", "BayesB", "BayesBpi", "BayesC", "BayesCpi", "BayesR".
#' @param map (optional, only for GWAS) the map information of genotype, columns are: SNPs, chromosome, physical position. 
#' @param pi percentage of zero effect SNPs. For bayesR, it is a vector for groups of SNPs, the default is c(0.95, 0.02, 0.02, 0.01).
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
#' #' need to update

bayes <- 
function(
	y,
    X,
    model = c("BayesRR", "BayesA", "BayesLASSO", "BayesB", "BayesBpi", "BayesC", "BayesCpi", "BayesR"),
    map = NULL,
    pi = NULL,
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
		if(match.arg(model) == "BayesR"){
			pi <- c(0.95, 0.02, 0.02, 0.01)
		}else{
			pi <- 0.95
		}
	}
	X <- as.matrix(X); gc()
	switch(
		match.arg(model), 
		"BayesRR"={
			res = BayesRR(y=y, X=X, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"BayesA"={
			res = BayesA(y=y, X=X, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"BayesLASSO"={
			res = BayesLASSO(y=y, X=X, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"BayesB"={
			res = BayesB(y=y, X=X, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"BayesBpi"={
			res = BayesBpi(y=y, X=X, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"BayesC"={
			res = BayesC(y=y, X=X, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"BayesCpi"={
			res = BayesCpi(y=y, X=X, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		},
		"BayesR"={
			res = BayesR(y=y, X=X, pi=pi, fold=fold, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vara=vara, dfvara=dfvara, s2vara=s2vara, vare=vare, dfvare=dfvare, s2vare=s2vare, outfreq=outfreq, verbose=verbose)
		}
	)

	if(!is.null(windsize)){
		wppa <- res$wppa
		wgve <- res$wgve
		res <- head(res, length(res) - 2)
		res$gwas <- data.frame(windinfo, wppa, wgve)
	}
	return(res)
}
