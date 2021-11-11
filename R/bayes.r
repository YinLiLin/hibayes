#' Bayes model
#'
#' Bayes linear regression model using individual level data
#'
#' @param y vector of phenotype, NAs are not allowed.
#' @param Z numeric matrix of genotype with individuals in rows and markers in columns, NAs are not allowed.
#' @param X (optional) covariate matrix of all individuals, all values should be in digits, characters and NAs are not allowed.
#' @param model bayes model including: "BayesRR", "BayesA", "BayesLASSO", "BayesB", "BayesBpi", "BayesC", "BayesCpi", "BayesR".
#' @param map (optional, only for GWAS) the map information of genotype, columns are: SNPs, chromosome, physical position. 
#' @param pi percentage of zero effect SNPs. For bayesR, it is a vector for groups of SNPs, the default is c(0.95, 0.02, 0.02, 0.01).
#' @param fold percentage of variance explained for groups of SNPs, the default is c(0, 0.0001, 0.001, 0.01).
#' @param niter the number of MCMC iteration.
#' @param nburn the number of iterations to be discarded.
#' @param windsize window size in bp for GWAS, the default is NULL.
#' @param wppa the threshold of genetic variance explained by single window, the default is 0.01.
#' @param vg prior value of genetic variance.
#' @param dfvg the number of degrees of freedom for the distribution of genetic variance. 
#' @param s2vg scale parameter for the distribution of genetic variance.
#' @param ve prior value of residual variance.
#' @param dfve the number of degrees of freedom for the distribution of residual variance.
#' @param s2ve scale parameter for the distribution of residual variance.
#' @param outfreq frequency of information output on console, the default is 100.
#' @param seed seed for random sample.
#' @param verbose whether to print the iteration information.

#' @examples
#' set.seed(123)
#' n = 1000                        # number of observations
#' m = 100                         # number of SNPs
#' k = 10                          # number of causal SNPs
#' h2 = 0.5                        # heritability
#' Z = matrix(sample(c(0, 1, 2), n*m, prob=c(0.25, 0.5, 0.25), replace=TRUE), n)
#' Z = apply(Z, 2, function(x){x-mean(x)})
#' qtl = sort(sample(1:m, k))
#' betaTrue = array(0,m)
#' betaTrue[qtl] = rnorm(k)
#' g = Z%*%betaTrue
#' vg = var(g)
#' ve = (1-h2)/h2 * vg
#' y = g + rnorm(n,0,sqrt(ve))
#' 
#' ## fit model
#' fit = bayes(y=y, Z=Z, model="BayesCpi", niter=1000, nburn=500)
#' 
#' ## check results
#' cor(fit$g, betaTrue)
#' plot(fit$g, betaTrue)
#' 
#' #####load plink binary file#####
#' # bfile_path = system.file("extdata", "example", package = "hibayes")
#' # data = read_plink(bfile_path, out=tempfile())
#' # pheno = data$pheno
#' # geno = data$geno
#' # map = data$map
#' ## For GS/GP
#' # fit = bayes(y=pheno[,1], Z=geno, model="BayesR", niter=20000, nburn=10000)
#' ## For GWAS
#' # fit = bayes(y=pheno[,1], Z=geno, map=map, windsize=1e6, model="BayesR", niter=20000, nburn=10000)

#' @export

bayes <- 
function(
    y,
	Z,
    X = NULL,
	R = NULL,
    model = c("BayesB", "BayesA", "BayesLASSO", "BayesRR", "BayesBpi", "BayesC", "BayesCpi", "BayesR"),
    map = NULL,
    pi = NULL,
    fold = NULL,
    niter = 20000,
    nburn = 12000,
    windsize = NULL,
    wppa = 0.01,
    vg = NULL,
    dfvg = NULL,
    s2vg = NULL,
    ve = NULL,
    dfve = NULL,
    s2ve = NULL,
    outfreq = 100,
    seed = 666666,
    verbose = TRUE
){
	if(!is.null(windsize)){
		if(is.null(map)){
			stop("map information must be provided.")
		}else{
			if(ncol(map) != 3)	stop(" Only 3 columns are required in map.")
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
		windchr <- unique(chr)[match(tapply(map[, 1], windindx, unique), unique(sort(map[,1])))]
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
	if(!is.numeric(y)){y <- as.matrix(y)[,1,drop=TRUE]}
	if(!is.matrix(Z)){Z <- as.matrix(Z); gc()}
	if(!is.null(X)){
		if(!is.matrix(X))	X <- as.matrix(X)
		X_is_num <- apply(X, 2, is.numeric)
		if(!all(X_is_num))	stop("covariates must be a numeric matrix, please use 'model.matrix' to convert.")
	}
	if(!is.null(R)){
		if(!is.matrix(R))	R <- as.matrix(R)
		R <- apply(R, 2, as.character)
	}
	switch(
		match.arg(model), 
		"BayesRR"={
			res = BayesRR(y=y, X=Z, C=X, R=R, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, verbose=verbose)
		},
		"BayesA"={
			res = BayesA(y=y, X=Z, C=X, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, verbose=verbose)
		},
		"BayesLASSO"={
			res = BayesLASSO(y=y, X=Z, C=X, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, verbose=verbose)
		},
		"BayesB"={
			res = BayesB(y=y, X=Z, C=X, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, verbose=verbose)
		},
		"BayesBpi"={
			res = BayesBpi(y=y, X=Z, C=X, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, verbose=verbose)
		},
		"BayesC"={
			res = BayesC(y=y, X=Z, C=X, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, verbose=verbose)
		},
		"BayesCpi"={
			res = BayesCpi(y=y, X=Z, C=X, pi=pi, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, verbose=verbose)
		},
		"BayesR"={
			res = BayesR(y=y, X=Z, C=X, pi=pi, fold=fold, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, verbose=verbose)
		}
	)

	if(!is.null(windsize)){
		WPPA <- res$wppa
		WGVE <- res$wgve
		res <- head(res, length(res) - 2)
		res$gwas <- data.frame(windinfo, WPPA, WGVE)
	}
	return(res)
}
