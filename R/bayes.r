#' Bayes model
#'
#' Bayes linear regression model using individual level data
#' \deqn{y = X \beta + R r + M \alpha + e}
#' where \eqn{\beta} is a vector of estimated coefficient for covariates, and \eqn{r} is a vector of environmental random effects. \eqn{M} is a matrix of genotype covariate, \eqn{\alpha} is a vector of estimated marker effect size. \eqn{e} is a vector of residuals.
#' 
#' @param y vector of phenotype, use 'NA' for the missings. The number and order of individuals of y, M, X, R should be exactly the same.
#' @param M numeric matrix of genotype with individuals in rows and markers in columns, NAs are not allowed.
#' @param X (optional) covariate matrix of all individuals, all values should be in digits, characters are not allowed, please use 'model.matrix.lm' function to prepare it.
#' @param R (optional) environmental random effects matrix of all individuals, NAs are not allowed for the individuals with phenotypic value.
#' @param model bayes model including: "BayesB", "BayesA", "BayesL", "BayesRR", "BayesBpi", "BayesC", "BayesCpi", "BayesR", "BSLMM".
#' \itemize{
#' \item "BayesRR": Bayes Ridge Regression, all SNPs have non-zero effects and share the same variance, equals to RRBLUP or GBLUP. 
#' \item "BayesA": all SNPs have non-zero effects, and take different variance which follows an inverse chi-square distribution. 
#' \item "BayesB": only a small proportion of SNPs (1-Pi) have non-zero effects, and take different variance which follows an inverse chi-square distribution. 
#' \item "BayesBpi": the same with "BayesB", but 'Pi' is not fixed. 
#' \item "BayesC": only a small proportion of SNPs (1-Pi) have non-zero effects, and share the same variance. 
#' \item "BayesCpi": the same with "BayesC", but 'Pi' is not fixed. 
#' \item "BayesL": BayesLASSO, all SNPs have non-zero effects, and take different variance which follows an exponential distribution.
#' \item "BSLMM": all SNPs have non-zero effects, and take the same variance, but a small proportion of SNPs have additional shared variance. 
#' \item "BayesR": only a small proportion of SNPs have non-zero effects, and the SNPs are allocated into different groups, each group has the same variance. 
#' }
#' @param map (optional, only for GWAS) the map information of genotype, at least 3 columns are: SNPs, chromosome, physical position. 
#' @param Pi vector, the proportion of zero effect and non-zero effect SNPs, the first value must be the proportion of non-effect markers.
#' @param fold proportion of variance explained for groups of SNPs, the default is c(0, 0.0001, 0.001, 0.01).
#' @param niter the number of MCMC iteration.
#' @param nburn the number of iterations to be discarded.
#' @param windsize window size in bp for GWAS, the default is NULL.
#' @param windnum fixed number of SNPs in a window for GWAS, if it is specified, 'windsize' will be invalid, the default is NULL.
#' @param vg prior value of genetic variance.
#' @param dfvg the number of degrees of freedom for the distribution of genetic variance. 
#' @param s2vg scale parameter for the distribution of genetic variance.
#' @param ve prior value of residual variance.
#' @param dfve the number of degrees of freedom for the distribution of residual variance.
#' @param s2ve scale parameter for the distribution of residual variance.
#' @param lambda value of ridge regression for inverting a matrix.
#' @param outfreq frequency of collecting the estimated parameters and printing on console. Note that smaller frequency may have higher accuracy of estimated parameters, but would result in more time and memory for collecting process, on contrary, bigger frequency may have an negative effect on accuracy of estimations.
#' @param seed seed for random sample.
#' @param threads number of threads used for OpenMP.
#' @param verbose whether to print the iteration information on console.
#'
#' @references
#' Meuwissen, Theo HE, Ben J. Hayes, and Michael E. Goddard. "Prediction of total genetic value using genome-wide dense marker maps." Genetics 157.4 (2001): 1819-1829. \cr 
#' de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D., and Calus, M. P. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. Genetics, 193(2), 327-345. \cr 
#' Habier, David, et al. "Extension of the Bayesian alphabet for genomic selection." BMC bioinformatics 12.1 (2011): 1-12. \cr 
#' Yi, Nengjun, and Shizhong Xu. "Bayesian LASSO for quantitative trait loci mapping." Genetics 179.2 (2008): 1045-1055. \cr 
#' Zhou, Xiang, Peter Carbonetto, and Matthew Stephens. "Polygenic modeling with Bayesian sparse linear mixed models." PLoS genetics 9.2 (2013): e1003264. \cr 
#' Moser, Gerhard, et al. "Simultaneous discovery, estimation and prediction analysis of complex traits using a Bayesian mixture model." PLoS genetics 11.4 (2015): e1004969. \cr 
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$mu}{the regression intercept}
#' \item{$pi}{estimated proportion of zero effect and non-zero effect SNPs}
#' \item{$beta}{estimated coefficients for all covariates}
#' \item{$r}{estimated environmental random effects}
#' \item{$Vr}{estimated variance for all environmental random effect}
#' \item{$Vg}{estimated genetic variance}
#' \item{$Ve}{estimated residual variance}
#' \item{$h2}{estimated heritability (h2 = Vg / (Vr + Vg + Ve))}
#' \item{$alpha}{estimated effect size of all markers}
#' \item{$g}{genomic estimated breeding value}
#' \item{$e}{residuals of the model}
#' \item{$pip}{the frequency for markers to be included in the model during MCMC iteration, known as posterior inclusive probability (PIP)}
#' \item{$gwas}{WPPA is defined to be the window posterior probability of association, it is estimated by counting the number of MCMC samples in which \deqn{\alpha} is nonzero for at least one SNP in the window}
#' \item{$MCMCsamples}{the collected samples of posterior estimation for all the above parameters across MCMC iterations}
#' }
#'
#' @examples
#' # Load the example data attached in the package
#' pheno_file_path = system.file("extdata", "pheno.txt", package = "hibayes")
#' pheno = read.table(pheno_file_path, header=TRUE)
#' bfile_path = system.file("extdata", "geno", package = "hibayes")
#' data = read_plink(bfile_path, out=tempfile())
#' fam = data$fam
#' geno = data$geno
#' map = data$map
#' 
#' # Adjust the order of phenotype by genotype id
#' geno.id = fam[, 2]
#' pheno = pheno[match(geno.id, pheno[, 1]), ]
#' 
#' # Add fixed effects, covariates, and random effect
#' X <- model.matrix.lm(~as.numeric(scale)+as.factor(sex), data=pheno, na.action = "na.pass")
#' X <- X[, -1] #remove the intercept
#' # then fit the model as: fit = bayes(..., X=X, R=pheno[,c("group")], ...)
#' 
#' # For GS/GP
#' fit = bayes(y=pheno[, 2], M=geno, model="BayesR", niter=200, nburn=100)
#' \donttest{
#' # For GWAS
#' fit = bayes(y=pheno[, 2], M=geno, map=map, windsize=1e6, model="BayesCpi")
#' }
#' 
#' # The standard deviation of unknow parameters can be obtained from the list 'MCMCsamples':
#' snp_effect_sd = apply(fit$MCMCsamples$alpha, 1, sd)    # get the SD of estimated SNP effects for markers
#' gebv_pev = apply(fit$MCMCsamples$g, 1, var)    # get the prediction error variance (PEV) of estimated breeding values
#' 
#' @export

bayes <- 
function(
    y,
	M,
    X = NULL,
	R = NULL,
    model = c("BayesCpi", "BayesA", "BayesL", "BSLMM", "BayesR", "BayesB", "BayesC", "BayesBpi", "BayesRR"),
    map = NULL,
    Pi = NULL,
    fold = NULL,
    niter = 20000,
    nburn = 12000,
    windsize = NULL,
	windnum = NULL,
    vg = NULL,
    dfvg = NULL,
    s2vg = NULL,
    ve = NULL,
    dfve = NULL,
    s2ve = NULL,
	lambda = 0.0,
    outfreq = NULL,
    seed = 666666,
	threads = 4,
    verbose = TRUE
){
	set.seed(seed)
	model <- match.arg(model)
	if(!is.null(windsize) | !is.null(windnum)){
		if(model == "BayesA" || model == "BayesRR" || model == "BayesL")
			stop(paste0("can not implement GWAS analysis for the method: ", model))
		if(is.null(map)){
			stop("map information must be provided.")
		}else{
			if(ncol(map) < 3)	stop("At least 3 columns in map.")
		}
		if(any(is.na(map[,2]))){
			stop("NAs are not allowed in chromosome.")
		}
		if(any((map[,2]) == 0)){
			stop("0 is not allowed in chromosome.")
		}
		if(any(is.na(map[,3]))){
			stop("NAs are not allowed in physical position.")
		}
		if(any((map[,3]) == 0)){
			stop("0 is not allowed in physical position.")
		}
		map <- as.matrix(map[,c(2, 3)])
		chr <- map[, 1]
		suppressWarnings(pos_num <- as.numeric(map[, 2]))
		if(any(is.na(pos_num)))	stop("Characters are not allowed in physical position.")
		suppressWarnings(max.chr <- max(as.numeric(map[, 1]), na.rm=TRUE))
		if(is.infinite(max.chr))	max.chr <- 0
		suppressWarnings(map.xy.index <- which(!as.numeric(map[, 1]) %in% c(0 : max.chr)))
		if(length(map.xy.index) != 0){
			chr.xy <- unique(map[map.xy.index, 1])
			for(i in 1:length(chr.xy)){
				map[map[, 1] == chr.xy[i], 1] <- max.chr + i
			}
		}
		map <- apply(map, 2, as.numeric)
		# map <- matrix(as.numeric(map), nrow(map))
		chr <- chr[order(map[,1])]
		if(!is.null(windnum)){
			if(nrow(map) < windnum)	stop("Number of markers specified in a window is larger than the total number of markers.")
			windindx <- cutwind_by_num(map[,1], map[,2], windnum)
		}else{
			if(max(map[,2]) < windsize)	stop("Maximum of physical position is smaller than wind size.")
			windindx <- cutwind_by_bp(map[,1], map[,2], windsize)
		}
		windrange <- do.call(rbind, tapply(map[, 2], windindx, range))
		windsnpN <- tapply(map[, 2], windindx, length)
		windchr <- unique(chr)[match(tapply(map[, 1], windindx, unique), unique(sort(map[,1])))]
		windinfo <- data.frame(paste0("wind", 1:max(windindx)), windchr, windsnpN, windrange)
		colnames(windinfo) <- c("Wind", "Chr", "N", "Start", "End")
	}else{
		windindx <- NULL
	}
	if(is.null(outfreq) || outfreq <= 0){
		outfreq <- ifelse(niter > 1000, niter %/% 1000, 1)
	}
	if(outfreq >= (niter - nburn))	stop("bad setting for out frequency.")
	if(is.null(Pi)){
		if(model == "BayesR"){
			Pi <- c(0.95, 0.02, 0.02, 0.01)
			if(is.null(fold))	fold <- c(0, 0.0001, 0.001, 0.01)
		}else{
			Pi <- c(0.95, 0.05)
		}
	}
	if(!is.numeric(y)){
		y <- as.matrix(y)[,1,drop=TRUE]
		if(is.character(y))	stop("y is not a vector of digital values.")
	}
	yNA <- is.na(y)
	if(sum(yNA) != 0){
		if(verbose)	cat(sum(yNA), "'NA' have been detected in y\n")
		y <- y[!yNA]
	}
	if(!is.matrix(M)){M <- as.matrix(M); gc()}
	if(nrow(M) != length(yNA))	stop("number of individuals not match between 'y' and 'M'.")
	if(!is.null(X)){
		if(!is.matrix(X))	X <- as.matrix(X)
		if(nrow(X) != length(yNA))	stop("number of individuals not match between 'y' and 'X'.")
		X <- X[!yNA, , drop=FALSE]
		X_is_num <- apply(X, 2, is.numeric)
		if(!all(X_is_num))	stop("covariates must be a numeric matrix, please use 'model.matrix' to convert.")
		if(!all(apply(X, 2, function(x){unix <- unique(x); if(length(unix) == 1 && unix == 1){FALSE}else{TRUE}})))	stop("please remove intercept from covariates.")
		if(!all(apply(X, 2, function(x){length(unique(x)) > 1})))	stop("covariates should not be a constant.")
	}
	if(!is.null(R)){
		if(!is.matrix(R))	R <- as.matrix(R)
		if(nrow(R) != length(yNA))	stop("number of individuals not match between 'y' and 'R'.")
		R <- R[!yNA, , drop=FALSE]
		R <- apply(R, 2, as.character)
	}
	if(sum(yNA) != 0){
		Mp <- M[yNA, , drop=FALSE]
		M <- M[!yNA, , drop=FALSE]; gc()
	}else{
		Mp <- NULL;
	}
	if(model == "BSLMM"){
		G <- make_grm(M, lambda=lambda, inverse=TRUE, verbose=verbose)
		indx <- c(1:nrow(M))
		res <- BayesK(y=y, X=M, model=model, Pi=Pi, K=G, K_index=indx, fold=fold, C=X, R=R, niter=niter, nburn=nburn, windindx=windindx, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, threads=threads, verbose=verbose)
	}else{
		res <- Bayes(y=y, X=M, model=model, Pi=Pi, fold=fold, C=X, R=R, niter=niter, nburn=nburn, windindx=windindx, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, threads=threads, verbose=verbose)
	}
	res$MCMCsamples$g <- matrix(0, length(yNA), ncol(res$MCMCsamples$alpha))
	res$MCMCsamples$g[!yNA, ] <- M %*% res$MCMCsamples$alpha
	if(!is.null(Mp)){
		res$MCMCsamples$g[yNA, ] <- Mp %*% res$MCMCsamples$alpha
	}
	res$g <- apply(res$MCMCsamples$g, 1, mean)

	e <- rep(NA, length(yNA))
	e[!yNA] <- res$e
	res$e <- e

	if(!is.null(windsize) | !is.null(windnum)){
		WPPA <- res$gwas
		res$gwas <- data.frame(windinfo, WPPA)
	}

	return(res)
}
