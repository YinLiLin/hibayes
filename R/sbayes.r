#' SBayes model
#'
#' Bayes linear regression model using summary level data
#'
#' @param sumstat matrix of summary data, details refer to https://cnsgenomics.com/software/gcta/#COJO.
#' @param ldm dense or sparse matrix, ld for reference panel (m * m, m is the number of SNPs). NOTE that the order of SNPs should be consistent with summary data.
#' @param method bayes methods including: "BayesB", "BayesA", "BayesL", "BayesRR", "BayesBpi", "BayesC", "BayesCpi", "BayesR", "CG".
#' \itemize{
#' \item "BayesRR": Bayes Ridge Regression, all SNPs have non-zero effects and share the same variance, equals to RRBLUP or GBLUP. 
#' \item "BayesA": all SNPs have non-zero effects, and take different variance which follows an inverse chi-square distribution. 
#' \item "BayesB": only a small proportion of SNPs (1-Pi) have non-zero effects, and take different variance which follows an inverse chi-square distribution. 
#' \item "BayesBpi": the same with "BayesB", but 'Pi' is not fixed. 
#' \item "BayesC": only a small proportion of SNPs (1-Pi) have non-zero effects, and share the same variance. 
#' \item "BayesCpi": the same with "BayesC", but 'Pi' is not fixed. 
#' \item "BayesL": BayesLASSO, all SNPs have non-zero effects, and take different variance which follows an exponential distribution. 
#' \item "BayesR": only a small proportion of SNPs have non-zero effects, and the SNPs are allocated into different groups, each group has the same variance. 
#' \item "CG": conjugate gradient algorithm with assigned lambda. 
#' }
#' @param map (optional, only for GWAS) the map information of genotype, at least 3 columns are: SNPs, chromosome, physical position. 
#' @param Pi vector, the proportion of zero effect and non-zero effect SNPs, the first value must be the proportion of non-effect markers.
#' @param lambda value or vector, the ridge regression value for each SNPs.
#' @param fold percentage of variance explained for groups of SNPs, the default is c(0, 0.0001, 0.001, 0.01).
#' @param niter the number of MCMC iteration.
#' @param nburn the number of iterations to be discarded.
#' @param thin the number of thinning after burn-in. Note that smaller thinning frequency may have higher accuracy of estimated parameters, but would result in more memory for collecting process, on contrary, bigger frequency may have negative effect on accuracy of estimations.
#' @param windsize window size in bp for GWAS, the default is 1e6.
#' @param windnum fixed number of SNPs in a window for GWAS, if it is specified, 'windsize' will be invalid, the default is NULL.
#' @param vg prior value of genetic variance.
#' @param dfvg the number of degrees of freedom for the distribution of genetic variance. 
#' @param s2vg scale parameter for the distribution of genetic variance.
#' @param ve prior value of residual variance.
#' @param dfve the number of degrees of freedom for the distribution of residual variance.
#' @param s2ve scale parameter for the distribution of residual variance.
#' @param printfreq frequency of collecting the estimated parameters and printing on console. Note that smaller frequency may have higher accuracy of estimated parameters, but would result in more time and memory for collecting process, on contrary, bigger frequency may have an negative effect on accuracy of estimations.
#' @param seed seed for random sample.
#' @param threads number of threads used for OpenMP.
#' @param verbose whether to print the iteration information on console.
#'
#' @details
#' \itemize{
#' 	   \item{if any one of the options 'windsize' and 'windnum' is specified, the GWAS results will be returned, and the 'map' information must be provided, in which the physical positions should be all in digital values.}
#'     \item{the 'windsize' or 'windnum' option only works for the methods of which the assumption has a proportion of zero effect markers, e.g., BayesB, BayesBpi, BayesC, BayesCpi, BSLMM, and BayesR.}
#' }
#'
#' @return
#' the function returns a 'blrMod' object containing
#' \describe{
#' \item{$pi}{estimated proportion of zero effect and non-zero effect SNPs}
#' \item{$Vg}{estimated genetic variance}
#' \item{$Ve}{estimated residual variance}
#' \item{$h2}{estimated heritability (h2 = Vg / (Vg + Ve))}
#' \item{$alpha}{estimated effect size of all markers}
#' \item{$pip}{the frequency for markers to be included in the model during MCMC iteration, also known as posterior inclusive probability (PIP)}
#' \item{$gwas}{WPPA is defined to be the window posterior probability of association, it is estimated by counting the number of MCMC samples in which \deqn{\alpha} is nonzero for at least one SNP in the window}
#' \item{$MCMCsamples}{the collected samples of posterior estimation for all the above parameters across MCMC iterations}
#' }
#' 
#' @references
#' Lloyd-Jones, Luke R., et al. "Improved polygenic prediction by Bayesian multiple regression on summary statistics." Nature communications 10.1 (2019): 1-11.
#' 
#' @examples
#' bfile_path = system.file("extdata", "demo", package = "hibayes")
#' bin = read_plink(bfile_path, threads=1)
#' fam = bin$fam
#' geno = bin$geno
#' map = bin$map
#' 
#' sumstat_path = system.file("extdata", "demo.ma", package = "hibayes")
#' sumstat = read.table(sumstat_path, header=TRUE)
#' head(sumstat)
#' 
#' \donttest{
#' # computate ld variance covariance matrix
#' ## construct genome wide full variance-covariance matrix
#' ldm1 <- ldmat(geno, threads=4)	
#' ## construct genome wide sparse variance-covariance matrix
#' # ldm2 <- ldmat(geno, chisq=5, threads=4)	
#' ## construct chromosome wide full variance-covariance matrix
#' # ldm3 <- ldmat(geno, map, ldchr=FALSE, threads=4)	
#' ## construct chromosome wide sparse variance-covariance matrix
#' # ldm4 <- ldmat(geno, map, ldchr=FALSE, chisq=5, threads=4)
#' 
#' # if the order of SNPs in genotype is not consistent with the order in sumstat file, 
#' # prior adjusting is necessary.
#' indx = match(map[, 1], sumstat[, 1])
#' sumstat = sumstat[indx, ]
#' 
#' # fit model
#' fit = sbrm(sumstat=sumstat, ldm=ldm1, method="BayesCpi", Pi = c(0.95, 0.05), 
#' 	niter=20000, nburn=12000, seed=666666, map=map, windsize=1e6, threads=1)
#' 
#' # overview of the returned results
#' summary(fit)
#' 
#' # get the SD of estimated SNP effects for markers
#' summary(fit)$alpha
#' }
#'
#' @export

sbrm <- 
function(
    sumstat,
    ldm,
    method = c("BayesB", "BayesA", "BayesL", "BayesRR", "BayesBpi", "BayesC", "BayesCpi", "BayesR", "CG"),
    map = NULL,
    Pi = NULL,
    lambda = NULL,
    fold = NULL,
    niter = NULL,
    nburn = NULL,
	thin = 5,
    windsize = NULL,
	windnum = NULL,
    vg = NULL,
    dfvg = NULL,
    s2vg = NULL,
    ve = NULL,
    dfve = NULL,
    s2ve = NULL,
    printfreq = 100,
    seed = 666666,
	threads = 4,
    verbose = TRUE
){
	if(is.matrix(ldm)){
		sparse = FALSE
	}else if(is(ldm, "dgCMatrix")){
		sparse = TRUE
	}else{
		stop("Unrecognized type of ldm.")
	}
	set.seed(seed)
	method <- match.arg(method)
	if(!is.null(windsize) || !is.null(windnum)){
		if(method == "BayesA" || method == "BayesRR" || method == "BayesL")
			stop(paste0("can not implement GWAS analysis for the method: ", method))
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
		colnames(windinfo) <- c("WIND", "CHR", "NUM", "START", "END")
	}else{
		windindx <- NULL
	}
	if(is.null(niter)){
		niter <- ifelse(method == "BayesR", 50000, 20000)
	}
	if(is.null(nburn)){
		nburn <- ifelse(method == "BayesR", 30000, 12000)
	}
	if(thin >= (niter - nburn))	stop("bad setting for collecting frequency 'thin'.")
	if(printfreq <= 0)	verbose <- FALSE
	if(is.null(Pi)){
		if(match.arg(method) == "BayesR"){
			Pi <- c(0.95, 0.02, 0.02, 0.01)
			if(is.null(fold))	fold <- c(0, 0.0001, 0.001, 0.01)
		}else{
			Pi <- c(0.95, 0.05)
		}
	}
	# if(ncol(sumstat) != 8)	stop("Inappropriate summary data format.")
	# if(any(sumstat[, 4] <=0 | sumstat[, 4] >= 1))	stop("Frequency of allele should be at range of (0, 1), please remove the SNPs whose allele frequency are out of this range.")
	# if(any(is.na(sumstat[, c(4, 5, 6)])))	stop("'NA' is not allowed for BETA and SE.")
	sumstat <- sumstat[, c(4, 5, 6, 8)]
	sumstat <- data.matrix(sumstat)
	if(method != "CG"){
		if(sparse){
			res <- SBayesS(sumstat=sumstat, ldm=ldm, model=method, Pi=Pi, fold=fold, niter=niter, nburn=nburn, windindx=windindx, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=printfreq, threads=threads, verbose=verbose)
		}else{
			res <- SBayesD(sumstat=sumstat, ldm=ldm, model=method, Pi=Pi, fold=fold, niter=niter, nburn=nburn, windindx=windindx, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=printfreq, threads=threads, verbose=verbose)
		}
	}else{
		if(!is.null(lambda)){
			if(length(lambda) == 1){
				lambda = rep(lambda, nrow(sumstat))
			}else if(length(lambda) != nrow(sumstat)){
				stop("length of lambda should be equal to the number of SNPs.")
			}
		}
		if(sparse){
			res <- conjgt_spa(sumstat=sumstat, ldm=ldm, lambda = lambda, outfreq=printfreq, verbose=verbose)
		}else{
			res <- conjgt_den(sumstat=sumstat, ldm=ldm, lambda = lambda, outfreq=printfreq, verbose=verbose)
		}
	}
	if((!is.null(windsize) | !is.null(windnum)) & match.arg(method) != "CG"){
		WPPA <- res$gwas
		res$gwas <- data.frame(windinfo, WPPA)
	}
	res$call <- paste0("b ~ nD","\U207b\U00b9", "V", "\U03b1", " + e")
	attr(res$call, "model") <- paste0("Summary level Bayesian model fit by [", method, "]")
	class(res) <- "blrMod"
	return(res)
}
