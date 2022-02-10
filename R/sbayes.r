#' SBayes model
#'
#' Bayes linear regression model using summary level data
#'
#' @param sumstat matrix of summary data, details refer to https://cnsgenomics.com/software/gcta/#COJO.
#' @param ldm dense or sparse matrix, ld for reference panel (m * m, m is the number of SNPs). NOTE that the order of SNPs should be consistent with summary data.
#' @param model bayes model including: "BayesB", "BayesA", "BayesL", "BayesRR", "BayesBpi", "BayesC", "BayesCpi", "BayesR", "CG".
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
#' @param windsize window size in bp for GWAS, the default is 1e6.
#' @param windnum fixed number of SNPs in a window for GWAS, if it is specified, 'windsize' will be invalid, the default is NULL.
#' @param vg prior value of genetic variance.
#' @param dfvg the number of degrees of freedom for the distribution of genetic variance. 
#' @param s2vg scale parameter for the distribution of genetic variance.
#' @param ve prior value of residual variance.
#' @param dfve the number of degrees of freedom for the distribution of residual variance.
#' @param s2ve scale parameter for the distribution of residual variance.
#' @param outfreq frequency of information output on console, the default is 100.
#' @param seed seed for random sample.
#' @param threads number of threads used for OpenMP.
#' @param verbose whether to print the iteration information.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$pi}{estimated proportion of zero effect and non-zero effect SNPs}
#' \item{$vg}{estimated genetic variance}
#' \item{$ve}{estimated residual variance}
#' \item{$alpha}{estimated effect size of all markers}
#' \item{$pip}{the frequency for markers to be included in the model during MCMC iteration, also known as posterior inclusive probability (PIP)}
#' \item{$gwas}{WPPA is defined to be the window posterior probability of association, it is estimated by counting the number of MCMC samples in which \deqn{\alpha} is nonzero for at least one SNP in the window}
#' }
#' 
#' @references
#' Lloyd-Jones, Luke R., et al. "Improved polygenic prediction by Bayesian multiple regression on summary statistics." Nature communications 10.1 (2019): 1-11.
#' 
#' @examples
#' bfile_path = system.file("extdata", "geno", package = "hibayes")
#' data = read_plink(bfile_path, out=tempfile())
#' geno = data$geno
#' map = data$map
#' head(map)
#' sumstat_path = system.file("extdata", "geno.ma", package = "hibayes")
#' sumstat = read.table(sumstat_path, header=TRUE)
#' head(sumstat)
#' 
#' \donttest{
#' # computate ld variance covariance matrix
#' ldm1 = ldmat(geno, threads=4)   #chromosome wide full ld matrix
#' 
#' # if the order of SNPs in genotype is not consistent with the order in sumstat file, 
#' # prior adjusting is necessary.
#' indx = match(map[, 1], sumstat[, 1])
#' sumstat = sumstat[indx, ]
#' 
#' # fit model
#' fit = sbayes(sumstat=sumstat, ldm=ldm1, model="BayesR")
#' }
#'
#' @export

sbayes <- 
function(
    sumstat,
    ldm,
    model = c("BayesB", "BayesA", "BayesL", "BayesRR", "BayesBpi", "BayesC", "BayesCpi", "BayesR", "CG"),
    map = NULL,
    Pi = NULL,
    lambda = NULL,
    fold = NULL,
    niter = 20000,
    nburn = 14000,
    windsize = NULL,
	windnum = NULL,
    vg = NULL,
    dfvg = NULL,
    s2vg = NULL,
    ve = NULL,
    dfve = NULL,
    s2ve = NULL,
    outfreq = 100,
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
	model <- match.arg(model)
	if(!is.null(windsize) || !is.null(windnum)){
		if(model == "BayesA" || model == "BayesRR" || model == "BayesL")
			stop(paste0("can not implement GWAS analysis for the method: ", model))
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
	if(is.null(Pi)){
		if(match.arg(model) == "BayesR"){
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
	if(model != "CG"){
		if(sparse){
			res <- SBayesS(sumstat=sumstat, ldm=ldm, model=model, Pi=Pi, fold=fold, niter=niter, nburn=nburn, windindx=windindx, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, threads=threads, verbose=verbose)
		}else{
			res <- SBayesD(sumstat=sumstat, ldm=ldm, model=model, Pi=Pi, fold=fold, niter=niter, nburn=nburn, windindx=windindx, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, threads=threads, verbose=verbose)
		}
	}else{
		if(!is.null(lambda)){
			if(length(lambda) == 1){
				lambda = rep(lambda, nrow(sumstat))
			}else if(length(lambda) != nrow(sumstat)){
				stop("length of lambda should equal to the length of SNPs.")
			}
		}
		if(sparse){
			res <- conjgt_spa(sumstat=sumstat, ldm=ldm, lambda = lambda, outfreq=outfreq, verbose=verbose)
		}else{
			res <- conjgt_den(sumstat=sumstat, ldm=ldm, lambda = lambda, outfreq=outfreq, verbose=verbose)
		}
	}
	if((!is.null(windsize) | !is.null(windnum)) & match.arg(model) != "CG"){
		WPPA <- res$gwas
		res$gwas <- data.frame(windinfo, WPPA)
	}

	return(res)
}
