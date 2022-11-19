#' Bayes model
#'
#' Bayes linear regression model using individual level data
#' \deqn{y = X \beta + R r + M \alpha + e}
#' where \eqn{\beta} is a vector of estimated coefficient for covariates, and \eqn{r} is a vector of environmental random effects. \eqn{M} is a matrix of genotype covariate, \eqn{\alpha} is a vector of estimated marker effect size. \eqn{e} is a vector of residuals.
#' 
#' @param formula a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a ‘~’ operator and the terms, separated by ‘+’ operators, on the right. Random-effects terms are distinguished by vertical bars (1|’) separating expressions for design matrices from grouping factors.
#' @param data the data frame containing the variables named in 'formula', NOTE that the first column in 'data' should be the individual id.
#' @param M numeric matrix of genotype with individuals in rows and markers in columns, NAs are not allowed.
#' @param M.id vector of id for genotyped individuals, NOTE that no need to adjust the order of id to be the same between 'data' and 'M', the package will do it automatically.
#' @param method bayes methods including: "BayesB", "BayesA", "BayesL", "BayesRR", "BayesBpi", "BayesC", "BayesCpi", "BayesR", "BSLMM".
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
#' @param thin the number of thinning after burn-in. Note that smaller thinning frequency may have higher accuracy of estimated parameters, but would result in more memory for collecting process, on contrary, bigger frequency may have negative effect on accuracy of estimations.
#' @param windsize window size in bp for GWAS, the default is NULL.
#' @param windnum fixed number of SNPs in a window for GWAS, if it is specified, 'windsize' will be invalid, the default is NULL.
#' @param dfvr the number of degrees of freedom for the distribution of environmental variance. 
#' @param s2vr scale parameter for the distribution of environmental variance.
#' @param vg prior value of genetic variance.
#' @param dfvg the number of degrees of freedom for the distribution of genetic variance. 
#' @param s2vg scale parameter for the distribution of genetic variance.
#' @param ve prior value of residual variance.
#' @param dfve the number of degrees of freedom for the distribution of residual variance.
#' @param s2ve scale parameter for the distribution of residual variance.
#' @param lambda value of ridge regression for inverting a matrix.
#' @param printfreq frequency of printing iterative details on console.
#' @param seed seed for random sample.
#' @param threads number of threads used for OpenMP.
#' @param verbose whether to print the iteration information on console.
#'
#' @details
#' \itemize{
#' 	   \item{the fixed effects and covariates in 'formula' must be in factors and numeric, respectively. if not, please remember to use 'as.factor' and 'as.numeric' to transform.}
#'     \item{the package has the automatical function of taking the intersection and adjusting the order of id between 'data' and the genotype 'M', thus the first column in 'data' should be the individual id.}
#'     \item{if any one of the options 'windsize' and 'windnum' is specified, the GWAS results will be returned, and the 'map' information must be provided, in which the physical positions should be all in digital values.}
#'     \item{the 'windsize' or 'windnum' option only works for the methods of which the assumption has a proportion of zero effect markers, e.g., BayesB, BayesBpi, BayesC, BayesCpi, BSLMM, and BayesR.}
#' }
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
#' the function returns a 'blrMod' object containing
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
#' pheno_file_path = system.file("extdata", "demo.phe", package = "hibayes")
#' pheno = read.table(pheno_file_path, header=TRUE)
#' 
#' bfile_path = system.file("extdata", "demo", package = "hibayes")
#' bin = read_plink(bfile_path)
#' fam = bin$fam
#' geno = bin$geno
#' map = bin$map
#' 
#' # For GS/GP
#' ## no environmental effects:
#' fit = bayes(T1~1, data=pheno, M=geno, M.id=fam[,2], method="BayesCpi",
#' 	niter=2000, nburn=1200, thin=5)
#' 
#' ## overview of the returned results
#' summary(fit)
#' 
#' \donttest{
#'
#' ## add fixed effects or covariates:
#' fit = bayes(T1~sex+season+day+bwt, data=pheno, M=geno, M.id=fam[,2],
#' 	method="BayesCpi")
#'  
#' ## add environmental random effects:
#' fit = bayes(T1~sex+(1|loc)+(1|dam), data=pheno, M=geno, M.id=fam[,2],
#' 	method="BayesCpi")
#' 
#' # For GWAS
#' fit = bayes(T1~sex+bwt+(1|dam), data=pheno, M=geno, M.id=fam[,2],
#' 	method="BayesCpi", map=map, windsize=1e6)
#' }
#' 
#' # get the SD of estimated SNP effects for markers
#' summary(fit)$alpha
#' # get the SD of estimated breeding values
#' summary(fit)$g
#' 
#' @export

bayes <- 
function(
	formula,
	data = NULL,
	M = NULL,
	M.id = NULL,
	method = c("BayesCpi", "BayesA", "BayesL", "BSLMM", "BayesR", "BayesB", "BayesC", "BayesBpi", "BayesRR"),
    map = NULL,
    Pi = NULL,
    fold = NULL,
    niter = 20000,
    nburn = 12000,
	thin = 5,
    windsize = NULL,
	windnum = NULL,
	dfvr = NULL,
	s2vr = NULL,
    vg = NULL,
    dfvg = NULL,
    s2vg = NULL,
    ve = NULL,
    dfve = NULL,
    s2ve = NULL,
	lambda = 0.0,
    printfreq = 100,
    seed = 666666,
	threads = 4,
    verbose = TRUE
){

	set.seed(seed)

	if(!inherits(formula, "formula"))	stop("not a standard formula.")
	if(is.null(data))	stop("no data assigned.")
	if(ncol(data) < 2)	stop("the first column in 'data' should be the individual id.")
	if(is.null(M))	stop("no genotype data.")
	if(is.null(M.id))	stop("please assign the individuals id to 'M.id'.")
	if(length(M.id) != nrow(M))	stop("number of individuals mismatched in 'M' and 'M.id'.")
	# if(nrow(data) != nrow(M))	stop("mismatched number of individuals between 'data' and 'M'.")

	M.id <- as.character(as.matrix(M.id)[, 1, drop=TRUE])
	if(length(intersect(as.character(data[, 1, drop=TRUE]), M.id)) == 0)
		stop("no shared individuals between 'M.id' and the first column in 'data'.")
	
	data <- data[match(M.id, as.character(data[, 1, drop=TRUE])), ]

	myformula  <- paste(formula[2], formula[3], sep='~')
	rand_term <- unlist(str_extract_all(myformula, "(?<=(\\+ |~)\\(1 \\| )[:\\w\\d]+(?=\\))"))
	R <- NULL
	for(r in rand_term){
		split_str = unlist(strsplit(r,":"))
		if(length(split_str) != 1){
			Ri <- as.matrix(apply(data[, split_str], 1, function(x){if(sum(is.na(x)) != 0){NA}else{paste(x, collapse = ":")}}))
		}else{
			Ri <- as.matrix(as.character(data[, split_str]))
		}
		R <- cbind(R, Ri)
	}

	fixed_formula  <- str_replace_all(myformula, pattern = "(\\+ |(?<=~))\\(1 \\| [:\\w\\d]+\\)", replacement = "" )

	fixed_formula  <- str_replace_all(fixed_formula, pattern = "~ *\\+ ", replacement = "~" )
	fixed_formula  <- str_replace_all(fixed_formula, pattern = "~ *\\- ", replacement = "~" )
	fixed_formula  <- str_replace_all(fixed_formula, pattern = "~ *$", replacement = "~1" )

	# warning
	warn_pattern = "(. |~)\\(.*? \\| .*?\\)"
	if (str_detect(fixed_formula, pattern = warn_pattern)) {
		stop(paste0("Invalid random effects expression '", 
			paste(unlist(str_extract_all(fixed_formula, warn_pattern)), collapse = "', '"), 
			"',\n  it should be in the format '(1 | x)' or '+ (1 | x1:x2:...:xn)'."))
	}

	fixed_formula <- formula(fixed_formula)

	# print(rand_term)
	# print(fixed_formula)

	yNA <- union(attr(model.frame(fixed_formula, data = data), "na.action"), attr(na.omit(R), "na.action"))
	yNA <- union(yNA, which(is.na(data[, as.character(formula[2]), drop = TRUE])))
	yNA <- c(1:length(M.id)) %in% yNA
	if(all(yNA))	stop("no effective data left.")
	
	X <- model.matrix(fixed_formula, data = data[!yNA, , drop = FALSE])
	X <- X[, !apply(X, 2, function(x){all(x == 1)}), drop = FALSE]
	if(!ncol(X))	X <- NULL
	R <- R[!yNA, , drop = FALSE]
	
	# print(head(X))
	# print(head(R))

	method <- match.arg(method)
	if(!is.null(windsize) | !is.null(windnum)){
		if(method == "BayesA" || method == "BayesRR" || method == "BayesL")
			stop(paste0("can not implement GWAS analysis for the method: ", method))
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
	if(thin >= (niter - nburn))	stop("bad setting for collecting frequency 'thin'.")
	if(printfreq <= 0)	verbose <- FALSE
	if(is.null(Pi)){
		if(method == "BayesR"){
			Pi <- c(0.95, 0.02, 0.02, 0.01)
			if(is.null(fold))	fold <- c(0, 0.0001, 0.001, 0.01)
		}else{
			Pi <- c(0.95, 0.05)
		}
	}
	
	y <- data[!yNA, as.character(formula[2]), drop = TRUE]
	if(is.character(y))	stop(paste0("'", as.character(formula[2]), "'", " is not a vector of digital values."))

	if(!is.matrix(M)){M <- as.matrix(M); gc()}

	if(sum(yNA) != 0){
		Mp <- M[yNA, , drop=FALSE]
		M <- M[!yNA, , drop=FALSE]; gc()
	}else{
		Mp <- NULL;
	}
	if(method == "BSLMM"){
		G <- make_grm(M, lambda=lambda, inverse=TRUE, verbose=verbose)
		indx <- c(1:nrow(M))
		res <- BayesK(y=y, X=M, model=method, Pi=Pi, K=G, K_index=indx, fold=fold, C=X, R=R, niter=niter, nburn=nburn, thin=thin, windindx=windindx, dfvr=dfvr, s2vr=s2vr, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=printfreq, threads=threads, verbose=verbose)
	}else{
		res <- Bayes(y=y, X=M, model=method, Pi=Pi, fold=fold, C=X, R=R, niter=niter, nburn=nburn, thin=thin, windindx=windindx, dfvr=dfvr, s2vr=s2vr, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=printfreq, threads=threads, verbose=verbose)
	}

	if(!is.null(res[["beta"]]))	names(res[["beta"]]) <- colnames(X)
	if(!is.null(res[["Vr"]]))	names(res[["Vr"]]) <- rand_term
	if(!is.null(res[["r"]]))	attr(res[["r"]], "nlevel") <- apply(R, 2, function(x){length(unique(x))})

	res$MCMCsamples[["g"]] <- matrix(0, length(M.id), ncol(res$MCMCsamples$alpha))
	res$MCMCsamples[["g"]][!yNA, ] <- M %*% res$MCMCsamples$alpha
	if(!is.null(Mp)){
		res$MCMCsamples[["g"]][yNA, ] <- Mp %*% res$MCMCsamples$alpha
	}
	res[["g"]] <- data.frame(id = M.id, gebv = apply(res$MCMCsamples[["g"]], 1, mean))

	res[["e"]] <- data.frame(id = M.id[!yNA], e = res[["e"]])

	if(!is.null(windsize) | !is.null(windnum)){
		WPPA <- res$gwas
		res$gwas <- data.frame(windinfo, WPPA)
	}
	res$call <- paste(formula[2], paste(formula[3], "M", sep=" + "), sep=' ~ ')
	attr(res$call, "model") <- paste0("Individual level Bayesian model fit by [", method, "]")
	class(res) <- "blrMod"
	return(res)
}
