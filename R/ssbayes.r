#' Single-step Bayes model
#'
#' Single-step Bayes linear regression model using individual level data and pedigree information
#' \deqn{y = X \beta + R r + M \alpha + U \epsilon + e}
#' where \eqn{y} is the vector of phenotypic values for both genotyped and non-genotyped individuals, \eqn{\beta} is a vector of estimated coefficient for covariates, \eqn{M} contains the genotype (\eqn{M_2}) for genotyped individuals and the imputed genotype (\eqn{M_1 = A_{12}A_{22}^{-1}M_2}) for non-genotyped individuals, \eqn{\epsilon} is the vector of genotype imputation error, \eqn{e} is a vector of residuals.
#'
#' @param formula a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a ‘~’ operator and the terms, separated by ‘+’ operators, on the right. Random-effects terms are distinguished by vertical bars (1|’) separating expressions for design matrices from grouping factors.
#' @param data the data frame containing the variables named in 'formula', NOTE that the first column in 'data' should be the individual id.
#' @param M numeric matrix of genotype with individuals in rows and markers in columns, NAs are not allowed.
#' @param M.id vector of id for genotype.
#' @param pedigree matrix of pedigree, 3 columns limited, the order of columns shoud be "id", "sir", "dam".
#' @param method bayes methods including: "BayesB", "BayesA", "BayesL", "BayesRR", "BayesBpi", "BayesC", "BayesCpi", "BayesR".
#' \itemize{
#' \item "BayesRR": Bayes Ridge Regression, all SNPs have non-zero effects and share the same variance, equals to RRBLUP or GBLUP. 
#' \item "BayesA": all SNPs have non-zero effects, and take different variance which follows an inverse chi-square distribution. 
#' \item "BayesB": only a small proportion of SNPs (1-Pi) have non-zero effects, and take different variance which follows an inverse chi-square distribution. 
#' \item "BayesBpi": the same with "BayesB", but 'Pi' is not fixed. 
#' \item "BayesC": only a small proportion of SNPs (1-Pi) have non-zero effects, and share the same variance. 
#' \item "BayesCpi": the same with "BayesC", but 'Pi' is not fixed. 
#' \item "BayesL": BayesLASSO, all SNPs have non-zero effects, and take different variance which follows an exponential distribution.
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
#' @param maf the effects of markers whose MAF is lower than the threshold will not be estimated.
#' @param dfvr the number of degrees of freedom for the distribution of environmental variance. 
#' @param s2vr scale parameter for the distribution of environmental variance.
#' @param vg prior value of genetic variance.
#' @param dfvg the number of degrees of freedom for the distribution of genetic variance. 
#' @param s2vg scale parameter for the distribution of genetic variance.
#' @param ve prior value of residual variance.
#' @param dfve the number of degrees of freedom for the distribution of residual variance.
#' @param s2ve scale parameter for the distribution of residual variance.
#' @param printfreq frequency of printing iterative details on console.
#' @param seed seed for random sample.
#' @param threads number of threads used for OpenMP.
#' @param verbose whether to print the iteration information on console.
#'
#' @return
#' the function returns a  a 'blrMod' object containing
#' \describe{
#' \item{$J}{coefficient for genotype imputation residuals}
#' \item{$Veps}{estimated variance of genotype imputation residuals}
#' \item{$epsilon}{genotype imputation residuals}
#' \item{$mu}{the regression intercept}
#' \item{$pi}{estimated proportion of zero effect and non-zero effect SNPs}
#' \item{$beta}{estimated coefficients for all covariates}
#' \item{$r}{estimated environmental random effects}
#' \item{$Vr}{estimated variance for all environmental random effect}
#' \item{$Vg}{estimated genetic variance}
#' \item{$Ve}{estimated residual variance}
#' \item{$h2}{estimated heritability (h2 = Vg / (Vr + Vg + Ve))}
#' \item{$g}{data.frame, the first column is the list of individual id, the second column is the genomic estimated breeding value for all individuals, including genotyped and non-genotyped.}
#' \item{$alpha}{estimated effect size of all markers}
#' \item{$e}{residuals of the model}
#' \item{$pip}{the frequency for markers to be included in the model during MCMC iteration, also known as posterior inclusive probability (PIP)}
#' \item{$gwas}{WPPA is defined to be the window posterior probability of association, it is estimated by counting the number of MCMC samples in which \deqn{\alpha} is nonzero for at least one SNP in the window}
#' \item{$MCMCsamples}{the collected samples of posterior estimation for all the above parameters across MCMC iterations}
#' }
#'
#' @references
#' Fernando, Rohan L., Jack CM Dekkers, and Dorian J. Garrick. "A class of Bayesian methods to combine large numbers of genotyped and non-genotyped animals for whole-genome analyses." Genetics Selection Evolution 46.1 (2014): 1-13. \cr 
#' Henderson, C.R.: A simple method for computing the inverse of a numerator relationship matrix used in prediction of breeding values. Biometrics 32(1), 69-83 (1976).
#'
#' @examples
#' # Load the example data attached in the package
#' pheno_file_path = system.file("extdata", "demo.phe", package = "hibayes")
#' pheno = read.table(pheno_file_path, header=TRUE)
#' 
#' bfile_path = system.file("extdata", "demo", package = "hibayes")
#' bin = read_plink(bfile_path, threads=1)
#' fam = bin$fam
#' geno = bin$geno
#' map = bin$map
#' 
#' pedigree_file_path = system.file("extdata", "demo.ped", package = "hibayes")
#' ped = read.table(pedigree_file_path, header=TRUE)
#' 
#' # For GS/GP
#' ## no environmental effects:
#' fit = ssbrm(T1~1, data=pheno, M=geno, M.id=fam[,2], pedigree=ped,
#' 	method="BayesCpi", niter=1000, nburn=600, thin=5, printfreq=100, threads=1)
#' 
#' ## overview of the returned results
#' summary(fit)
#' 
#' \donttest{
#' 
#' ## add fixed effects or covariates:
#' fit = ssbrm(T1~sex+bwt, data=pheno, M=geno, M.id=fam[,2], pedigree=ped,
#' 	method="BayesCpi")
#' 
#' ## add environmental random effects:
#' fit = ssbrm(T1~(1|loc)+(1|dam), data=pheno, M=geno, M.id=fam[,2],
#' 	pedigree=ped, method="BayesCpi")
#' 
#' # For GWAS
#' fit = ssbrm(T1~sex+bwt+(1|dam), data=pheno, M=geno, M.id=fam[,2],
#' 	pedigree=ped, method="BayesCpi", map=map, windsize=1e6)
#' }
#' 
#' # get the SD of estimated SNP effects for markers
#' summary(fit)$alpha
#' # get the SD of estimated breeding values
#' summary(fit)$g
#' 
#' @export

ssbrm <- 
function(
	formula,
	data = NULL,
	M = NULL,
	M.id = NULL,
	pedigree = NULL,
    method = c("BayesCpi", "BayesA", "BayesL", "BayesR", "BayesB", "BayesC", "BayesBpi", "BayesRR"),
	map = NULL,
    Pi = NULL,
    fold = NULL,
    niter = 20000,
    nburn = 12000,
	thin = 5,
    windsize = NULL,
	windnum = NULL,
	maf = 0.01,
	dfvr = NULL,
	s2vr = NULL,
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
	set.seed(seed)

	if(!inherits(formula, "formula"))	stop("not a standard formula.")
	if(is.null(data))	stop("no data assigned.")
	if(ncol(data) < 2)	stop("the first column in 'data' should be the individual id.")
	if(is.null(M))	stop("no genotype data.")
	if(is.null(M.id))	stop("please assign the individuals id to 'M.id'.")
	if(length(M.id) != nrow(M))	stop("number of individuals mismatched in 'M' and 'M.id'.")
	if(is.null(pedigree))	stop("pedigree should be provided for single-step bayesian model.")
	# if(nrow(data) != nrow(M))	stop("mismatched number of individuals between 'data' and 'M'.")

	method <- match.arg(method)
	if(!is.null(windsize) || !is.null(windnum)){
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

	yNA <- union(attr(model.frame(fixed_formula, data = data), "na.action"), attr(na.omit(R), "na.action"))
	yNA <- union(yNA, which(is.na(data[, as.character(formula[2]), drop = TRUE])))
	yNA <- c(1:nrow(data)) %in% yNA
	if(all(yNA))	stop("no effective data left.")
	if(verbose && sum(yNA))	cat(sum(yNA), " individuals have been removed due to missings.\n")
	
	if(!is.matrix(M)){M <- as.matrix(M); gc()}
	M.id <- as.character(as.matrix(M.id)[, 1, drop=TRUE])
	
	p <- apply(M, 2, function(x){p <- mean(x) / 2; return(min(c(p, 1 - p)))})
	if(sum(p < maf)){M[, p < maf] <- 0}
	if(ncol(pedigree) != 3)	stop("3 columns ('id', 'sir', 'dam') are required in pedigree.")
	ped <- as.matrix(pedigree)
	ped <- apply(ped, 2, as.character)
	ped.id <- unique(as.character(ped))
	Msub.id <- M.id[!M.id %in% ped.id]
	if(length(Msub.id) == length(M.id))	stop("no shared individuals between 'M.id' and 'pedigree'.")
	if(length(Msub.id)){
		ped <- rbind(ped, cbind(Msub.id, "0", "0"))
		ped.id <- c(Msub.id, ped.id)
	}
	if(all(ped.id %in% M.id))	stop("all individuals have been genotyped, no necessaries to fit single-step bayes model.")

	y.id <- as.character(data[!yNA, 1, drop=TRUE])
	ysub.id <- y.id[!y.id %in% ped.id]
	if(length(ysub.id) == length(y.id))	stop("no shared individuals between 'data' and 'pedigree'.")
	if(length(ysub.id)){
		if(verbose)	cat(length(ysub.id), " individuals cannot be found in genotype or pedigree\n")
		yNA[match(ysub.id, data[, as.character(formula[2]), drop = TRUE])] <- TRUE
		y.id <- as.character(data[!yNA, 1, drop=TRUE])
	}
	y <- data[!yNA, as.character(formula[2]), drop = TRUE]
	X <- model.matrix(fixed_formula, data = data[!yNA, , drop = FALSE])
	X <- X[, !apply(X, 2, function(x){all(x == 1)}), drop = FALSE]
	if(!ncol(X))	X <- NULL
	R <- R[!yNA, , drop = FALSE]

	pednew <- make_ped(ped[, 1], ped[, 2], ped[, 3], verbose)
	ped.id <- pednew[[1]]
	Ai <- make_Ainv(pednew[[2]], pednew[[3]], verbose)
	rm(ped, pednew); gc()
	g.indx <- match(M.id, ped.id)
	Mn.id <- ped.id[-g.indx]
	Ai.nn <- Ai[-g.indx, -g.indx]
	Ai.ng <- Ai[-g.indx,  g.indx]
	rm(Ai); gc();
	if(verbose)	cat("Linear solver for sparse matrix\n")
	A.ng <- solve(Ai.nn, -Ai.ng)
	rm(Ai.ng); gc();
	if(verbose)	cat("Start to impute genotype for", length(Mn.id), "individuals\n")
	# Mn <- geno_impute(A.ng, M, threads)
	Mn <- as.matrix(A.ng %*% M);
	J <- rep(-1, nrow(M))
	Jn <- as.vector(A.ng %*% J)
	rm(A.ng); gc();
	if(verbose)	cat("Impute genotype successfully\n")
	y.M.id <- M.id[M.id %in% y.id]
	y.Mn.id <- Mn.id[Mn.id %in% y.id]
	y.id.comb <- c(y.M.id, y.Mn.id)
	y.indx <- match(y.id.comb, y.id)
	y <- y[y.indx]
	if(!is.null(X)){X <- X[y.indx, , drop=FALSE]}
	if(!is.null(R)){R <- R[y.indx, , drop=FALSE]}
	y.Mn.indx <- match(y.Mn.id, Mn.id)
	y.M <- rbind(M[M.id %in% y.id, ], Mn[Mn.id %in% y.id, ])
	y.J <- c(J[M.id %in% y.id], Jn[Mn.id %in% y.id])

	res <- Bayes(y=y, X=y.M, model=method, Pi=Pi, fold=fold, C=X, R=R, epsl_y_J=y.J, epsl_Gi=Ai.nn, epsl_index=y.Mn.indx, niter=niter, thin=thin, nburn=nburn, windindx=windindx, dfvr=dfvr, s2vr=s2vr, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=printfreq, threads=threads, verbose=verbose)
	rm(y.M, y.J, Ai.nn); gc()

	if(length(y.Mn.indx)){
		res$MCMCsamples[["g"]] <- as.matrix(c(J, Jn)) %*% res$MCMCsamples[["J"]] + rbind(M %*% res$MCMCsamples$alpha, Mn %*% res$MCMCsamples$alpha + res$MCMCsamples$epsilon)
		epsilon <- data.frame(id = Mn.id, epsilon = res$epsilon)
		res$epsilon <- epsilon
	}else{
		warning("all phenotypic individuals have genotype information, thus can't fit imputation errors.")
		res$MCMCsamples[["g"]] <- rbind(M %*% res$MCMCsamples$alpha, Mn %*% res$MCMCsamples$alpha)
	}

	if(!is.null(res[["beta"]]))	names(res[["beta"]]) <- colnames(X)
	if(!is.null(res[["Vr"]]))	names(res[["Vr"]]) <- rand_term
	if(!is.null(res[["r"]]))	attr(res[["r"]], "nlevel") <- apply(R, 2, function(x){length(unique(x))})

	res[["g"]] <- data.frame(id = c(M.id, Mn.id), gebv = apply(res$MCMCsamples[["g"]], 1, mean))

	e <- rep(NA, length(y))
	e[match(y.id.comb, y.id)] <- res[["e"]]
	res[["e"]] <- data.frame(id = y.id, e = e)

	if(!is.null(windsize) | !is.null(windnum)){
		WPPA <- res$gwas
		res$gwas <- data.frame(windinfo, WPPA)
	}
	res$call <- paste(formula[2], paste(formula[3], "J", "M[pedigree]", sep=" + "), sep=' ~ ')
	attr(res$call, "model") <- paste0("Single-step Bayesian model fit by [", method, "]")
	class(res) <- "blrMod"
	return(res)
}
