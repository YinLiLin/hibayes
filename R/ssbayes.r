#' Single-step Bayes model
#'
#' Single-step Bayes linear regression model using individual level data and pedigree information
#' \deqn{y = X \beta + R r + M \alpha + U \epsilon + e}
#' where \eqn{y} is the vector of phenotypic values for both genotyped and non-genotyped individuals, \eqn{\beta} is a vector of estimated coefficient for covariates, \eqn{M} contains the genotype (\eqn{M_2}) for genotyped individuals and the imputed genotype (\eqn{M_1 = A_{12}A_{22}^{-1}M_2}) for non-genotyped individuals, \eqn{\epsilon} is the vector of genotype imputation error, \eqn{e} is a vector of residuals.
#'
#' @param y vector of phenotype, use 'NA' for the missings.
#' @param y.id vector of id for phenotype.
#' @param M numeric matrix of genotype with individuals in rows and markers in columns, NAs are not allowed.
#' @param M.id vector of id for genotype.
#' @param P matrix of pedigree, 3 columns limited, the order of columns shoud be "id", "sir", "dam".
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
#' @param outfreq frequency of collecting the estimated parameters and printing on console. Note that smaller frequency may have higher accuracy of estimated parameters, but would result in more time and memory for collecting process, on contrary, bigger frequency may have an negative effect on accuracy of estimations.
#' @param seed seed for random sample.
#' @param threads number of threads used for OpenMP.
#' @param verbose whether to print the iteration information on console.
#'
#' @return
#' the function returns a list containing
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
#' pheno_file_path = system.file("extdata", "pheno.txt", package = "hibayes")
#' pheno = read.table(pheno_file_path, header=TRUE)
#' pedigree_file_path = system.file("extdata", "ped.txt", package = "hibayes")
#' ped = read.table(pedigree_file_path, header=TRUE)
#' bfile_path = system.file("extdata", "geno", package = "hibayes")
#' data = read_plink(bfile_path, out=tempfile())
#' fam = data$fam
#' geno = data$geno
#' map = data$map
#' 
#' # NOTE: for ssbayes model, there is no NEED to adjust the order of id in different files
#' geno.id = fam[, 2]
#' pheno.id = pheno[, 1]
#' 
#' # Add fixed effects, covariates, and random effect
#' X <- model.matrix.lm(~as.numeric(scale)+as.factor(sex), data=pheno, na.action = "na.pass")
#' X <- X[, -1] #remove the intercept
#' # then fit the model as: fit = ssbayes(..., X=X, R=pheno[,c("group")], ...)
#' 
#' \donttest{
#' # For GS/GP
#' fit = ssbayes(y=pheno[, 2], y.id=pheno.id, M=geno, M.id=geno.id, P=ped, 
#' 			model="BayesR", niter=200, nburn=100)
#' # For GWAS
#' fit = ssbayes(y=pheno[, 2], y.id=pheno.id, M=geno, M.id=geno.id, P=ped, 
#' 			map=map, windsize=1e6, model="BayesCpi")
#' 
#' # The standard deviation of unknow parameters can be obtained from the list 'MCMCsamples':
#' snp_effect_sd = apply(fit$MCMCsamples$alpha, 1, sd)    # get the SD of estimated SNP effects for markers
#' gebv_pev = apply(fit$MCMCsamples$g, 1, var)    # get the prediction error variance (PEV) of estimated breeding values
#' }
#' 
#' @export

ssbayes <- 
function(
    y,
	y.id,
	M,
	M.id,
	P,
	X = NULL,
	R = NULL,
    model = c("BayesCpi", "BayesA", "BayesL", "BayesR", "BayesB", "BayesC", "BayesBpi", "BayesRR"),
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
    outfreq = NULL,
    seed = 666666,
	threads = 4,
    verbose = TRUE
){
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
	y.id <- as.character(as.matrix(y.id)[, 1, drop=TRUE])
	M.id <- as.character(as.matrix(M.id)[, 1, drop=TRUE])
	if(!is.numeric(y)){
		y <- as.matrix(y)[, 1, drop=TRUE]
		if(is.character(y))	stop("y is not a vector of digital values.")
	}
	if(length(y) != length(y.id))	stop("number of individuals not match between 'y' and 'y.id'.")
	e <- rep(NA, length(y))
	yNA <- is.na(y)
	ytmp.id <- y.id
	if(sum(yNA) != 0){
		if(verbose)	cat(sum(yNA), "'NA' have been detected from y\n")
		y <- y[!yNA]
		y.id <- y.id[!yNA]
	}
	if(!is.matrix(M)){M <- as.matrix(M); gc()}
	if(nrow(M) != length(M.id))	stop("number of individuals not match between 'M' and 'M.id'.")
	if(ncol(P) != 3)	stop("3 columns ('id', 'sir', 'dam') are required for pedigree.")
	ped <- as.matrix(P)
	ped <- apply(ped, 2, as.character)
	ped.id <- unique(as.character(ped))
	Msub.id <- M.id[!M.id %in% ped.id]
	if(length(Msub.id) == length(M.id))	stop("no shared individuals between 'M.id' and 'P'.")
	if(length(Msub.id)){
		ped <- rbind(ped, cbind(Msub.id, "0", "0"))
		ped.id <- c(Msub.id, ped.id)
	}
	ysub.id <- y.id[!y.id %in% ped.id]
	if(length(ysub.id) == length(y.id))	stop("no shared individuals between 'y.id' and 'P'.")
	if(length(ysub.id)){
		if(verbose)	cat(length(ysub.id), " individuals can not be found in genotype or pedigree\n")
		y.id <- y.id[y.id %in% ped.id]
		y <- y[y.id %in% ped.id]
	}
	if(all(ped.id %in% M.id))	stop("all individuals have been genotyped, no necessaries to fit single-step bayes model.")
	indx <- match(y.id, ytmp.id)

	if(!is.null(X)){
		if(!is.matrix(X))	X <- as.matrix(X)
		if(nrow(X) != length(yNA))	stop("number of individuals not match between 'y' and 'X'.")
		X <- X[indx, , drop=FALSE]
		X_is_num <- apply(X, 2, is.numeric)
		if(!all(X_is_num))	stop("covariates must be a numeric matrix, please use 'model.matrix' to convert.")
		if(!all(apply(X, 2, function(x){unix <- unique(x); if(length(unix) == 1 && unix == 1){FALSE}else{TRUE}})))	stop("please remove intercept from covariates.")
		if(!all(apply(X, 2, function(x){length(unique(x)) > 1})))	stop("covariates should not be a constant.")
	}
	if(!is.null(R)){
		if(!is.matrix(R))	R <- as.matrix(R)
		if(nrow(R) != length(yNA))	stop("number of individuals not match between 'y' and 'R'.")
		R <- R[indx, , drop=FALSE]
		R <- apply(R, 2, as.character)
	}

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

	res <- Bayes(y=y, X=y.M, model=model, Pi=Pi, fold=fold, C=X, R=R, epsl_y_J=y.J, epsl_Gi=Ai.nn, epsl_index=y.Mn.indx, niter=niter, nburn=nburn, windindx=windindx, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, threads=threads, verbose=verbose)
	rm(y.M, y.J, Ai.nn); gc()

	if(length(y.Mn.indx)){
		res$MCMCsamples$g <- as.matrix(c(J, Jn)) %*% res$MCMCsamples$J + rbind(M %*% res$MCMCsamples$alpha, Mn %*% res$MCMCsamples$alpha + res$MCMCsamples$epsilon)
		epsilon <- data.frame(id = Mn.id, epsilon = res$epsilon)
		res$epsilon <- epsilon
	}else{
		warning("all phenotypic individuals have genotype information, thus can't fit imputation errors.")
		res$MCMCsamples$g <- rbind(M %*% res$MCMCsamples$alpha, Mn %*% res$MCMCsamples$alpha)
	}
	res$g <- data.frame(id = c(M.id, Mn.id), gebv = apply(res$MCMCsamples$g, 1, mean))

	e[match(y.id.comb, ytmp.id)] <- res$e
	res$e <- data.frame(id = ytmp.id, e = e)

	if(!is.null(windsize) | !is.null(windnum)){
		WPPA <- res$gwas
		res$gwas <- data.frame(windinfo, WPPA)
	}

	return(res)
}
