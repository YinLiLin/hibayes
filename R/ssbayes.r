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
#' \item "BayesB": only a small proportion of SNPs (1-pi) have non-zero effects, and take different variance which follows an inverse chi-square distribution. 
#' \item "BayesBpi": the same with "BayesB", but 'pi' is not fixed. 
#' \item "BayesC": only a small proportion of SNPs (1-pi) have non-zero effects, and share the same variance. 
#' \item "BayesCpi": the same with "BayesC", but 'pi' is not fixed. 
#' \item "BayesL": BayesLASSO, all SNPs have non-zero effects, and take different variance which follows an exponential distribution.
#' \item "BSLMM": all SNPs have non-zero effects, and take the same variance, but a small proportion of SNPs have additional shared variance. 
#' \item "BayesR": only a small proportion of SNPs have non-zero effects, and the SNPs are allocated into different groups, each group has the same variance. 
#' }
#' @param map (optional, only for GWAS) the map information of genotype, at least 3 columns are: SNPs, chromosome, physical position. 
#' @param pi vector, the proportion of zero effect and non-zero effect SNPs, the first value must be the proportion of non-effect markers.
#' @param fold proportion of variance explained for groups of SNPs, the default is c(0, 0.0001, 0.001, 0.01).
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
#' @param threads number of threads used for OpenMP.
#' @param verbose whether to print the iteration information.
#'
#' @return
#' the function returns a list containing
#' \describe{
#' \item{$mu}{the regression intercept}
#' \item{$pi}{estimated proportion of zero effect and non-zero effect SNPs}
#' \item{$beta}{estimated coefficients for all covariates}
#' \item{$r}{estimated environmental random effects}
#' \item{$vr}{estimated variance for all environmental random effect}
#' \item{$vg}{estimated genetic variance}
#' \item{$ve}{estimated residual variance}
#' \item{$alpha}{estimated effect size of all markers}
#' \item{$modfreq}{the frequency for markers to be included in the model during MCMC iteration, also known as posterior inclusive probability (PIP)}
#' \item{$g}{data.frame, the first column is the list of individual id, the second column is the genomic estimated breeding value for all individuals, including genotyped and non-genotyped.}
#' \item{$gwas}{WPPA is defined to be the window posterior probability of association, it is the ratio of the number of iterations that \eqn{Pw} (the proportion of the total genetic variance explained by the window \eqn{w}) > 1% divided by the total number of MCMC iterations, WGVE is the explained genetic variance for each window}
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
#' # For GS/GP
#' fit = ssbayes(y=pheno[, 2], y.id=pheno.id, M=geno, M.id=geno.id, P=ped, 
#' 				 model="BayesR", niter=200, nburn=100, outfreq=10)
#' # For GWAS
#' # fit = ssbayes(y=pheno[, 2], y.id=pheno.id, M=geno, M.id=geno.id, P=ped, 
#' # 			  map=map, windsize=1e6, model="BayesCpi")
#' 
#' # Add fixed effects, covariates, and random effect
#' X <- model.matrix.lm(~as.numeric(scale)+as.factor(sex), data=pheno, na.action = "na.pass")
#' X <- X[, -1] #remove the intercept
#' # fit = ssbayes(..., X=X, R=pheno[,c("group")], ...)
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
    model = c("BayesCpi", "BayesA", "BayesL", "BSLMM", "BayesR", "BayesB", "BayesC", "BayesBpi", "BayesRR"),
	map = NULL,
    pi = NULL,
    fold = NULL,
    niter = 20000,
    nburn = 14000,
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
	threads = 4,
    verbose = TRUE
){
	if(!is.null(windsize)){
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
		map <- matrix(as.numeric(map), nrow(map))
		chr <- chr[order(map[,1])]
		if(max(map[,2]) < windsize)	stop("Maximum of physical position is smaller than wind size.")
		windindx <- cutwind(map[,1], map[,2], windsize)
		windrange <- do.call(rbind, tapply(map[, 2], windindx, range))
		windsnpN <- tapply(map[, 2], windindx, length)
		windchr <- unique(chr)[match(tapply(map[, 1], windindx, unique), unique(sort(map[,1])))]
		windinfo <- data.frame(paste0("wind", 1:max(windindx)), windchr, windsnpN, windrange)
		colnames(windinfo) <- c("Wind", "Chr", "N", "Start", "End")
	}else{
		windindx <- NULL
	}
	set.seed(seed)
	model <- match.arg(model)
	if(is.null(pi)){
		if(model == "BayesR"){
			pi <- c(0.95, 0.02, 0.02, 0.01)
			if(is.null(fold))	fold <- c(0, 0.0001, 0.001, 0.01)
		}else{
			pi <- c(0.95, 0.05)
		}
	}
	y.id <- as.character(y.id)
	M.id <- as.character(M.id)
	if(!is.numeric(y)){
		y <- as.matrix(y)[,1,drop=TRUE]
		if(is.character(y))	stop("y is not a vector of digital values.")
	}
	if(length(y) != length(y.id))	stop("number of individuals not match between 'y' and 'y.id'.")
	yNA <- is.na(y)
	ytmp.id <- y.id
	if(sum(yNA) != 0){
		if(verbose)	cat(sum(yNA), "'NA' have been removed from y\n")
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
		if(verbose)	cat(length(ysub.id), " individuals that can't be found in genotype or pedigree have been removed\n")
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
	Mn <- geno_impute(A.ng, M, threads)
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

	res = Bayes(y=y, X=y.M, model=model, pi=pi, fold=fold, C=X, R=R, epsl_y_J=y.J, epsl_Gi=Ai.nn, epsl_index=y.Mn.indx, niter=niter, nburn=nburn, windindx=windindx, wppa=wppa, vg=vg, dfvg=dfvg, s2vg=s2vg, ve=ve, dfve=dfve, s2ve=s2ve, outfreq=outfreq, threads=threads, verbose=verbose)
	rm(y.M, y.J, Ai.nn); gc()

	if(length(y.Mn.indx)){
		g <- c(J, Jn) * res$eps_J + c(as.vector(M %*% res$alpha), as.vector(Mn %*% res$alpha)) + c(rep(0, nrow(M)), res$eps_R)
	}else{
		g <- c(as.vector(M %*% res$alpha), as.vector(Mn %*% res$alpha))
	}
	g <- data.frame(id = c(M.id, Mn.id), g = g)
	
	res = tail(res, length(res) - 2)
	if(!is.null(windsize)){
		WPPA <- res$wppa
		WGVE <- res$wgve
		res <- head(res, length(res) - 2)
		res$g <- g
		res$gwas <- data.frame(windinfo, WPPA, WGVE)
	}else{
		res$g <- g	
	}
	return(res)
}
