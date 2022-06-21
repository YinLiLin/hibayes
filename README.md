# hibayes 
[![GitHub issues](https://img.shields.io/github/issues/YinLiLin/hibayes?color=green)](https://github.com/YinLiLin/hibayes/issues/new)  [![CRAN Version](https://www.r-pkg.org/badges/version/hibayes?color=yellow)](https://CRAN.R-project.org/package=hibayes) [![](https://img.shields.io/badge/GitHub-1.1.0-blueviolet.svg)]() ![](http://cranlogs.r-pkg.org/badges/grand-total/hibayes?color=red) [![](https://cranlogs.r-pkg.org/badges/last-month/hibayes)](https://CRAN.R-project.org/package=hibayes) <a href="https://hits.seeyoufarm.com"/><img src="https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FYinLiLin%2Fhibayes"/></a>
## Individual-Level, Summary-Level and Single-Step Bayesian Regression Models for Genomic Prediction and Genome-Wide Association Studies

**```hibayes```** (say 'Hi' to Bayes) is an user-friendly [R](https://www.r-project.org) package to fit 3 types of Bayesian models using **[individual-level](#1-individual-level-bayesian-model)**, **[summary-level](#2-summary-level-bayesian-model)**, and **[individual plus pedigree-level](#3-single-step-bayesian-model)** (single-step) data for both Genomic prediction/selection (GS) and Genome-Wide Association Study (GWAS), it was designed to estimate joint effects and genetic parameters for a complex trait, including:  
**(1)** fixed effects and coefficients of covariates  
**(2)** environmental random effects, and its corresponding variance  
**(3)** genetic variance  
**(4)** residual variance  
**(5)** heritability  
**(6)** genomic estimated breeding values (*GEBV*) for both genotyped and non-genotyped individuals  
**(7)** SNP effect size  
**(8)** phenotype/genetic variance explained (*PVE*) for single or multiple SNPs  
**(9)** posterior probability of association of the genomic window (*WPPA*)  
**(10)** posterior inclusive probability (*PIP*)  
The functions are not limited, we will keep on going in enriching **```hibayes```** with more features.

**```hibayes```** is written in C++ by aid of Rcpp and RcppArmadillo, some time-consuming functions are enhanced with [LAPACK](http://www.netlib.org/lapack/) package, it is recommended to run **```hibayes```** in [**MRO**](https://mran.microsoft.com) instead of **R**, as the BLAS/LAPACK library can be accelerated automatically in multi-threads by MKL library, which would significantly reduce computation time. 

***If you have any bug reports or questions, please feed back :point_right:[here](https://github.com/YinLiLin/hibayes/issues/new):point_left:.***

## :toolbox: Software tools for genetic analyses and genomic breeding

<table>
    <tr>
	<td><g-emoji class="g-emoji" alias="postbox" fallback-src="https://github.githubassets.com/images/icons/emoji/unicode/1f4ee.png">📮</g-emoji> <strong><a href="https://github.com/xiaolei-lab/rMVP">rMVP</a></strong>: Efficient and easy-to-use GWAS tool.</td>
        <td><g-emoji class="g-emoji" alias="four_leaf_clover" fallback-src="https://github.githubassets.com/images/icons/emoji/unicode/1f340.png">🍀</g-emoji> <strong><a href="https://github.com/xiaolei-lab/SIMER">SIMER</a></strong>: data simulation for life science and breeding.</td>
    </tr>
    <tr>
        <td><g-emoji class="g-emoji" alias="mailbox" fallback-src="https://github.githubassets.com/images/icons/emoji/unicode/1f4eb.png">📫</g-emoji> <strong><a href="https://www.hiblup.com/" rel="nofollow">HIBLUP</a></strong>: Versatile and easy-to-use GS toolbox.</td>
        <td><g-emoji class="g-emoji" alias="mountain_snow" fallback-src="https://github.githubassets.com/images/icons/emoji/unicode/1f3d4.png">🏔️</g-emoji> <strong><a href="http://iswine.iomics.pro/pig-iqgs/iqgs/index" rel="nofollow">ISWINE</a></strong>: an omics knowledgebase for swine.</td>
    </tr>
    <tr>
        <td><g-emoji class="g-emoji" alias="biking_man" fallback-src="https://github.githubassets.com/images/icons/emoji/unicode/1f6b4-2642.png">🚴&zwj;♂️</g-emoji> <strong><a href="https://github.com/YinLiLin/KAML">KAML</a></strong>: Advanced GS method for complex traits.</td>
        <td><g-emoji class="g-emoji" alias="bar_chart" fallback-src="https://github.githubassets.com/images/icons/emoji/unicode/1f4ca.png">📊</g-emoji> <strong><a href="https://github.com/YinLiLin/CMplot">CMplot</a></strong>: A drawing tool for genetic analyses.</td>
    </tr>
    <tr>
        <td><g-emoji class="g-emoji" alias="swimmer" fallback-src="https://github.githubassets.com/images/icons/emoji/unicode/1f3ca.png">🏊</g-emoji> <strong><a href="https://github.com/YinLiLin/hibayes">hibayes</a></strong>: A Bayesian-based GWAS and GS tool.</td>
        <td></td>
    </tr>
</table>

## Installation
The stable version of **```hibayes```** can be accessed from CRAN, type the following script to install:
```r
> install.packages("hibayes")
```
After installed successfully, type *```library(hibayes)```* to use.  
The latest version of **```hibayes```**  in development can be installed from GitHub as following, please ensure **```devtools```** has been installed prior to installing **```hibayes```**.
```r
> devtools::install_github("YinLiLin/hibayes")
```
## Citing the package
Yin LL, Zhang HH, Li XY, Zhao SH, Liu XL. [hibayes: An R Package to Fit Individual-Level, Summary-Level and Single-Step Bayesian Regression Models for Genomic Prediction and Genome-Wide Association Studies](https://www.biorxiv.org/content/10.1101/2022.02.12.480230v2), ***bioRxiv*** (2022), doi: 10.1101/2022.02.12.480230.
## Usage
### 1. Individual level Bayesian model
To fit individual level Bayesian model (*```bayes```*), at least the phenotype(***n***), numeric genotype (***n*** * ***m***, ***n*** is the number of individuals, ***m*** is the number of SNPs) should be provided. Users can load the phenotype and genotype data that coded by other softwares by *```read.table```* to fit model, note that 'NA' is not allowed in genotype data:
```r
> pheno = read.table("your_pheno.txt")
> geno = read.table("your_geno.txt") # genotype should be coded in digits (either 0, 1, 2 or -1, 0, 1 is acceptable)
> geno.id = read.table("your_genoid.txt")
> # the order of individuals should be exactly the same between phenotype and genotype
> pheno = pheno[match(geno.id[, 1], pheno[, 1]), ]   # supposing the first column is the individual id
```
Additionally, we pertinently provide a function *```read_plink```* to load [PLINK binary files](http://zzz.bwh.harvard.edu/plink/binary.shtml) into memory. For example, load the attached tutorial data in **```hibayes```**:
```r
> bfile_path = system.file("extdata", "example", package = "hibayes")
> data = read_plink(bfile=bfile_path, mode="A", threads=4)
> # bfile: the prefix of binary files
> # mode: "A" (additive) or "D" (dominant)
> pheno = data$fam
> nrow(pheno) # number of individuals
[1] 4798
> head(pheno)
1 F1 Ind1 0 0 0  -2.0464930
2 F2 Ind2 0 0 0  -7.7857772
3 F3 Ind3 0 0 0   3.4085493
4 F4 Ind4 0 0 0 -10.3651915
5 F5 Ind5 0 0 0   0.4604928
6 F6 Ind6 0 0 0  -9.7750056
> geno = data$geno
> dim(geno) # number of individuals and markers
[1] 4798 7385
> geno[1:5,1:5]
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    0    0    0    1
[2,]    0    0    0    0    0
[3,]    0    0    0    0    1
[4,]    1    1    0    1    0
[5,]    1    0    1    1    0
> map = data$map
> head(map)
   SNP Chr     Pos A1 A2
1 snp1   1 1198554  T  C
2 snp2   1 1720354  A  G
3 snp3   1 1825948  A  C
4 snp4   1 3428453  G  A
5 snp5   1 4195032  T  C
6 snp6   1 4412357  T  C
```
In this function, missing genotype will be replaced by the major genotype of each allele. **```hibayes```** will code the genotype ***A1A1*** as 2, ***A1A2*** as 1, and ***A2A2*** as 0, where ***A1*** is the first allele of each marker in *\*.bim* file, therefore the estimated effect size is on ***A1*** allele, users should pay attention to it when a process involves marker effect. 

By default, the memory-mapped files are directed into work directory, users could redirect to new path as following:
```r
> data <- read_plink(bfile=bfile_path, out="./test")
> # directly use the genotype for the next time, no need to use 'read_plink' again:
> geno <- attach.big.matrix("./test.desc")
> map <- read.table("./test.map", header=TRUE)
```
For **fixed effects** and **covariates**, please use *```model.matrix.lm()```* to make the model matrix prior to fitting models:
```r
> # For fixed effects, use 'as.factor', eg. 'sex'. 
> # For covariates, use 'as.numeric', eg. 'weight'.
> X <- model.matrix.lm(~as.factor(sex)+as.numeric(weight), data=pheno, na.action = "na.pass")
> X <- X[, -1] #remove the intercept
```
We can also fit the interactions between environmental effects, for example:
```r
> X <- model.matrix.lm(~as.factor(sex)+as.numeric(weight)+
   as.factor(group:location), data=pheno, na.action = "na.pass")
> X <- X[, -1] #remove the intercept
```
For **random effects**, no needs to convert, just pick them out from the phenotype data, eg. 'group', 'location':
```r
> R <- pheno[, c("group", "location")]
```
Then add it into different models like:
```r
fit <- bayes(..., X=X, R=R, ...)    # bayes model
fit <- ssbayes(..., X=X, R=R, ...)  # single-step bayes model
```
Following methods are available currently, including:
 - ***"BayesRR":*** Bayesian Ridge Regression, all SNPs have non-zero effects and share the same variance, equals to RRBLUP or GBLUP. 
 - ***"BayesA":*** all SNPs have non-zero effects, and take different variance which follows an inverse chi-square distribution. 
 - ***"BayesB":*** only a small proportion of SNPs (1-Pi) have non-zero effects, and take different variance which follows an inverse chi-square distribution. 
 - ***"BayesBpi":*** the same with "BayesB", but 'Pi' is not fixed. 
 - ***"BayesC":*** only a small proportion of SNPs (1-Pi) have non-zero effects, and share the same variance. 
 - ***"BayesCpi":*** the same with "BayesC", but 'Pi' is not fixed. 
 - ***"BayesL":*** BayesLASSO, all SNPs have non-zero effects, and take different variance which follows an exponential distribution.
 - ***"BSLMM":*** all SNPs have non-zero effects, and take the same variance, but a small proportion of SNPs have additional shared variance. 
 - ***"BayesR":*** only a small proportion of SNPs have non-zero effects, and the SNPs are allocated into different groups, each group has the same variance. 

Type *```?bayes```* to see details of all parameters.

#### (a) Gemonic prediction/selection
```r
> fit <- bayes(y = pheno[, 6], M = geno, model = "BayesCpi", niter = 20000, 
		Pi = c(0.95, 0.05), nburn = 12000)
> str(fit)	# overview of the returns
List of 10
 $ Vg         : num 3.91
 $ Ve         : num 24.9
 $ h2         : num 0.136
 $ mu         : num 0.405
 $ alpha      : num [1:7385, 1] 0.00 -3.17e-04 2.56e-04 -2.36e-05 0.00 ...
 $ pi         : num [1:2, 1] 0.99762 0.00238
 $ g          : num [1:4798] -2.98 -5.45 1.52 -2.34 -3.26 ...
 $ e          : num [1:4798] 0.527 -2.744 1.479 -8.43 3.317 ...
 $ pip        : num [1:7385, 1] 0.000125 0.0005 0.000625 0.0005 0.000375 ...
 $ MCMCsamples:List of 7
  ..$ Vg   : num [1, 1:400] 3.68 3.83 3.88 3.9 3.93 ...
  ..$ Ve   : num [1, 1:400] 25.3 24.6 24.5 24.3 25.4 ...
  ..$ h2   : num [1, 1:400] 0.127 0.135 0.137 0.138 0.134 ...
  ..$ mu   : num [1, 1:400] 0.4576 0.4069 0.4688 -0.0538 -0.0587 ...
  ..$ alpha: num [1:7385, 1:400] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ pi   : num [1:2, 1:400] 0.99648 0.00352 0.99881 0.00119 0.99812 ...
  ..$ g    : num [1:4798, 1:400] -2.02 -6.15 1.51 -2.42 -2.28 ...
> SNPeffect <- fit$alpha	# get the estimated SNP effects for markers
> gebv <- fit$g		# get the genomic estimated breeding values (GEBV) for all individuals
> pve <- apply(as.matrix(geno), 2, var) * (fit$alpha^2) / var(pheno[, 6])    # the phenotypic variance explained (pve) for each SNPs
```
Note that the standard deviation of all estimated unknow parameters ('p') can be obtained by the function *```apply(fit$MCMCsamples$p, 1, sd)```*, for example:
```r
> SNPeffect_SD <- apply(fit$MCMCsamples$alpha, 1, sd)	# get the SD of estimated SNP effects for markers
> gebv_pev <- apply(fit$MCMCsamples$g, 1, var)	# get the prediction error variance (PEV) of estimated breeding values
```
View the results by [CMplot](https://github.com/YinLiLin/R-CMplot) package:
```r
> source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
> CMplot(cbind(map[,1:3], SNPeffect), type="h", plot.type="m", LOG10=FALSE, ylab="SNP effect")
```
<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/hibayes/master/figure/SNPeff.jpg">
<img src="figure/SNPeff.jpg" height="350px" width="900px">
</a>
</p>

```r
> highlight <- map[pve>0.001,1]
> CMplot(cbind(map[, 1:3], 100 * pve), type = "h", plot.type = "m",
	LOG10 = FALSE, ylab = "Phenotypic variance explained (%)",
	highlight = highlight, highlight.text = highlight)
```
<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/hibayes/master/figure/PVE.jpg">
<img src="figure/PVE.jpg" height="340px" width="900px">
</a>
</p>

#### (b) Gemone-Wide association study
**WPPA** is defined to be the window posterior probability of association, it is estimated by counting the number of MCMC samples in which the effect size is nonzero for at least one SNP in the window. To run GWAS, *```map```* should be provided, every marker should have clear physical position for the downstream genome cutting, and also the argument *```windsize```* or *```windnum```* should be specified, the argument *```windsize```* is used to control the size of the windows, the number  of markers in a window is not fixed. Contrarily, the argument *```windnum```*, e.g. windnum = 10, can be used to control the fixed number of markers in a window, the size for the window is not fixed for this case.
```r
> fit <- bayes(y = pheno[, 6], M = geno, model = "BayesCpi", niter = 20000,
		Pi = c(0.95, 0.05), nburn = 12000, seed = 666666,
		map = map, windsize = 1e6)
> gwas <- fit$gwas
> head(gwas)
   Wind Chr N   Start     End     WPPA
1 wind1   1 3 1198554 1825948 0.001250
2 wind2   1 1 3428453 3428453 0.000500
3 wind3   1 8 4195032 4916148 0.013250
4 wind4   1 7 5109162 5881216 0.007500
5 wind5   1 3 6705835 6952985 0.003125
6 wind6   1 7 7075618 7863025 0.004375
```
View the results by [CMplot](https://github.com/YinLiLin/R-CMplot) package:
```r
> highlight <- gwas[(1 - gwas[, "WPPA"]) < 0.01, 1]
> CMplot(cbind(gwas[, c(1, 2, 4)], 1 - gwas[, "WPPA"]), type = "h",
	plot.type = "m", LOG10 = TRUE, threshold = 0.01, ylim = c(0, 5),
	ylab = expression(-log[10](1 - italic(WPPA))), highlight = highlight,
	highlight.col = NULL, highlight.text = highlight)
```
<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/hibayes/master/figure/WPPA.jpg">
<img src="figure/WPPA.jpg" height="350px" width="900px">
</a>
</p>

One can also derive the association significance from the posterior inclusive probability (**PIP**) of each SNP for certain genome region in whole MCMC procedure.
```r
> data <- cbind(map[, 1:3], (1 - fit[["pip"]]))
> chr5 <- data[data[, 2] == 5, ]
> # visualize the results
> CMplot(chr5, plot.type = "m", width = 9, height = 5, threshold = 0.01,
	ylab = expression(-log[10](1 - italic(PIP))), LOG10 = TRUE,
	amplify = FALSE)
```
<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/hibayes/master/figure/PIP.jpg">
<img src="figure/PIP.jpg" height="380px" width="750px">
</a>
</p>

-----

### 2. Summary level Bayesian model
To fit summary level data based Bayesian model (*```sbayes```*), the variance-covariance matrix calculated from the reference panel (can be done by **```hibayes```**), and summary data in [COJO](https://cnsgenomics.com/software/gcta/#COJO) file format should be provided. Specially, if the summary data is derived from reference panel, means that all data come from the same population, then summary data level based Bayesian model equals to the individual level Bayesian model. 

The available methods for *```sbayes```* include ***"BayesRR", "BayesA", "BayesLASSO", "BayesB", "BayesBpi", "BayesC", "BayesCpi", "BayesR", "CG" (conjuction gradient)***. For 'CG' method, parameter *```lambda```* should be assigned with *```m * (1 / h2 - 1)```*, where ***m*** is the total number of SNPs and ***h2*** is the heritability that can be estimated from LD score regression analysis using the summary data.

#### Step1: construct full/sparse LD variance-covariance matrix
Sparse matrix could significantly reduce the memory cost by setting some of elements of full matrix to zero, on condition that *```n*r^2 < chisq```*, where ***n*** is the number of individuals, ***r*** is the LD correlation of pairs of SNPs, some low LD values would be replaced by 0.
```r
> # load reference panel
> bfile_path = system.file("extdata", "geno", package = "hibayes")
> data = read_plink(bfile_path)
> geno = data$geno
> map = data$map
> # construct LD variance-covariance matrix
> ldm1 = ldmat(geno, threads=4)   #chromosome wide full ld matrix
> ldm2 = ldmat(geno, chisq=5, threads=4)   #chromosome wide sparse ld matrix
> ldm3 = ldmat(geno, map, ldchr=FALSE, threads=4)   #chromosome block ld matrix
> ldm4 = ldmat(geno, map, ldchr=FALSE, chisq=5, threads=4)   #chromosome block + sparse ld matrix
```
From ```ldm1``` to ```ldm4```, the memory cost less, but the model stability of *```sbayes```* would be worse.

#### Step2: fit SBayes model
If the order of SNPs in variance-covariance matrix is not consistent with the order in summary data file, prior adjusting is necessary:
```r
> sumstat_path = system.file("extdata", "geno.ma", package = "hibayes")
> sumstat = read.table(sumstat_path, header=TRUE)
> head(sumstat)
   SNP A1 A2    MAF    BETA     SE      P NMISS
1 snp1  G  A 0.3000  0.1783 0.3215 0.5813    60
2 snp2  T  G 0.3667  0.1451 0.2735 0.5978    60
3 snp3  A  G 0.3167  0.3815 0.3363 0.2613    60
4 snp4  C  A 0.3417  0.3699 0.3286 0.2649    60
5 snp5  T  G 0.3250  0.5380 0.3522 0.1321    60
6 snp6  T  G 0.3000 -0.2677 0.3346 0.4270    60
> sumstat = sumstat[match(map[,1], sumstat[,1]), ]  # match the order of SNPs
```
Note that **```hibayes```** only use the 'BETA', 'SE' and 'NMISS' columns.  
Type *```?sbayes```* to see details of all parameters.
#### (a) Gemonic prediction/selection
```r
> fit = sbayes(sumstat=sumstat, ldm=ldm1, model="BayesCpi", niter=20000, nburn=12000)
```
#### (b) Gemone-Wide association study
```r
> fit = sbayes(sumstat=sumstat, ldm=ldm1, map=map, model="BayesCpi", windsize=1e6, niter=20000, nburn=12000)
```
Note that the standard deviation of all estimated unknow parameters ('p') can be obtained by the function *```apply(fit$MCMCsamples$p, 1, sd)```*, for example:
```r
> SNPeffect <- fit$alpha	# get the estimated SNP effects for markers
> SNPeffect_SD <- apply(fit$MCMCsamples$alpha, 1, sd)	# get the SD of estimated SNP effects for markers
```
-----

### 3. Single-step Bayesian model
To fit single-step Bayesian model (*```ssbayes```*), at least the phenotype(***n1***, the number of phenotypic individuals), numeric genotype (***n2*** * ***m***, ***n2*** is the number of genotyped individuals, ***m*** is the number of SNPs), and pedigree information (***n3*** * ***3***, the three columns are "id" "sir" "dam" orderly) should be provided, ***n1***, ***n2***, ***n3*** can be different, all the individuals in pedigree will be predicted, including genotyped and non-genotyped, therefore the total number of predicted individuals depends on the number of unique individuals in pedigree.  
For example, load the attached tutorial data in **```hibayes```**:
```r
> # load phenotype file
> pheno_file_path = system.file("extdata", "pheno.txt", package = "hibayes")
> pheno = read.table(pheno_file_path, header=TRUE)
> nrow(pheno) # number of individuals
[1] 100
> head(pheno)
    id          y scale group sex
1 ind1 -0.5796816  0.77    g1   m
2 ind2 -2.0224628 -1.02    g2   m
3 ind3 -1.4807132  0.52    g1   f
4 ind4 -3.0303065 -1.05    g4   m
5 ind5  2.1881874  2.06    g3   m
6 ind6 -3.2110719 -1.94    g4   m
> # load pedigree file
> pedigree_file_path = system.file("extdata", "ped.txt", package = "hibayes")
> ped = read.table(pedigree_file_path, header=TRUE)
> head(ped)
     id  sire   dam
1 ind20  <NA>  <NA>
2 ind21 ind17 ind12
3 ind22  ind3 ind20
4 ind23  ind4 ind16
5 ind24  ind1 ind14
6 ind25  ind5 ind13
> # load genotype file
> bfile_path = system.file("extdata", "geno", package = "hibayes")
> data = read_plink(bfile=bfile_path, mode="A", threads=4)
> # bfile: the prefix of binary files
> # mode: "A" (addtive) or "D" (dominant)
> fam = data$fam
> geno = data$geno
> map = data$map
> dim(geno) # number of genotyped individuals and markers
[1]   60 1000
> # get the of phenotype and genotype id
> geno.id = fam[, 2]
> pheno.id = pheno[, 1]
```
For fixed effects, covariates, and environmental random effects, please refer to the chapter of [*bayes*](#1-individual-level-bayesian-model) model.
```r
> X <- model.matrix.lm(~as.factor(sex)+as.numeric(scale), data=pheno, na.action = "na.pass")
> X <- X[, -1]	# fixed effects and covariates
> R <- pheno[, c("group")]	# environmental random effects
```
***NOTE:*** for *```ssbayes```* model, there is no NEED to adjust the order of id in different files.

The available methods for *```ssbayes```* model are consistent with *```bayes```* model, except for "BSLMM". Type *```?ssbayes```* to see details of all parameters.

#### (a) Gemonic prediction/selection
```r
> fit = ssbayes(y=pheno[, 2], y.id=pheno.id, M=geno, M.id=geno.id, P=ped, 
				X=X, R=R, model="BayesR", niter=20000, nburn=12000)
> str(fit)	# overview of the returns
List of 16
 $ Vr         : num [1, 1] 0.464
 $ Vg         : num 0.568
 $ Ve         : num 0.283
 $ h2         : num 0.433
 $ mu         : num -0.947
 $ beta       : num [1:2, 1] 0.0798 1.0948
 $ alpha      : num [1:1000, 1] 0.00379 -0.00721 0.00186 -0.00282 0.00267 ...
 $ pi         : num [1:2, 1] 0.893 0.107
 $ Veps       : num 0.146
 $ J          : num 0.208
 $ epsilon    :'data.frame':	40 obs. of  2 variables:
  ..$ id     : chr [1:40] "ind20" "ind17" "ind12" "ind3" ...
  ..$ epsilon: num [1:40] 0.00261 0.34582 -0.23766 -0.1511 -0.32559 ...
 $ r          :'data.frame':	5 obs. of  2 variables:
  ..$ Levels    : chr [1:5] "g1" "g2" "g3" "g4" ...
  ..$ Estimation: num [1:5] -0.22552 0.00772 -0.27953 -0.4405 1.19692
 $ g          :'data.frame':	100 obs. of  2 variables:
  ..$ id  : chr [1:100] "ind41" "ind42" "ind43" "ind44" ...
  ..$ gebv: num [1:100] -0.705 -0.521 -0.23 -0.23 0.417 ...
 $ e          :'data.frame':	100 obs. of  2 variables:
  ..$ id: chr [1:100] "ind1" "ind2" "ind3" "ind4" ...
  ..$ e : num [1:100] -0.145 -0.0251 0.0337 -0.2476 0.4869 ...
 $ pip        : num [1:1000, 1] 0.092 0.1008 0.0925 0.0959 0.0877 ...
 $ MCMCsamples:List of 13
  ..$ Vr     : num [1, 1:400] 0.407 0.472 0.491 0.537 0.454 ...
  ..$ Vg     : num [1, 1:400] 0.505 0.691 0.713 0.484 0.51 ...
  ..$ Ve     : num [1, 1:400] 0.372 0.285 0.205 0.296 0.194 ...
  ..$ h2     : num [1, 1:400] 0.393 0.477 0.506 0.367 0.441 ...
  ..$ mu     : num [1, 1:400] -1.252 -1.515 -0.661 -1.256 -0.839 ...
  ..$ beta   : num [1:2, 1:400] -0.1551 1.0734 0.0364 1.1592 0.1838 ...
  ..$ alpha  : num [1:1000, 1:400] 0 0.0513 0 0 0 ...
  ..$ pi     : num [1:2, 1:400] 0.822 0.178 0.829 0.171 0.799 ...
  ..$ Veps   : num [1, 1:400] 0.0936 0.1755 0.153 0.118 0.038 ...
  ..$ J      : num [1, 1:400] -0.00933 0.35316 0.62347 0.66544 0.85299 ...
  ..$ epsilon: num [1:40, 1:400] 0.2484 0.2465 -0.0985 -0.3612 0.2507 ...
  ..$ r      : num [1:5, 1:400] 0.14696 0.62815 0.31242 -0.00248 1.5989 ...
  ..$ g      : num [1:100, 1:400] -0.3187 -0.7824 0.3202 0.0458 -0.0824 ...
```
#### (b) Gemone-Wide association study
```r
> fit = ssbayes(y=pheno[, 2], y.id=pheno.id, M=geno, M.id=geno.id, P=ped, 
				X=X, R=R, map=map, windsize=1e6, model="BayesCpi")
```
Note that the standard deviation of all estimated unknow parameters ('p') can be obtained by the function *```apply(fit$MCMCsamples$p, 1, sd)```*, for example:
```r
> SNPeffect <- fit$alpha	# get the estimated SNP effects for markers
> SNPeffect_SD <- apply(fit$MCMCsamples$alpha, 1, sd)	# get the SD of estimated SNP effects for markers
> gebv <- fit$g		# get the genomic estimated breeding values (GEBV) for all individuals
> gebv_pev <- apply(fit$MCMCsamples$g, 1, var)	# get the prediction error variance (PEV) of GEBV
```
## Citing the methods in package
For *```bayes```* model, please cite following papers:
```
1. Meuwissen, Theo HE, Ben J. Hayes, and Michael E. Goddard. "Prediction of total genetic value using genome-wide dense marker maps." Genetics 157.4 (2001): 1819-1829.
2. de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D., and Calus, M. P. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. Genetics, 193(2), 327-345.
3. Habier, David, et al. "Extension of the Bayesian alphabet for genomic selection." BMC bioinformatics 12.1 (2011): 1-12.
4. Yi, Nengjun, and Shizhong Xu. "Bayesian LASSO for quantitative trait loci mapping." Genetics 179.2 (2008): 1045-1055.
5. Zhou, Xiang, Peter Carbonetto, and Matthew Stephens. "Polygenic modeling with Bayesian sparse linear mixed models." PLoS genetics 9.2 (2013): e1003264.
6. Moser, Gerhard, et al. "Simultaneous discovery, estimation and prediction analysis of complex traits using a Bayesian mixture model." PLoS genetics 11.4 (2015): e1004969.
```
For *```sbayes```* model, please cite following papers:
```
Lloyd-Jones, Luke R., et al. "Improved polygenic prediction by Bayesian multiple regression on summary statistics." Nature communications 10.1 (2019): 1-11.
```
For *```ssbayes```* model, please cite following papers:
```
1. Fernando, Rohan L., Jack CM Dekkers, and Dorian J. Garrick. "A class of Bayesian methods to combine large numbers of genotyped and non-genotyped animals for whole-genome analyses." Genetics Selection Evolution 46.1 (2014): 1-13.
2. Henderson, C.R.: A simple method for computing the inverse of a numerator relationship matrix used in prediction of breeding values. Biometrics 32(1), 69-83 (1976).
```
