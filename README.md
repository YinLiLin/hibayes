# hibayes [![](https://img.shields.io/badge/Issues-%2B-brightgreen.svg)](https://github.com/YinLiLin/hibayes/issues/new) [![](https://img.shields.io/badge/Release-v0.99.1-darkred.svg)](https://github.com/YinLiLin/hibayes)
## Individual and summary level data based BAYES model for Genome-Wide Association and Genomic Prediction

```hibayes``` is an user-friendly tool ([R](https://www.r-project.org) version of [GCTB](http://cnsgenomics.com/software/gctb/#Overview)) to fit BAYES model using individual-level data and summary-level data for both Genomic prediction/selection and Genome-Wide association study, it was desighed to estimate joint effects and genetic parameters for a complex trait, including 1)genetic variance, 2)residual variance, 3)heritability, 4)joint distribution of effect size, 5)phenotype/genetic variance explained (PVE) for single or multiple SNPs, and 6)posterior probability of association of the genomic window (WPPA). The functions are not limited, we will keep on going in enriching ```hibayes``` with more features.

```hibayes``` is developed by [Lilin Yin](https://github.com/YinLiLin) with the support of [Haohao Zhang](https://github.com/hyacz), [Xiaolei Liu](https://github.com/XiaoleiLiuBio), [Jian Zeng](http://researchers.uq.edu.au/researcher/14033) and [Jian Yang](https://researchers.uq.edu.au/researcher/2713). If you have any bug reports or questions, please send an email to Lilin Yin ([ylilin@163.com](mailto:ylilin@163.com)).

## Installation
```hibayes``` can be installed from GitHub as folowing, please ensure ```devtools``` has been installed prior to ```hibayes```.
```r
> devtools::install_github("YinLiLin/hibayes")
```
After installed successfully, type ```library(hibayes)``` to use. The package is on its way to R CRAN, and would be coming soon.

## Usage
### Individual level bayes model
To fit individual level bayes model, the phenotype(n), numeric genotype (n * m, n is the number of individuals, m is the number of SNPs) should be provided. Uses can load the data that coded by other softwares by 'read.table' to fit model. Additionally, we pertinently provide a function ```read_plink``` to load [PLINK binary files](http://zzz.bwh.harvard.edu/plink/binary.shtml) into memory. For example, load the attached tutorial data in ```hibayes```:
```r
> bfile_path = system.file("extdata", "example", package = "hibayes")
> data = read_plink(bfile=bfile_path, mode="A", threads=4)
  ## bfile: the prefix of binary files
  ## model: "A" (addtive) or "D" (dominant)
> pheno = data$pheno
> geno = data$geno
> map = data$map
```

Total 8 bayes models are available currently, including:
 - ***"BayesRR" (ridge regression):*** all SNPs have non-zero effects and share the same variance, equals to GBLUP.
 - ***"BayesA":*** all SNPs have non-zero effects but use different variance which follows an inverse chi-square distribution.
 - ***"BayesLASSO":*** all SNPs have non-zero effects but use different variance which follows an exponential distribution.
 - ***"BayesB":*** only a small part of SNPs (1-pi) have non-zero effects but use different variance which follows an inverse chi-square distribution.
 - ***"BayesBpi":*** the same with "BayesB", but 'pi' is not fixed.
 - ***"BayesC":*** only a small part of SNPs (1-pi) have non-zero effects and share the same variance.
 - ***"BayesCpi":*** the same with "BayesC", but 'pi' is not fixed.
 - ***"BayesR":*** only a small part of SNPs have non-zero effects, but the SNPs are allocated into different groups, each group has the same variance.
 
Type ```?bayes``` to see details of all parameters.

#### (1) Gemonic prediction/selection
```r
> fit <- bayes(y=pheno[, 1], X=geno, pi=0.95, model="BayesB", niter=20000, nburn=10000, outfreq=10, verbose=TRUE)
> SNPeffect <- fit$g
> gebv <- geno %*% SNPeffect    # calculate the estimated genomic breeding value
> pve <- apply(as.matrix(geno), 2, var) * fit$g^2    # the phenotypic variance explained by each SNPs
> nonZeroRate <- fit$nzrate    # the rate of stepping into non-zero effects in MCMC iteration for each SNPs
```
View the results by [CMplot](https://github.com/YinLiLin/R-CMplot) package:
```r
> source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
> CMplot(cbind(map, SNPeffect), type="h", plot.type="m", LOG10=FALSE, ylab="SNP effect")
```
<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/Figure/3.jpg">
<img src="Figure/3.jpg" height="385px" width="900px">
</a>
</p>

```r
> CMplot(cbind(map, pve), type="p", plot.type="m", LOG10=FALSE, ylab="Phenotypic variance explained")
```
<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/Figure/3.jpg">
<img src="Figure/3.jpg" height="385px" width="900px">
</a>
</p>

#### (2) Gemone-Wide association study
```r
> fit <- bayes(y=pheno[, 1], X=geno, map=map, windsize=1e6, model="BayesR", niter=20000, nburn=10000, outfreq=10)
```


## Not done yet
