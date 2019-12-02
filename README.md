# hibayes [![](https://img.shields.io/badge/Issues-%2B-brightgreen.svg)](https://github.com/YinLiLin/hibayes/issues/new) [![](https://img.shields.io/badge/Release-v0.99.1-darkred.svg)](https://github.com/YinLiLin/hibayes)
## Individual and summary level data based BAYES model for Genome-Wide Association and Genomic Prediction

```hibayes``` is an user-friendly tool ([R](https://www.r-project.org) version of [GCTB](http://cnsgenomics.com/software/gctb/#Overview)) to fit BAYES model using individual-level data and summary-level data for both Genomic prediction/selection and Genome-Wide association study, it was desighed to estimate joint effects and genetic parameters for a complex trait, including 1)genetic variance, 2)residual variance, 3)heritability, 4)joint distribution of effect size, 5)phenotype/genetic variance explained (PVE) for single or multiple SNPs, and 6)posterior probability of association of the genomic window (WPPA). The functions are not limited, we will keep on going in enriching ```hibayes``` with more features.

```hibayes``` is developed by [Lilin Yin](https://github.com/YinLiLin) with the support of [Haohao Zhang](https://github.com/hyacz), [Xiaolei Liu](https://github.com/XiaoleiLiuBio), [Jian Zeng](http://researchers.uq.edu.au/researcher/14033) and [Jian Yang](https://researchers.uq.edu.au/researcher/2713). If you have any bug reports or questions, please feed back :point_right:[here](https://github.com/YinLiLin/hibayes/issues/new):point_left:.

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
> geno = as.matrix(data$geno)
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
> pve <- apply(geno,2,var) * (fit$g^2) / var(pheno[,1])    # the phenotypic variance explained for each SNPs
> nonZeroRate <- fit$nzrate    # the rate of stepping into non-zero effects in MCMC iteration for each SNPs
```
View the results by [CMplot](https://github.com/YinLiLin/R-CMplot) package:
```r
> source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
> CMplot(cbind(map, SNPeffect), type="h", plot.type="m", LOG10=FALSE, ylab="SNP effect")
```
<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/hibayes/master/figure/1.jpg">
<img src="figure/1.jpg" height="385px" width="900px">
</a>
</p>

```r
> highlight <- map[pve>0.001,1]
> CMplot(cbind(map,nonZeroRate), type="h", plot.type="m", LOG10=FALSE, ylab="Phenotypic variance explained (%)",
        highlight=highlight, highlight.col=NULL)
```
<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/hibayes/master/figure/2.jpg">
<img src="figure/2.jpg" height="385px" width="900px">
</a>
</p>

#### (2) Gemone-Wide association study
**WPPA** is defined to be the window posterior probability of association ([Fernando and Garrick (2013)](https://link.springer.com/protocol/10.1007/978-1-62703-447-0_10)), it is the ratio of the number of iterations that ***Pw*** (the proportion of the total genetic variance explained by the window ***w***) > 1% divided by the total number of MCMC iterations.
```r
> fit <- bayes(y=pheno[, 1], X=geno, map=map, windsize=1e6, wppa=0.01, model="BayesR", niter=20000, nburn=10000, outfreq=10)
> gwas <- fit$gwas
> head(gwas)
   WIND CHR NUM   START     END wppa         wgve
1 wind1   1   3 1198554 1825948    0 4.500114e-05
2 wind2   1   1 3428453 3428453    0 8.286156e-06
3 wind3   1   8 4195032 4916148    0 6.749507e-05
4 wind4   1   7 5109162 5881216    0 2.985705e-05
5 wind5   1   3 6705835 6952985    0 2.150887e-05
6 wind6   1   7 7075618 7863025    0 6.133333e-05
```
View the results by [CMplot](https://github.com/YinLiLin/R-CMplot) package:
```r
> highlight <- gwas[gwas$wppa>0.95, 1]
> CMplot(gwas[,c(1,2,4,6)], type="h", plot.type="m", LOG10=FALSE, ylab="WPPA", ylim=c(0,1.2), 
        highlight=highlight, highlight.col=NULL, highlight.text=highlight)
```
<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/hibayes/master/figure/3.jpg">
<img src="figure/3.jpg" height="385px" width="900px">
</a>
</p>

## Not done yet
