#!/usr/bin/env Rscript

DATASETS_ROOT <- '/Users/hyacz/Data/GS2019/'

args = commandArgs(trailingOnly=TRUE)

library(hibayes)
library(VIGoR)

runVIGoR <- function(method, datasets, trait, fold) {
    bin_path <- paste0(DATASETS_ROOT, datasets, "/", datasets, ".", trait, ".qc")
    print(paste0("reading ", bin_path, " ..."))
    bin <- read_plink(bin_path)

    fam  <- bin[["fam"]]
    geno <- bin[["geno"]]
    map  <- bin[["map"]]
    X <- geno[]

    #  BL     = matrix(c(1, 1), nrow=1),
    #  EBL    = matrix(c(0.1, 0.1, 1, 0.1), nrow=1),
    #  BayesA = matrix(c(5, 0.01, 1), nrow=1),
    #  BayesB = matrix(c(5, 0.1, 0.01), nrow=1),
    #  BayesC = matrix(c(5, 0.1, 0.01), nrow=1),
    #  BRR    = matrix(c(5, 0.01, 1), nrow=1),
    #  BLUP   = matrix(c(5, 0.3), nrow=1))

    H <- NULL
    if (method == "BayesA") {
        H <- hyperpara(X, 0.5, "BayesA")
    } else if (method == "BayesB") {
        H <- hyperpara(X, 0.5, "BayesB", 0.05)
    } else if (method == "BayesC") {
        H <- hyperpara(X, 0.5, "BayesC", 0.05)
    }
    ETA <- list(list(model = method, X = X, H = H))
    
    for (f in fold) {
        Y   <- fam[, 6 + f]
        print(paste("Running: ", method, datasets, trait, f))
        time1 <- Sys.time()
        Result <- vigor(Y, ETA)
        time2 <- Sys.time()
        print(paste("=== FINISH === "))
        # calc Accuracy
        test.id <- is.na(Y)
        ebv <- geno[test.id, ] %*% Result$ETA[[1]]$Beta
        Result$pcc <- cor(fam[test.id, 6], ebv)
        Result$time <- time2 - time1
        saveRDS(Result, file = paste0("VIGoR/", paste(method, datasets, trait, f, sep='_'), ".rds"))
    }
}

# test if there is at least one argument: if not, return an error
if (length(args) == 3) {
    print(args)
    method   <- args[1]
    dataset  <- args[2]
    trait    <- args[3]
    runVIGoR(method, dataset, trait, c(1:20))

} else {
    stop("runVIGoR.R <method> <datasets> <trait> or runVIGoR.R", call.=FALSE)
}

