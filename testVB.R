library(hibayes)

DATASETS_ROOT <- '/Users/hyacz/Data/GS2019/'
set.seed(1)

runHibayes <- function(method, datasets, trait, fold) {
    fold <- as.numeric(fold)
    bin_path <- paste0(DATASETS_ROOT, datasets, "/", datasets, ".", trait, ".qc")
    print(paste0("reading ", bin_path, " ..."))
    bin <- read_plink(bin_path)

    fam  <- bin[["fam"]]
    geno <- bin[["geno"]]
    map  <- bin[["map"]]

    pheno <- fam
    colnames(pheno) <- paste0("V", 1:length(pheno))
    mask <- is.na(pheno[paste0('V', fold + 6)])

    X <- geno[]

    time1 <- Sys.time()
    print(paste0("Start ..."))
    fit <- hibayes:::VariationalBayes(
        y = pheno[!mask, 6 + fold],
        X = X[!mask, ],
        model_str = method,
        Pi = c(0.95, 0.05),
        max_iteration = 500,
        block_size = 50000,
    )

    time2 <- Sys.time()
    fit$time <- time2 - time1
    
    # mask <- is.na(pheno[paste0('V', fold + 6)])
    
    p <- pheno["V6"][mask]
    g <- X[mask, ] %*% fit$beta

    boxplot(fit$gamma)
    fit$pcc <- cor(p, g)
    print(paste(trait, "Cor:", fit$pcc))

    vigor_res <- readRDS(paste0('VIGoR/', paste("BayesC", datasets, trait, fold, sep='_'), '.rds'))
    plot(fit$beta, vigor_res$ETA[[1]]$Beta)
    abline(0, 1)
    # saveRDS(fit, file = paste0('hibayes/', paste(method, datasets, trait, fold, sep='_'), '.rds'))
    return(str(fit))
}

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) == 4) {
    print(args)
    method   <- args[1]
    dataset  <- args[2]
    trait    <- args[3]
    fold     <- args[4]
    runHibayes(method, dataset, trait, fold)

} else {
    stop("runHibayes.R <method> <datasets> <trait> <fold>", call.=FALSE)
}