summary.blrMod <- function(object, ...) {

    is_bayes <- grepl("Individual", attr(object$call, "model"))
    is_sbayes <- grepl("Summary", attr(object$call, "model"))
    is_ssbayes <- grepl("Single-step", attr(object$call, "model"))

    res <- list()
    res$call <- object$call

    if(is_bayes || is_ssbayes){
        fixed <- matrix(NA, length(object$beta) + length(object$J) + 1, 2)
        rownames(fixed) <- c("(Intercept)", ifelse(is.null(object$J), NULL, "J"), names(object$beta))
        # colnames(fixed) <- c("Estimate", "Std. Error", "t value")
        colnames(fixed) <- c("Estimate", "SD")
        fixed[, 1] <- c(object$mu, object$J, object$beta)
        fixed[, 2] <- c(apply(object$MCMCsamples$mu, 1, sd), 
            ifelse(is.null(object$MCMCsamples$J), NULL, apply(object$MCMCsamples$J, 1, sd)),
            ifelse(is.null(object$MCMCsamples$beta), NULL, apply(object$MCMCsamples$beta, 1, sd))
        )
        # fixed[, 3] <- fixed[, 1] / fixed[, 2]
        res$fixed <- fixed
    }

    envirR <- matrix(NA, length(object$Vr) + 1, 2)
    rownames(envirR) <- c(names(object$Vr), "Residual")
    colnames(envirR) <- c("Variance", "SD")
    envirR[, 1] <- c(object$Vr, object$Ve)
    envirR[, 2] <- c(
        ifelse(is.null(object$MCMCsamples$Vr), NULL, apply(object$MCMCsamples$Vr, 1, sd)),
        apply(object$MCMCsamples$Ve, 1, sd)
    )
    res$envirR <- envirR
    if(!is.null(object$r)){
        res$r <- object$r
        res$r$SD <- apply(object$MCMCsamples$r, 1, sd)
    }

    geneR <- matrix(NA, 1 + 1 + length(object$pi), 2)
    rownames(geneR) <- c("Vg", "h2", ifelse(is.null(object$Veps), NULL, "Vε"), paste0("Pi",1:length(object$pi)))
    colnames(geneR) <- c("Estimate", "Std. Error")
    geneR[, 1] <- c(object$Vg, object$h2, object$Veps, object$pi)
    geneR[, 2] <- c(apply(object$MCMCsamples$Vg, 1, sd), 
                    apply(object$MCMCsamples$h2, 1, sd),
                    ifelse(is.null(object$MCMCsamples$Veps), NULL, apply(object$MCMCsamples$Veps, 1, sd)),
                    apply(object$MCMCsamples$p1, 1, sd))
    res$geneR <- geneR

    if(!is.null(object$g)){
        res$g <- object$g
        res$g$SD <- apply(object$MCMCsamples$g, 1, sd)
    }

    res$alpha <- data.frame(object$alpha, SD = apply(object$MCMCsamples$alpha, 1, sd))

    if(!is.null(object$e))  res$residuals <- object$e

    class(res) <- "summary.hibayes"
    res
}

print.summary.hibayes <- function(x, ...) {
    
    cat(attr(x$call, "model"), "\n")
    cat("Formula:", x$call, "\n")

    if(!is.null(x$residuals)){
        cat("\nResiduals:\n")
        print(summary(x$residuals[, 2], digits=5)[-4])
    }

    digits <- max(3, getOption("digits") - 3)

    if(!is.null(x$fixed)){
        cat("\nFixed effects:\n")
        printCoefmat(x$fixed, digits=digits)
    }

    cat("\nEnvironmental random effects:\n")
    printCoefmat(x$envirR, digits=digits)
    if(!is.null(x$residuals)){
        cat("Number of obs:", nrow(x$residuals))
        if(nrow(x$envirR) > 1){
            cat(", group: ")
            cat(paste(rownames(x$envirR)[-nrow(envirR)], attr(x$r, "nlevel"), sep = ", ", collapse="; "))
            cat("\n")
        }
    }

    cat("\nGenetic random effects:\n")
    printCoefmat(x$geneR, digits=digits)
    cat("Number of markers:", nrow(x$residuals), "predicted Individuals:", nrow(x$g))

    cat("\nMarker effects:\n")
    print(summary(x$alpha[, 2], digits=6)[-4])

    invisible(x)
}
