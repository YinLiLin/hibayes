
summary.blrMod <- function(object, ...) {

    is_bayes <- grepl("Individual", attr(object$call, "model"))
    is_sbayes <- grepl("Summary", attr(object$call, "model"))
    is_ssbayes <- grepl("Single-step", attr(object$call, "model"))

    res <- list()
    res$call <- object$call

    if(is_bayes || is_ssbayes){
        coef <- matrix(NA, length(object$beta) + length(object[["J"]]) + 1, 2)
        rownames(coef) <- c("(Intercept)", 
            if(is.null(object[["J"]])) NULL else {"J"},
            names(object$beta)
        )
        # colnames(coef) <- c("Estimate", "Std. Error", "t value")
        colnames(coef) <- c("Estimate", "SD")
        coef[, 1] <- c(object$mu, object[["J"]], object$beta)
        coef[, 2] <- c(apply(object$MCMCsamples$mu, 1, sd), 
            if(is.null(object$MCMCsamples[["J"]])) NULL else {apply(object$MCMCsamples[["J"]], 1, sd)},
            if(is.null(object$MCMCsamples$beta)) NULL else {apply(object$MCMCsamples$beta, 1, sd)}
        )
        # coef[, 3] <- coef[, 1] / coef[, 2]
        res$beta <- coef
    }

    envirR <- matrix(NA, length(object[["Vr"]]) + 1, 2)
    rownames(envirR) <- c(names(object[["Vr"]]), "Residual")
    colnames(envirR) <- c("Variance", "SD")
    envirR[, 1] <- c(object[["Vr"]], object$Ve)
    envirR[, 2] <- c(
        if(is.null(object$MCMCsamples[["Vr"]])) NULL else {apply(object$MCMCsamples[["Vr"]], 1, sd)},
        apply(object$MCMCsamples$Ve, 1, sd)
    )
    res$VER <- envirR
    if(!is.null(object[["r"]])){
        res[["r"]] <- object[["r"]]
        res[["r"]]$SD <- apply(object$MCMCsamples[["r"]], 1, sd)
    }

    geneR <- matrix(NA, 1 + 1 + length(object$Veps) + length(object[["pi"]]), 2)

    rownames(geneR) <- c("Vg", "h2", 
        if(is.null(object$Veps)) NULL else {paste0("V", "\U03b5")},
        paste0("\U03c0",1:length(object[["pi"]]))
    )
    colnames(geneR) <- c("Estimate", "SD")
    geneR[, 1] <- c(object$Vg, object$h2, object$Veps, object[["pi"]])
    geneR[, 2] <- c(apply(object$MCMCsamples$Vg, 1, sd), 
                    apply(object$MCMCsamples$h2, 1, sd),
                    if(is.null(object$MCMCsamples$Veps)) NULL else {apply(object$MCMCsamples$Veps, 1, sd)},
                    apply(object$MCMCsamples[["pi"]], 1, sd))
    res$VGR <- geneR

    if(!is.null(object[["g"]])){
        res[["g"]] <- object[["g"]]
        res[["g"]]$SD <- apply(object$MCMCsamples[["g"]], 1, sd)
    }

    res$alpha <- data.frame(Effect = object$alpha, SD = apply(object$MCMCsamples$alpha, 1, sd))

    if(!is.null(object[["e"]]))  res$e <- object[["e"]]

    class(res) <- "summary.blrMod"
    res
}

print.summary.blrMod <- function(x, ...) {
    
    cat(attr(x$call, "model"), "\n")
    cat("Formula:", x$call, "\n")

    if(!is.null(x$e)){
        cat("\nResiduals ($e):\n")
        print(summary(x$e[, 2], digits=5)[-4])
    }

    digits <- max(3, getOption("digits") - 3)

    if(!is.null(x$beta)){
        cat("\nFixed effects ($beta):\n")
        printCoefmat(x$beta, digits=digits)
    }

    cat("\nEnvironmental random effects ($VER, $r):\n")
    printCoefmat(x$VER, digits=digits)
    if(!is.null(x$e)){
        cat("Number of obs:", nrow(x$e))
        if(nrow(x$VER) > 1){
            cat(", group: ")
            cat(paste(rownames(x$VER)[-nrow(x$VER)], attr(x[["r"]], "nlevel"), sep = ", ", collapse="; "))
        }
        cat("\n")
    }

    cat("\nGenetic random effects ($VGR, $g):\n")
    printCoefmat(x$VGR, digits=digits)
    cat("Number of markers:", nrow(x$alpha), ", predicted individuals:", ifelse(is.null(x[["g"]]), 0, nrow(x[["g"]])), "\n")

    cat("\nMarker effects ($alpha):\n")
    print(summary(x$alpha[, 1], digits=6)[-4])

    invisible(x)
}
