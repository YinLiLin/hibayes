library(randRotation)

X <- randorth(1000)
h2 <- 0.3

qtn.idx <- sample.int(1000, size=50)
effects <- rnorm(50) * 10

bv = X[, qtn.idx] %*% effects
y = bv + rnorm(1000, sqrt(var(bv) * (1 - h2) / h2))

X_train = X[1:800, ]
X_test  = X[801:1000, ]

y_train = y[1:800]
y_test  = y[801:1000]

library(hibayes)

fit <- hibayes:::VariationalBayes(
    y = y_train,
    X = X_train,
    model_str = "BayesC",
    Pi = c(0.95, 0.05),
    max_iteration = 500,
    threads = 1,
    block_size = 1,
)

y_hat_vi <- X_test %*% fit$beta


# ==================================================
# pheno <- data.frame(ID=paste0("ID", 1:800),Y=y_train)

# fit <- ibrm(
#     as.formula('Y ~ 1'),
#     data = pheno, 
#     M = X_train,
#     M.id = pheno$ID,
#     method = "BayesC"
#     )

# y_hat_mcmc <- X_test %*% fit$alpha

# ==================================================
colnames(X_train) <- paste0("V", 1:ncol(X_train))

alpha <- c()
for(i in 1:ncol(X_train)) {
    alpha <- rbind(alpha, lm(y_train ~ X_train[, i])$coefficients[2])
}

print(paste0("VI: ", cor(y_hat_vi, y_test)))
print(paste0("MCMC: ", cor(y_hat_mcmc, y_test)))
print(paste0("LM: ", cor(X_test %*% alpha, y_test)))