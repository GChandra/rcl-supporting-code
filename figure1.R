library(glmnet)
library(MASS)
library(ggplot2)
library(dplyr)

# set seed for reproducibility
set.seed(4239)

# set simulation parameters
n <- 300
p <- 100
s <- 5
S <- sample(1:p, s)

# set correlation for the block-diagonal matrix
rho <- 0.8

# Set this to the appropriate working directory on your computer, where all
# downloaded supporting code is stored.
wd <- "~/supporting_code"
setwd(wd)

# needed for the helper function, choose_matrix, used below
source("run_sim.R")

# set covariance matrices and generate random variables
Sigma_x <- choose_matrix(sig2=1, rho=rho, p=p, structure="block")
Sigma_u <- 0.4*diag(p)
X <- mvrnorm(n=n, mu=rep(0,p), Sigma=Sigma_x)
W <- X + mvrnorm(n=n, mu=rep(0,p), Sigma=Sigma_u)

# set true (sparse) coefficients and generate response
beta <- rep(0, p); beta[S] <- 1
y <- X%*%beta + rnorm(n, sd=sqrt(0.1))

# fit True Lasso and Naive Lasso
fit.TL <- cv.glmnet(X, y, intercept=FALSE, family="gaussian")
fit.NL <- cv.glmnet(W, y, intercept=FALSE, family="gaussian")

# obtain coefficient estimates from True and Naive Lasso, excluding intercept
beta.TL <- coef(fit.TL, s=fit.TL$lambda.1se)[2:(p+1)]
beta.NL <- coef(fit.NL, s=fit.NL$lambda.1se)[2:(p+1)]

# Store results in a data frame, and remove estimates < 0.1 before plotting.
# This is done to keep the figure clean.
res <- data.frame(
  Estimator = rep(c("Naive", "True"), each=p),
  ind = rep(1:p, 2),
  est = c(beta.NL, beta.TL)
) %>% filter(abs(est) >= 0.1)

# Plot simulation results
ggplot(data=res, aes(x=ind, y=est, shape=Estimator, group=Estimator)) +
  geom_point(size=3) +
  geom_hline(aes(yintercept = 1, linetype="True coefficient value")) +
  scale_linetype_manual(name = NULL, values=2) +
  xlim(0, 100) +
  xlab("predictor index") +
  ylab("coefficient estimate") +
  theme_minimal()
