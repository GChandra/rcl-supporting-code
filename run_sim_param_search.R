library(MASS)
library(Matrix)
library(glmnet)
library(hdme)
library(matrixcalc)
library(glasso)

run_sim <- function(n, p, runs, Sig_xx, Sig_uu, sig2e, mu=0, mu_x, s_0, k, family,
                    beta_rand=FALSE, beta_fixed=NULL, opt_kappa, estimator, nlambda, of){
  if (length(mu_x) == 1){
    mu_x <- rep(mu_x, p)
  }

  # set fields to store results
  tp <- fp <- rep(0, nlambda)
  l1err <- l1err_S0 <- rep(NA, runs)
  
  # keep track of time
  runs_start <- proc.time()
  upd_start <- runs_start
  
  for (r in 1:runs){
    # print time every 25 runs
    if (r %% 25 == 0){
      upd_time <- proc.time()
      cat("\nprogress update:\nrun ", r, "\n",
              round((upd_time-upd_start)[3], 2), " seconds since last update",
          file=of, append=TRUE)
      upd_start <- proc.time()
    }
    
    # generate active set and true coefficient vector
    S_0 <- sample(1:p, s_0)
    beta_0 <- rep(0, p)
    if (beta_rand){
      if (family == "gaussian"){
        beta_0[S_0] <- rnorm(s_0, mean=0, sd=2)
      } else {
        beta_0[S_0] <- rnorm(s_0, mean=0, sd=5)
      }
    } else {
      if (is.null(beta_fixed)){
        if (family == "gaussian"){
          beta_0[S_0] <- 1
        } else {
          beta_0[S_0] <- sample(c(-0.5, 0.5), size=s_0, replace=TRUE)
        }
      } else {
        if (family == "gaussian"){
          beta_0[S_0] <- beta_fixed
        } else {
          beta_0[S_0] <- sample(c(-beta_fixed, beta_fixed), size=s_0, replace=TRUE)
        }
      }
    }
    
    # generate true predictors
    X <- mvrnorm(n, mu=mu_x, Sigma=Sig_xx)
    
    # generate response
    if (family == "binomial"){
      mu <- -mean(X%*%beta_0)
      
      z <- mu + X %*% beta_0
      pr <- 1/(1+exp(-z))
      y <- rbinom(n, 1, pr)
    } else if (family == "gaussian"){
      e <- rnorm(n, 0, sqrt(sig2e))
      y <- mu + X %*% beta_0 + e
      y0 <- y - mean(y)
    }
    
    if (k==1){
      # generate error-prone predictors and compute reliability matrix.
      # if k = 1, we assume the matrices to estimate are known.
      U <- mvrnorm(n, mu=rep(0,p), Sigma=Sig_uu)
      W_bar <- X + U
      W_bar0 <- scale(W_bar, center=TRUE, scale=FALSE)
      
      SCov_uu <- Sig_uu
      Lambda <- solve(Sig_xx + Sig_uu) %*% Sig_xx
      X_hat <- matrix( rep( t(colMeans(W_bar)) %*% (diag(p) - Lambda), n ),
                       byrow=TRUE,
                       nrow=n ) +
        W_bar %*% Lambda
    } else {
      # generate error-prone predictors, store replicates in a list. Get average
      # of replicates as well.
      W_list <- list()
      W_bar <- matrix(0, n, p)
      for (i in 1:k){
        W_list[[i]] <- X + mvrnorm(n, mu=rep(0,p), Sigma=Sig_uu)
        W_bar <- W_bar + W_list[[i]]
      }
      W_bar <- W_bar / k
      
      # estimate measurement error variances (diagonal matrix)
      svar2u <- rep(0, p)
      for (i in 1:p){
        for (j in 1:k){
          svar2u[i] <- svar2u[i] + sum( (W_list[[j]][,i] - W_bar[,i])^2 )
        }
      }
      svar2u <- svar2u / (n*(k-1))
      SCov_uu <- diag(svar2u)
      
      W_bar0 <- scale(W_bar, center=TRUE, scale=FALSE)
      
      if (estimator == "RC"){
        # get glasso estimate of error-prone data covariance matrix.
        gl <- glasso(cov(W_bar0), rho=0.2)
        Lambda <- gl$wi %*% (gl$w - SCov_uu/k)
        # Compute reliability matrix
        X_hat <- matrix( rep( t(colMeans(W_bar)) %*% (diag(p) - Lambda), n ),
                         byrow=TRUE,
                         nrow=n ) +
          W_bar %*% Lambda
      }
    }
    
    X0 <- scale(X, center=TRUE, scale=FALSE)
    
    if (estimator == "True"){
      # True LASSO (for comparison)
      if (family == "gaussian") {
        l.min <- 0.1
        lambda_max <- 10
      } else {
        l.min <- 0.01
        lambda_max <- 1
      }
      lambda_seq <- exp( seq( log(lambda_max), log(l.min),
                              length.out=nlambda ) )
      # perform cross-validated Lasso over default regularization parameters
      cvfit <- cv.glmnet(X, y, intercept=TRUE, standardize=TRUE,
                         family=family, nfolds=10)
      # use 1-se lambda for linear, min lambda for logistic
      if (family == "gaussian"){
        mu_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[1]
        beta_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[2:(p+1)]
      } else {
        mu_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[1]
        beta_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[2:(p+1)]
      }

      # get l1 errors of coefficient estimates
      l1err_S0[r] <- sum(abs(beta_hat[S_0] - beta_0[S_0]))
      l1err[r] <- sum(abs(beta_hat - beta_0))
      
      # fit glmnet over lambda_seq for true/false positive curve.
      fit <- glmnet(X, y, intercept=TRUE, standardize=TRUE,
                       lambda=lambda_seq, family=family)
      for (l in 1:nlambda){
        tp[l] <- tp[l] + sum( coef(fit, s=fit$lambda[l])[2:(p+1)][S_0]!= 0 )
        fp[l] <- fp[l] + sum( coef(fit, s=fit$lambda[l])[2:(p+1)][-S_0] != 0 )
      }
    } else if (estimator == "Naive") {
      # Naive LASSO
      if (family == "gaussian") {
        l.min <- 0.1
        lambda_max <- 10
      } else {
        l.min <- 0.01
        lambda_max <- 1
      }
      lambda_seq <- exp( seq( log(lambda_max), log(l.min),
                              length.out=nlambda ) )
      
      cvfit <- cv.glmnet(W_bar, y, intercept=TRUE, standardize=FALSE,
                         family=family, nfolds=10)
      
      if (family == "gaussian"){
        mu_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[1]
        beta_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[2:(p+1)]
      } else {
        mu_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[1]
        beta_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[2:(p+1)]
      }
      
      l1err_S0[r] <- sum(abs(beta_hat[S_0] - beta_0[S_0]))
      l1err[r] <- sum(abs(beta_hat - beta_0))

      fit <- glmnet(W_bar, y, intercept=TRUE, standardize=TRUE,
                          lambda=lambda_seq, family=family)
      
      for (l in 1:nlambda){
        tp[l] <- tp[l] + sum( coef(fit, s=fit$lambda[l])[2:(p+1)][S_0] != 0 )
        fp[l] <- fp[l] + sum( coef(fit, s=fit$lambda[l])[2:(p+1)][-S_0] != 0 )
      }
    } else if (estimator == "RC") {
      if (family == "gaussian") {
        l.min <- 0.1
        lambda_max <- 10
        } else {
        l.min <- 0.01
        lambda_max <- 1
        }
      lambda_seq <- exp( seq( log(lambda_max), log(l.min),
                              length.out=nlambda ) )
      
      cvfit <- cv.glmnet(X_hat, y, intercept=TRUE, standardize=TRUE,
                         family=family, nfolds=10)
      
      if (family == "gaussian"){
        mu_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[1]
        beta_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[2:(p+1)]
      } else {
        mu_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[1]
        beta_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[2:(p+1)]
      }
      
      l1err_S0[r] <- sum(abs(beta_hat[S_0] - beta_0[S_0]))
      l1err[r] <- sum(abs(beta_hat - beta_0))

      fit <- glmnet(X_hat, y, intercept=TRUE, standardize=TRUE,
                           lambda=lambda_seq, family=family)
      for (l in 1:nlambda){
        tp[l] <- tp[l] + sum( coef(fit, s=fit$lambda[l])[2:(p+1)][S_0] != 0 )
        fp[l] <- fp[l] + sum( coef(fit, s=fit$lambda[l])[2:(p+1)][-S_0] != 0 )
      }
    } else { # Corrected (and CS) Lasso
      kappa <- seq(2e-3, 9, length.out=nlambda) * 2
      
      if (family == "gaussian"){
        cvfit <- tryCatch(
          {
            hdme::cv_corrected_lasso(W_bar, y, SCov_uu/k,
                                     n_folds = 10, tol=1e-4)
          },
          error = function(e){ # nl norm sometimes (rarely) 0, make a small value
            print(e)
            R_grid <- seq(2e-3, 2, length.out=nlambda)*0.01
            hdme::cv_corrected_lasso(W_bar, y, SCov_uu/k, radii=R_grid,
                                     n_folds = 10, tol=1e-4)
          }
        )
        
        ccl_fit <- hdme::corrected_lasso(W_bar, y, SCov_uu/k, family = "gaussian",
                                         radii = cvfit$radius_1se)
        beta_hat <- ccl_fit$betaCorr
      } else {
        ccl_fit <- hdme::corrected_lasso(W_bar, y, SCov_uu/k, family = family,
                                         radii = opt_kappa, tol=1e-4)
        beta_hat <- ccl_fit$betaCorr
      }
      l1err_S0[r] <- sum(abs(beta_hat[S_0] - beta_0[S_0]))
      l1err[r] <- sum(abs(beta_hat - beta_0))

      fit <- hdme::corrected_lasso(W_bar, y, SCov_uu/k, family=family,
                                   radii=kappa, tol=1e-4)
      tp <- tp + colSums( fit$betaCorr[S_0,] != 0 )
      fp <- fp + colSums( fit$betaCorr[-S_0,] != 0 )
    }
  }
  
  tp <- tp / runs
  fp <- fp / runs
  
  runs_time <- proc.time()
  cat("\nruns completed in: ", round((runs_time-runs_start)[3], 2),
          " seconds", file=of, append=TRUE)
  
  res <- c(tp, fp, l1err, l1err_S0)
  res <- round(res, 3)
  
  return (res)
}

# function for generating a covariance matrix based on desired structure.
choose_matrix <- function(sig2, rho, p, blk_sz = 20, structure){
  if(structure == "diag"){
    mat <- matrix(0, p, p)
    diag(mat) <- sig2
  }
  else if(structure == "block"){
    blk <- matrix(NA, blk_sz, blk_sz)
    for (i in 1:blk_sz){
      for (j in 1:blk_sz){
        blk[i,j] <- sig2 * rho^(abs(i-j))
      }
    }
    
    mat <- as.matrix( bdiag(replicate(p/blk_sz, blk, simplify=FALSE)) )
  }
  else{ mat <- matrix(NA, p, p) }
  
  return(mat)
}