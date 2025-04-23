library(MASS)
library(Matrix)
library(glmnet)
library(hdme)
library(glasso)
library(cvCovEst)

run_sim <- function(n, p, q, runs, Sig_xx, blk_sz=20, Sig_uu, sig2e, mu=0, mu_x,
                    mu_z, s_0, k, Sig_zz, Sig_xz, family, beta_rand=FALSE, S_0=NULL,
                    beta_fixed=NULL, opt_kappa, estimator, est_fun, nlambda, of){
  if (length(mu_x) == 1){
    mu_x <- rep(mu_x, p)
  }
  if (length(mu_z) == 1){
    mu_z <- rep(mu_z, q)
  }

  # set matrices to store results
  tp <- fp <- rep(0, nlambda)
  l1err_mu <- rep(NA, runs)
  l1err <- l1err_S0 <- rep(NA, runs)
  l1err_Z <- rep(NA, runs)
  
  runs_start <- proc.time()
  upd_start <- runs_start
  
  if (!beta_rand) {
    beta_0 <- rep(0, p)
    gamma_0 <- rep(0, q)
    if (is.null(beta_fixed)){
      if (family == "gaussian"){
        beta_0[S_0] <- 1
        gamma_0 <- sample(1:3, size=q, replace=TRUE)
      } else {
        beta_0[S_0] <- sample(c(-0.5, 0.5), size=s_0, replace=TRUE)
        gamma_0 <- sample(c(-0.5, 0.5), size=q, replace=TRUE)
      }
    } else {
      if (family == "gaussian"){
        beta_0[S_0] <- beta_fixed
        gamma_0 <- 1:q
      } else {
        beta_0[S_0] <- (S_0 < median(S_0))*1 + (S_0 > median(S_0))*-1
        beta_0[S_0] <- beta_0[S_0]*beta_fixed
        gamma_0 <- rep(-0.5, q)
        gamma_0[1:floor(q/2)] <- 0.5
      }
    }
  }
  
  for (r in 1:runs){
    if (r %% 25 == 0){
      upd_time <- proc.time()
      cat("\nprogress update:\nrun ", r, "\n",
              round((upd_time-upd_start)[3], 2), " seconds since last update",
          file=of, append=TRUE)
      upd_start <- proc.time()
    }
    
    if (beta_rand){
      beta_0 <- rep(0, p)
      S_0 <- sample(1:p, s_0)
      
      if (family == "gaussian"){
        beta_0[S_0] <- rnorm(s_0, mean=0, sd=2)
      } else {
        beta_0[S_0] <- rnorm(s_0, mean=0, sd=5)
      }
      
      gamma_0 <- rnorm(q, mean=0, sd=2)
    }
    
    pred <- MASS::mvrnorm(n, mu=c(mu_x, mu_z),
                          Sigma=cbind(
                            rbind(Sig_xx, t(Sig_xz)),
                            rbind(Sig_xz, Sig_zz)
                          ))
    
    X <- pred[,1:p]
    Z <- pred[,(p+1):(p+q)]
    
    if (family == "binomial"){
      mu <- -c(t(mu_x)%*%beta_0 + t(mu_z)%*%gamma_0)
      
      z <- mu + X %*% beta_0 + Z %*% gamma_0
      pr <- 1/(1+exp(-z))
      y <- rbinom(n, 1, pr)
    } else if (family == "gaussian"){
      e <- rnorm(n, 0, sqrt(sig2e))
      y <- mu + X %*% beta_0 + Z %*% gamma_0 + e
    }
    
    mu_z_hat <- colMeans(Z)
    
    if (k==1){
      U <- mvrnorm(n, mu=rep(0,p), Sigma=Sig_uu)
      W_bar <- X + U
      
      mu_w_hat <- colMeans(W_bar)
      
      SCov_uu <- Sig_uu
      Lambda <- solve(cbind(
        rbind(Sig_xx + Sig_uu, t(Sig_xz)),
        rbind(Sig_xz, Sig_zz)
      )) %*% rbind(Sig_xx, t(Sig_xz))
      
      X_hat <- matrix(rep(1,n), ncol=1) %*%
        ( t(mu_w_hat) - t(c(mu_w_hat, mu_z_hat)) %*% Lambda ) +
        cbind(W_bar, Z) %*% Lambda
    } else {
      W_list <- list()
      W_bar <- matrix(0, n, p)
      for (i in 1:k){
        W_list[[i]] <- X + mvrnorm(n, mu=rep(0,p), Sigma=Sig_uu)
        W_bar <- W_bar + W_list[[i]]
      }
      W_bar <- W_bar / k
      mu_w_hat <- colMeans(W_bar)
      
      svar2u <- rep(0, p)
      for (i in 1:p){
        for (j in 1:k){
          svar2u[i] <- svar2u[i] + sum( (W_list[[j]][,i] - W_bar[,i])^2 )
        }
      }
      svar2u <- svar2u / (n*(k-1))
      SCov_uu <- diag(svar2u)
      
      if (estimator == "RC"){
        if (est_fun == "glasso") {
          gl <- glasso( cov(cbind(W_bar, Z)), rho=0.2 )
          Lambda <- gl$wi %*% rbind(gl$w[1:p,1:p] - SCov_uu/k, gl$w[(p+1):(p+q), 1:p])
        } else if (est_fun == "known") {
          Lambda <- solve(cbind(
            rbind(Sig_xx + SCov_uu/k, t(Sig_xz)),
            rbind(Sig_xz, Sig_zz)
          )) %*% rbind(Sig_xx, t(Sig_xz))
        } else {
          stop("unknown covariance matrix estimator")
        }
        
        X_hat <- matrix(rep(1,n), ncol=1) %*%
          ( matrix(mu_w_hat, nrow=1) - t(c(mu_w_hat, mu_z_hat)) %*% Lambda ) +
          cbind(W_bar, Z) %*% Lambda
      }
    }
    
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
      
      cvfit <- cv.glmnet(cbind(X,Z), y, intercept=TRUE, standardize=TRUE,
                         family=family, nfolds=10)
      if (family == "gaussian"){
        mu_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[1]
        beta_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[2:(p+1)]
        gamma_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[(p+2):(p+q+1)]
      } else {
        mu_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[1]
        beta_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[2:(p+1)]
        gamma_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[(p+2):(p+q+1)]
      }
      
      l1err_mu[r] <- abs(mu_hat - mu)
      l1err_S0[r] <- sum(abs(beta_hat[S_0] - beta_0[S_0]))
      l1err[r] <- sum(abs(beta_hat - beta_0))
      l1err_Z[r] <- sum(abs(gamma_hat - gamma_0))
      
      fit <- glmnet(cbind(X,Z), y, intercept=TRUE, standardize=TRUE,
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
      
      cvfit <- cv.glmnet(cbind(W_bar,Z), y, intercept=TRUE, standardize=FALSE,
                         family=family, nfolds=10)
      
      if (family == "gaussian"){
        mu_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[1]
        beta_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[2:(p+1)]
        gamma_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[(p+2):(p+q+1)]
      } else {
        mu_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[1]
        beta_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[2:(p+1)]
        gamma_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[(p+2):(p+q+1)]
      }
      
      l1err_mu[r] <- abs(mu_hat - mu)
      l1err_S0[r] <- sum(abs(beta_hat[S_0] - beta_0[S_0]))
      l1err[r] <- sum(abs(beta_hat - beta_0))
      l1err_Z[r] <- sum(abs(gamma_hat - gamma_0))
      
      fit <- glmnet(cbind(W_bar,Z), y, intercept=TRUE, standardize=TRUE,
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
      
      cvfit <- cv.glmnet(cbind(X_hat,Z), y, intercept=TRUE, standardize=TRUE,
                         family=family, nfolds=10)
      
      if (family == "gaussian"){
        mu_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[1]
        beta_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[2:(p+1)]
        gamma_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.1se))[(p+2):(p+q+1)]
      } else {
        mu_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[1]
        beta_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[2:(p+1)]
        gamma_hat <- as.matrix(coef(cvfit, s=cvfit$lambda.min))[(p+2):(p+q+1)]
      }
      
      l1err_mu[r] <- abs(mu_hat - mu)
      l1err_S0[r] <- sum(abs(beta_hat[S_0] - beta_0[S_0]))
      l1err[r] <- sum(abs(beta_hat - beta_0))
      l1err_Z[r] <- sum(abs(gamma_hat - gamma_0))
      
      fit <- glmnet(cbind(X_hat,Z), y, intercept=TRUE, standardize=TRUE,
                    lambda=lambda_seq, family=family)
      for (l in 1:nlambda){
        tp[l] <- tp[l] + sum( coef(fit, s=fit$lambda[l])[2:(p+1)][S_0] != 0 )
        fp[l] <- fp[l] + sum( coef(fit, s=fit$lambda[l])[2:(p+1)][-S_0] != 0 )
      }
    } else {
      stop("not a valid estimator")
    }
  }
  
  tp <- tp / runs
  fp <- fp / runs
  
  runs_time <- proc.time()
  cat("\nruns completed in: ", round((runs_time-runs_start)[3], 2),
      " seconds", file=of, append=TRUE)
  
  res <- c(tp, fp, l1err_mu, l1err, l1err_S0, l1err_Z)
  res <- round(res, 3)
  
  return (res)
}

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