library(dplyr)
library(parallel)

run_row <- function(r, results, q, sig2x, sig2z, blk_sz=20, sig2e, Sig_xz,
                    mu=0, mu_x, mu_z, s_0,
                    res_col, runs, output_file){
  opt <- results[r,]
  
  cat("n", opt$n, "p ", opt$p, " | tNSR ", opt$tnsr, " | rho ", opt$rho, " | k ", opt$k,
      " | beta ", opt$beta, " | structure ", opt$structure, " | family ", opt$family,
      " | estimator ", opt$estimator, " | est fun ", opt$est_fun, " Lasso",
      file=output_file, append=TRUE)
  
  sig2u <- opt$tnsr / (1-opt$tnsr) * sig2x
  Sig_xx <- choose_matrix(sig2x, opt$rho, opt$p, blk_sz=blk_sz, structure=opt$structure)
  Sig_zz <- choose_matrix(sig2z, p=q, structure="diag")
  Sig_uu <- choose_matrix(opt$k*sig2u, 0, opt$p, structure="diag")
  
  S_0 <- floor(seq(from=1, to=opt$p, length.out=s_0+2))[-c(1,s_0+2)]
  
  res <- run_sim(n=opt$n, p=opt$p, q=q, runs=runs, Sig_xx=Sig_xx, blk_sz, Sig_uu=Sig_uu,
                 sig2e=sig2e, mu=mu, mu_x=mu_x, mu_z=mu_z, s_0=s_0, k=opt$k,
                 Sig_zz=Sig_zz, Sig_xz=Sig_xz,
                 family=as.character(opt$family),
                 beta_rand=FALSE, S_0=S_0, beta_fixed=opt$beta,
                 opt_kappa=opt$kappa, estimator=opt$estimator, est_fun=opt$est_fun,
                 nlambda=100, of=output_file)
  
  opt[(res_col+1):(res_col+length(res))] <- res
  
  cat("\nresults updated\n", file=output_file, append=TRUE)
  
  return ( opt )
}

save_dir <- "~/supporting_code"
output_file <- "output.txt"
cat("", file=output_file)

source("run_sim_mu_Z.R")

runs <- 500

opt_grid <- read.csv("options_mu_Z.csv", header=TRUE)

estimators <- data.frame(estimator=c("True", "Naive", "RC"))
est_funs <- data.frame(est_fun=c("known", "glasso"))

results <- merge(opt_grid, merge(estimators, est_funs), sort=FALSE)
# remove True, Naive, Corrected with anything other than "known"
results <- filter(
  results,
  (estimator=="RC" | est_fun=="known")
)
results <- arrange(results, structure, n, tnsr, family, estimator, est_fun)

res_len <- nrow(results)
res_col <- ncol(results)

sig2x <- 1
sig2z <- 1
sig2e <- 0.1

q <- 2
mu_x <- 5
mu_z <- 0
mu <- 100

s_0 <- 10

Sig_xz <- matrix(0, nrow=results$p[1], ncol=q)
Sig_xz[cbind(sample(1:results$p[1], q), 1:q)] <- 0.4

results <- mclapply(1:res_len, run_row,
                    results=results,
                    q=q,
                    sig2x=sig2x, sig2z=sig2z, blk_sz=20, sig2e=sig2e,
                    Sig_xz=Sig_xz, mu=mu, mu_x=mu_x, mu_z=mu_z,
                    s_0=s_0, res_col=res_col, runs=runs, output_file=output_file,
                    mc.cores=60)
results <- matrix(unlist(results), ncol = length(results[[1]]), byrow = TRUE)

write.csv(results, file=paste(save_dir, "results_mu_Z.csv", sep="/"))