library(dplyr)
library(parallel)

run_row <- function(r, results, sig2x, blk_sz=20, sig2e, mu=0, mu_x, s_0,
                    res_col, runs, output_file){
  opt <- results[r,]
  
  cat("n", opt$n, "p ", opt$p, " | tNSR ", opt$tnsr, " | rho ", opt$rho, " | k ", opt$k,
      " | beta ", opt$beta, " | structure ", opt$structure, " | family ", opt$family,
      " | estimator ", opt$estimator, " | est fun ", opt$est_fun, " Lasso",
      file=output_file, append=TRUE)
  
  # Set more simulation parameters
  sig2u <- opt$tnsr / (1-opt$tnsr) * sig2x
  Sig_xx <- choose_matrix(sig2x, opt$rho, opt$p, blk_sz=blk_sz, structure=opt$structure)
  Sig_uu <- choose_matrix(opt$k*sig2u, 0, opt$p, structure="diag")
  
  # Set active set indices to be equidistant between 1 and p, excluding 1 and p
  # themselves as indices.
  S_0 <- floor(seq(from=1, to=opt$p, length.out=s_0+2))[-c(1,s_0+2)]
  
  # Run simulations for this simulation option.
  res <- run_sim(n=opt$n, p=opt$p, runs=runs, Sig_xx=Sig_xx, blk_sz, Sig_uu=Sig_uu,
                 sig2e=sig2e, mu=mu, mu_x=mu_x, s_0=s_0, k=opt$k,
                 family=as.character(opt$family),
                 beta_rand=FALSE, S_0=S_0, beta_fixed=opt$beta, opt_kappa=opt$kappa,
                 estimator=opt$estimator, est_fun=opt$est_fun,
                 nlambda=100, of=output_file)
  
  opt[(res_col+1):(res_col+length(res))] <- res
  
  cat("\nresults updated\n", file=output_file, append=TRUE)
  
  return ( opt )
}

# Set this to the appropriate working directory on your computer, where all
# downloaded supporting code is stored.
setwd("~/supporting_code")

# Set save directory and file for printing progress updates.
save_dir <- getwd()
output_file <- "output.txt"
cat("", file=output_file)

# Source additional functions from the run_sim.R file.
source("run_sim.R")

# Set number of simulation runs.
runs <- 500

# Options contains the different simulation settings for which we will run
# our simulation study.
opt_grid <- read.csv("options.csv", header=TRUE)

# We want to run the simulations in options.csv for all estimators of interest.
# Expand the options to include an option for the estimator and the estimation
# function.
estimators <- data.frame(estimator=c("True", "Naive", "RC", "CS"))
est_funs <- data.frame(est_fun=c("known", "glasso", "tapering"))

results <- merge(opt_grid, merge(estimators, est_funs), sort=FALSE)
# remove True, Naive, CS Lasso with anything other than "known", since the
# other options are only used for the RC Lasso
results <- filter(
  results,
  (estimator=="RC" | est_fun=="known")
)
results <- arrange(results, structure, n, tnsr, family, estimator, est_fun)

res_len <- nrow(results)
res_col <- ncol(results)

# Set more simulation parameters
sig2x <- 1
sig2e <- 0.1

mu_x <- 0
mu <- 5

s_0 <- 10

# Run simulations in parallel between all options
results <- mclapply(1:res_len, run_row,
                    results=results,
                    sig2x=sig2x, blk_sz=20, sig2e=sig2e, mu=mu, mu_x=mu_x, s_0=s_0,
                    res_col=res_col, runs=runs, output_file=output_file,
                    mc.cores=60)
# Save results as a csv file
results <- matrix(unlist(results), ncol = length(results[[1]]), byrow = TRUE)
write.csv(results, file=paste(save_dir, "results.csv", sep="/"))