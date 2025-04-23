library(dplyr)

# Set this to the appropriate working directory on your computer, where all
# downloaded supporting code is stored.
save_dir <- "~/supporting_code"
setwd(save_dir)

output_file <- "output.txt"
cat("", file=output_file)

source("run_sim_param_search.R")

runs <- 500

# For tnsr subplot of Figure 3, read 'options_tnsr.csv'. For rho subplot, use
# 'options_rho.csv'. For beta subplot, use 'options_beta.csv'.
results <- read.csv("options_tnsr.csv", header=TRUE)

res_len <- nrow(results)
res_col <- ncol(results)

sig2x <- 1
sig2e <- 0.1

s_0 <- 10

for (r in 1:res_len){
  opt <- results[r,]
  
  cat("n", opt$n, "p ", opt$p, " | tNSR ", opt$tnsr, " | rho ", opt$rho, " | k ", opt$k,
      " | beta ", opt$beta, " | structure ", opt$structure, " | family ", opt$family,
      " | block size ", opt$blk_sz, " | estimator ", opt$estimator, " Lasso",
      file=output_file, append=TRUE)
  
  sig2u <- opt$tnsr / (1-opt$tnsr) * sig2x
  Sig_xx <- choose_matrix(sig2x, opt$rho, opt$p, structure=opt$structure)
  Sig_uu <- choose_matrix(opt$k*sig2u, 0, opt$p, structure="diag")
  
  res <- run_sim(n=opt$n, p=opt$p, runs=runs, Sig_xx=Sig_xx, Sig_uu=Sig_uu,
                 sig2e=sig2e, s_0=s_0, k=opt$k,
                 family=as.character(opt$family),
                 beta_rand=FALSE, beta_fixed=opt$beta, opt_kappa=opt$kappa,
                 estimator=opt$estimator, nlambda=100, of=output_file)
    
  results[r, (res_col+1):(res_col+length(res))] <- res
    
  cat("\nresults updated\n", file=output_file, append=TRUE)
}

# Save results.
# For tnsr subplot of Figure 3, save as 'param_search_tnsr.csv'. For rho subplot,
# save as 'param_search_rho.csv'. For beta subplot, save as
# 'param_search_beta.csv'.
write.csv(results, file=paste(save_dir, "param_search_tnsr.csv", sep="/"))