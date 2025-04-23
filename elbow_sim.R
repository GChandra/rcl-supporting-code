library(dplyr)
library(parallel)

run_row <- function(r, results, sig2x, blk_sz=20, sig2e, mu=0, s_0, d_grid,
                    res_col, runs, output_file){
  opt <- results[r,]
  
  cat("n", opt$n, "p ", opt$p, " | tNSR ", opt$tnsr, " | rho ", opt$rho, " | k ", opt$k,
      " | beta ", opt$beta, " | structure ", opt$structure, " | family ", opt$family,
      " | estimator ", opt$estimator, " Lasso",
      file=output_file, append=TRUE)
  
  sig2u <- opt$tnsr / (1-opt$tnsr) * sig2x
  Sig_xx <- choose_matrix(sig2x, opt$rho, opt$p, blk_sz=20, structure=opt$structure)
  Sig_uu <- choose_matrix(opt$k*sig2u, 0, opt$p, structure="diag")
  
  NZ <- get_elbow(n=opt$n, p=opt$p, runs=runs, Sig_xx=Sig_xx, Sig_uu=Sig_uu,
                  sig2e=sig2e, s_0=s_0, k=opt$k, mu=mu,
                  family=as.character(opt$family),
                  d_grid=d_grid, beta_rand=FALSE, beta_fixed=opt$beta,
                  nlambda=100, of=output_file)
  
  filename <- paste0(
    paste("elbow", opt$n, opt$p, opt$k,
          strsplit(as.character(opt$tnsr), split="\\.")[[1]][2], opt$structure,
          opt$family, sep="_"), ".png"
  )
  
  return (list(
    filename=filename,
    NZ=NZ
  ))
}

setwd("~/supporting_code")

source("run_sim.R")

save_dir <- getwd()
output_file <- "elbow_output.txt"
cat("", file=output_file)

opt_grid <- read.csv("elbow_options.csv", header=TRUE)

estimators <- "CS"
results <- opt_grid[rep(1:nrow(opt_grid), each=length(estimators)), ]
results["estimator"] <- c(rep(estimators, times=nrow(opt_grid)))
res_len <- nrow(results)
res_col <- ncol(results)

sig2x <- 1
sig2e <- 0.1

mu_x <- 0
mu <- 5

s_0 <- 10

runs <- 100

d_grid <- seq(from=0.1, to=5, by=0.1)

res <- mclapply(1:res_len, run_row,
                results=results,
                sig2x=sig2x, blk_sz=20, sig2e=sig2e, mu=mu, s_0=s_0,
                d_grid=d_grid, res_col=res_col, runs=runs,
                output_file=output_file, mc.cores=20)

res_csv <- matrix( unlist(lapply(res, function(x) x$NZ)),
                   nrow=length(res), ncol=length(d_grid), byrow=TRUE)
rownames(res_csv) <- unlist(lapply(res, function(x) x$filename))
write.csv(res_csv, file=paste(save_dir, "results_elbow.csv", sep="/"))

for (r in 1:length(res)){
  png(file=paste(save_dir, res[[r]]$filename, sep="/"),
      width=500, height=500)
  plot(d_grid, res[[r]]$NZ, xlim=rev(range(d_grid)), xlab="NZ")
  dev.off()
  cat("\nelbow plot saved\n", file=output_file, append=TRUE)
}
