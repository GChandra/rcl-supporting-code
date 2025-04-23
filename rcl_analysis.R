library(tidyverse)
library(rcl)
library(glasso)
library(glmnet)

set.seed(4012)

select <- dplyr::select

# define glasso estimator function for the rcl.
inv_fun <- function(M, rho, ...){
  return( glasso::glasso(M, rho=rho)$wi )
}

setwd("~/supporting_code")
y <- read.csv("asa24_tn.d071919.csv")
dat <- readRDS("preprocessed_dat.rds") %>% arrange(SUBJECT_ID)

# metabolites with estimated measurement error variance > 2, whch are filtered
# out
metabs_noisy <- c(21, 58, 63, 76, 96, 116, 183, 190, 195, 203, 234, 247, 249,
                  274, 348, 362, 371, 395)
dat <- dat[, setdiff(1:ncol(dat), metabs_noisy+2)]

merge_dat <- merge(dat, y, by.x="SUBJECT_ID", by.y="iid")

# nutrients to include
nutrients_inc <- c("tn_acar_asa24", "tn_alc_asa24", "tn_atoc_asa24", "tn_b12_add_asa24",
                   "tn_bcar_asa24", "tn_caff_asa24", "tn_calc_asa24", "tn_carb_asa24",
                   "tn_chole_asa24", "tn_choln_asa24", "tn_copp_asa24", "tn_cryp_asa24",
                   "tn_fdfe_asa24", "tn_fibe_asa24", "tn_iron_asa24", "tn_lyco_asa24",
                   "tn_lz_asa24", "tn_magn_asa24", "tn_mfat_asa24", "tn_mois_asa24",
                   "tn_niac_asa24", "tn_pfat_asa24", "tn_phos_asa24", "tn_pota_asa24",
                   "tn_prot_asa24", "tn_sele_asa24", "tn_sfat_asa24", "tn_sodi_asa24",
                   "tn_sugr_asa24", "tn_theo_asa24", "tn_vara_asa24", "tn_vb12_asa24",
                   "tn_vb1_asa24", "tn_vb2_asa24", "tn_vb6_asa24", "tn_vc_asa24",
                   "tn_vitd_asa24", "tn_vite_add_asa24", "tn_vk_asa24", "tn_zinc_asa24")
inc_col <- c("SUBJECT_ID", "task", nutrients_inc)

# average 24HR over month 0 to 6
y <- merge_dat %>% filter(VISIT==1) %>% select(all_of(inc_col)) %>%
  filter(task <= 3) %>%
  arrange(SUBJECT_ID) %>%
  group_by(SUBJECT_ID) %>%
  summarise(across( .cols=starts_with("tn_"),
                    .fns=list(mean=mean),
                    .names="{.col}")
  )

# remove subjects without both metabolite data and 24HR data.
dat <- dat %>% filter(SUBJECT_ID %in% y$SUBJECT_ID) %>%
  arrange(SUBJECT_ID)

y <- select(y, -SUBJECT_ID)

W <- split(select(dat, -c(SUBJECT_ID, VISIT)), dat$SUBJECT_ID)

# run data analysis with glasso with rho=0.2 (as with the simulation study).
rho <- 0.2
n <- length(W)
p <- ncol(W[[1]])
k <- rep(2, times=length(unique(dat$SUBJECT_ID)))

# For subjects with only 1 replicate, fill in the other with the same value as
# the first replicate. This allows us to easily apply the RC method without
# affecting the results.
for (i in 1:n){
  for (j in 1:p){
    W[[i]][which(is.na(W[[i]][,j])), j] <- mean(W[[i]][,j], na.rm=TRUE)
  }
}

# average of replicates
W_bar <- t(sapply(W, colMeans, na.rm=TRUE, simplify=TRUE))

# estimate Sig_u, a diagonal matrix
SCov_u <- lapply(1:n,
                 function(i) cov(W[[i]]) )
SCov_u <- diag(diag( Reduce('+', SCov_u) / (sum(k) - n) ))

# estimate other quantities
mu_w_hat <- colSums( diag(k) %*% W_bar ) / sum(k)

S_w <- glasso::glasso( cov(W_bar), rho=rho )
SCov_w <- S_w$w
SCov_x <- SCov_w - SCov_u/2

Lambda_hat <- S_w$wi %*% SCov_x

y0 <- c(as.matrix(y[,"tn_caff_asa24"]))

# Run Naive, RC, and CS Lasso and obtain coefficient estimates for each
me_list <- list(
  k=k,
  y=y0,
  W=W
)

cv.nl.fit <- glmnet::cv.glmnet(W_bar, y0, family="gaussian",
                               intercept=TRUE, standardize=TRUE, nfolds=10)
betaNL_hat <- as.matrix(coef(cv.nl.fit, s=cv.nl.fit$lambda.1se))[2:(p+1)]

cv.rcl.fit <- rcl::cv_rcl(me_list, Lambda=Lambda_hat,
                          family="gaussian", intercept=TRUE, standardize=TRUE, nfolds=10)
betaRC_hat <- as.matrix(coef(cv.rcl.fit, s=cv.rcl.fit$lambda.1se))[2:(p+1)]

betaNL_1 <- as.matrix(coef(cv.nl.fit, s=cv.nl.fit$lambda.1se))[2:(p+1)]
R <- 2*sum(abs(betaNL_1))
radii <- seq(from=1e-3*R, to=R, length.out=20)
cv.csl.fit <- hdme::cv_corrected_lasso(W_bar, y0, SCov_u/2,
                                       radii=radii)
csl.opt.fit <- hdme::corrected_lasso(W_bar, y0, SCov_u/2,
                                     family="gaussian", radii=cv.csl.fit$radius_1se)
betaCS_hat <- csl.opt.fit$betaCorr

# Store results in a data frame, which is used to make Table 2.
res <- data.frame(metabolite=colnames(W_bar),
                  NL=betaNL_hat, RC=betaRC_hat, CS=betaCS_hat) %>%
  arrange(desc(abs(NL)), desc(abs(RC)), desc(abs(CS)))