library(tidyverse)
library(ggpubr)

select <- dplyr::select

# Reproduce parameters from simulation studies
nlambda <- 100
runs <- 500
s_0 <- 10
p <- 1000

# Set this to the appropriate working directory on your computer, where all
# downloaded supporting code is stored.
wd <- "~/supporting_code"
setwd(wd)

# Give names to the first few columns. Needed for further parsing.
col_nms <- c("X", "n",	"p",	"tnsr"	, "rho", "beta",	"family",	"k",
             "structure",	"kappa",	"estimator",	"est_fun")

res1 <- read.csv("./results.csv", header=TRUE)
colnames(res1)[1:length(col_nms)] <- col_nms
# Rename options for cleanliness of plots
res1$structure <- ifelse(res1$structure=="diag", "Diagonal",
                         ifelse(res1$structure=="block", "Block", "NA"))
res1$family <- ifelse(res1$family=="gaussian", "Linear",
                      ifelse(res1$family=="binomial", "Logistic", "NA"))

# Remove True Lasso for the plots
res <- res1[res1$estimator != "True" & res1$estimator != "TRUE",]
# Reorganize results to long format
roc1 <- data.frame(ind=rep(1:nlambda, 2), Metric=rep(c("TP","FP"), each=nlambda),
                   t(as.matrix( res[, -(1:12)][, 1:(2*nlambda)] )))
# Rename columns
colnames(roc1)[-c(1:2)] <- apply(
  distinct(res, family, structure, estimator, est_fun),
  MARGIN=1, paste, collapse="_")
# Reshape again to long format based on linear/logistic setting, covariance
# structure, estimator, and estimator function.
roc2 <- roc1 %>% gather(key="estimator", value="Rate", -c(ind, Metric)) %>%
  reshape(idvar=c("ind", "estimator"), timevar="Metric",
          direction="wide")
# Separate varying options into separate fields
roc <- data.frame(separate_wider_delim(roc2, cols=estimator, delim="_",
                                       names=colnames(res)[c(7,9,11,12)]))
# Convert TP/FP incidence into TP/FP rates
roc$Rate.TP <- roc$Rate.TP / (s_0)
roc$Rate.FP <- roc$Rate.FP / (p-s_0)

# Plot results as subplots (Figure 2)
grp_names <- c("structure", "family")
p_roc <- ggplot(data=roc, aes(x=Rate.FP, y=Rate.TP,
                                         color=interaction(estimator, est_fun))) +
  geom_line() +
  geom_point(size=0.5) +
  theme(panel.background=element_rect(fill="white", color="white"),
        # remove the vertical grid lines
        panel.grid.major.x = element_line( linewidth=.1, color="white"),
        panel.grid.minor.x = element_line( linewidth=.1, color="white"),
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( linewidth=.1, color="black"),
        panel.grid.minor.y = element_line( linewidth=.1, color="black"),
        legend.position="top",
        legend.title=element_blank()
  )  +
  scale_color_manual(labels=c("RC.glasso" = "RC (glasso)",
                              "CS.known"="CS",
                              "Naive.known"="Naive",
                              "RC.known"="RC (known)",
                              "RC.tapering"="RC (tapering)"),
                     values=c("red", "orange", "green", "blue", "purple")) +
  xlim(0, 20/(p-s_0)) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  facet_wrap(facets=grp_names, scales="free")

################################## Tables for l1 error #########################
# The results below are used to generate Table 1.
# Extract l1 error results
t_l1err <- cbind(
  res1[, 1:12],
  res1[, -(1:12)][, (2*nlambda+1):(2*nlambda + runs)]
)
# Extract l1 error results in active set
t_l1errS0 <- cbind(
  res1[, 1:12],
  res1[, -(1:12)][, (2*nlambda+runs+1):(2*nlambda + 2*runs)]
)
# Obtain average and sd estimates of l1 and l1(S) results
l1err <- cbind( t_l1err[,1:12], mean=rowMeans(t_l1err[,-c(1:12)]),
                sd=apply(t_l1err[,-c(1:12)], 1, sd) )
l1err$mean <- round(l1err$mean, 3)
l1err$sd <- round(l1err$sd, 3)
l1err <- l1err %>% arrange(n, p, k, family, structure, estimator, est_fun)

l1errS0 <- cbind( t_l1errS0[,1:12], mean=rowMeans(t_l1errS0[,-c(1:12)]),
                  sd=apply(t_l1errS0[,-c(1:12)], 1, sd) )
l1errS0$mean <- round(l1errS0$mean, 3)
l1errS0$sd <- round(l1errS0$sd, 3)
l1errS0 <- l1errS0 %>% arrange(n, p, k, family, structure, estimator, est_fun)