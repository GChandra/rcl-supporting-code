library(tidyverse)
library(ggpubr)

select <- dplyr::select

wd <- "~/supporting_code"
setwd(wd)

# Give names to the first few columns. Needed for further parsing.
col_nms <- c("X", "n",	"p",	"tnsr"	, "rho", "beta",	"family",	"k",
             "structure",	"kappa",	"estimator",	"est_fun")

res1 <- read.csv("./results_mu_Z.csv", header=TRUE)
colnames(res1)[1:length(col_nms)] <- col_nms
res1$structure <- ifelse(res1$structure=="diag", "Diagonal",
                         ifelse(res1$structure=="block", "Block", "NA"))
res1$family <- ifelse(res1$family=="gaussian", "Linear",
                      ifelse(res1$family=="binomial", "Logistic", "NA"))

res1 <- res1 %>% 
  select(-c(X, p, k, rho, beta, kappa)) %>%
  filter(est_fun != "tapering")

# do not plot True results
res <- res1[res1$estimator != "True" & res1$estimator != "TRUE", ]

nlambda <- 100
runs <- 500
s_0 <- 10
p <- 1000

roc1 <- data.frame(ind=rep(1:nlambda, 2), Metric=rep(c("TP","FP"), each=nlambda),
                   t(as.matrix( res[, -(1:6)][, 1:(2*nlambda)] )))
colnames(roc1)[-c(1:2)] <- apply(
  distinct(res, family, structure, estimator, est_fun),
  MARGIN=1, paste, collapse="_")

roc2 <- roc1 %>% gather(key="estimator", value="Rate", -c(ind, Metric)) %>%
  reshape(idvar=c("ind", "estimator"), timevar="Metric",
          direction="wide")
roc <- data.frame(separate_wider_delim(roc2, cols=estimator, delim="_",
                                       names=colnames(res)[3:6]))
roc$Rate.TP <- roc$Rate.TP / (s_0)
roc$Rate.FP <- roc$Rate.FP / (p-s_0)

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
                              "Naive.known"="Naive",
                              "RC.known"="RC (known)",
                              "RC.tapering"="RC (tapering)"),
                     values=c("red", "green", "blue", "purple")) +
  xlim(0, 10/(p-s_0)) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  facet_wrap(facets=grp_names, scales="free")

################################## Tables for l1 error #########################
t_l1err_mu <- cbind(
  res1[, 1:6],
  res1[, -(1:6)][, (2*nlambda+1):(2*nlambda+runs)]
)
t_l1err <- cbind(
  res1[, 1:6],
  res1[, -(1:6)][, (2*nlambda+runs+1):(2*nlambda+2*runs)]
)
t_l1errS0 <- cbind(
  res1[, 1:6],
  res1[, -(1:6)][, (2*nlambda+2*runs+1):(2*nlambda+3*runs)]
)
t_l1err_Z <- cbind(
  res1[, 1:6],
  res1[, -(1:6)][, (2*nlambda+3*runs+1):(2*nlambda+4*runs)]
)

l1err_mu <- cbind( t_l1err_mu[,1:6], mean=rowMeans(t_l1err_mu[,-c(1:6)]),
                sd=apply(t_l1err_mu[,-c(1:6)], 1, sd) )
l1err_mu$mean <- round(l1err_mu$mean, 3)
l1err_mu$sd <- round(l1err_mu$sd, 3)
l1err_mu <- l1err_mu %>% arrange(n, tnsr, family, structure, estimator, est_fun)

l1err <- cbind( t_l1err[,1:6], mean=rowMeans(t_l1err[,-c(1:6)]),
                sd=apply(t_l1err[,-c(1:6)], 1, sd) )
l1err$mean <- round(l1err$mean, 3)
l1err$sd <- round(l1err$sd, 3)
l1err <- l1err %>% arrange(n, tnsr, family, structure, estimator, est_fun)

l1errS0 <- cbind( t_l1errS0[,1:6], mean=rowMeans(t_l1errS0[,-c(1:6)]),
                  sd=apply(t_l1errS0[,-c(1:6)], 1, sd) )
l1errS0$mean <- round(l1errS0$mean, 3)
l1errS0$sd <- round(l1errS0$sd, 3)
l1errS0 <- l1errS0[order(l1errS0$family),]
l1errS0 <- l1errS0 %>% arrange(n, tnsr, family, structure, estimator, est_fun)

l1err_Z <- cbind( t_l1err_Z[,1:6], mean=rowMeans(t_l1err_Z[,-c(1:6)]),
                  sd=apply(t_l1err_Z[,-c(1:6)], 1, sd) )
l1err_Z$mean <- round(l1err_Z$mean, 3)
l1err_Z$sd <- round(l1err_Z$sd, 3)
l1err_Z <- l1err_Z[order(l1err_Z$family),]
l1err_Z <- l1err_Z %>% arrange(n, tnsr, family, structure, estimator, est_fun)
