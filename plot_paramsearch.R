library(tidyverse)
library(ggpubr)

select <- dplyr::select

# Set this to the appropriate working directory on your computer, where all
# downloaded supporting code is stored.
wd <- "~/supporting_code"
setwd(wd)

################################## tnsr figure #################################
res <- read.csv("param_search_tnsr.csv", header=TRUE)
# Rename options for cleanliness of plots
res$structure <- ifelse(res$structure=="diag", "Diagonal",
                        ifelse(res$structure=="block", "Block", "NA"))
res$family <- ifelse(res$family=="gaussian", "Gaussian",
                     ifelse(res$family=="binomial", "Binomial", "NA"))

# filter just a few tnsr values to avoid crowding the figure.
res <- res %>% filter(tnsr %in% c(0.1, 0.3, 0.5, 0.8))
leg <- "attenuation factor" # legend title

# Reproduce parameters from simulation studies
nlambda <- 100
s_0 <- 10
p <- 1000

n_models <- nrow(distinct(res, estimator))
n_opts <- nrow(distinct(res, tnsr))
rng <- 1:(n_models*n_opts)
# Reorganize results to long format
roc1 <- data.frame(ind=rep(1:nlambda, 2), Metric=rep(c("TP","FP"), each=nlambda),
                   t(as.matrix( res[rng, -(1:10)][, 1:(2*nlambda)] )))
# Rename columns
colnames(roc1)[-c(1:2)] <- apply(
  distinct(res, tnsr, estimator),
  MARGIN=1, paste, collapse="_")
# Reshape again to long format.
roc2 <- roc1 %>% gather(key="estimator", value="Rate", -c(ind, Metric)) %>%
  reshape(idvar=c("ind", "estimator"), timevar="Metric",
          direction="wide")
# Separate varying options into separate fields
roc <- data.frame(separate_wider_delim(roc2, cols=estimator, delim="_",
                                       names=colnames(res)[c(9,10)]))
# Convert TP/FP incidence into TP/FP rates
roc$Rate.TP <- roc$Rate.TP / (s_0)
roc$Rate.FP <- roc$Rate.FP / (p-s_0)

# Store subplot for Figure 3
p_tnsr <- ggplot(data=roc, aes(x=Rate.FP, y=Rate.TP, color=tnsr, linetype=estimator)) +
  geom_line(linewidth=0.8) +
  guides(color = guide_legend(order = 2), 
         shape = guide_legend(order = 1)) +
  theme(panel.background=element_rect(fill="white", color="white"),
        # remove the vertical grid lines
        panel.grid.major.x = element_line( linewidth=.1, color="white"),
        panel.grid.minor.x = element_line( linewidth=.1, color="white"),
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( linewidth=.1, color="black"),
        panel.grid.minor.y = element_line( linewidth=.1, color="black")
  ) +
  labs(x="False Positive Rate", y="True Positive Rate",
       shape="estimator", color=leg)

################################### rho figure #################################
# Perform the same set of operations to generate subplot for rho
res <- read.csv("param_search_rho.csv", header=TRUE)
res$structure <- ifelse(res$structure=="diag", "Diagonal",
                        ifelse(res$structure=="block", "Block", "NA"))
res$family <- ifelse(res$family=="gaussian", "Gaussian",
                     ifelse(res$family=="binomial", "Binomial", "NA"))

res <- res %>% filter(rho %in% c(0.1, 0.5, 0.7, 0.8))
leg <- "rho"

nlambda <- 100
s_0 <- 10
p <- 1000

n_models <- nrow(distinct(res, estimator))
n_opts <- nrow(distinct(res, rho))
rng <- 1:(n_models*n_opts)
roc1 <- data.frame(ind=rep(1:nlambda, 2), Metric=rep(c("TP","FP"), each=nlambda),
                   t(as.matrix( res[rng, -(1:10)][, 1:(2*nlambda)] )))
colnames(roc1)[-c(1:2)] <- apply(
  distinct(res, rho, estimator),
  MARGIN=1, paste, collapse="_")

roc2 <- roc1 %>% gather(key="estimator", value="Rate", -c(ind, Metric)) %>%
  reshape(idvar=c("ind", "estimator"), timevar="Metric",
          direction="wide")
roc <- data.frame(separate_wider_delim(roc2, cols=estimator, delim="_",
                                       names=colnames(res)[c(9,10)]))

roc$Rate.TP <- roc$Rate.TP / (s_0)
roc$Rate.FP <- roc$Rate.FP / (p-s_0)

p_rho <- ggplot(data=roc, aes(x=Rate.FP, y=Rate.TP, color=rho, linetype=estimator)) +
  geom_line(linewidth=0.8) +
  guides(color = guide_legend(order = 2), 
         shape = guide_legend(order = 1)) +
  theme(panel.background=element_rect(fill="white", color="white"),
        # remove the vertical grid lines
        panel.grid.major.x = element_line( linewidth=.1, color="white"),
        panel.grid.minor.x = element_line( linewidth=.1, color="white"),
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( linewidth=.1, color="black"),
        panel.grid.minor.y = element_line( linewidth=.1, color="black")
  ) +
  labs(x="False Positive Rate", y="True Positive Rate",
       shape="estimator", color=leg)

################################## beta figure #################################
# Perform the same set of operations to generate subplot for rho
res <- read.csv("param_search_beta.csv", header=TRUE)
res$structure <- ifelse(res$structure=="diag", "Diagonal",
                        ifelse(res$structure=="block", "Block", "NA"))
res$family <- ifelse(res$family=="gaussian", "Gaussian",
                     ifelse(res$family=="binomial", "Binomial", "NA"))

res <- res %>% filter(beta %in% c(0.5, 1, 2, 10))
leg <- "beta"

nlambda <- 100
s_0 <- 10
p <- 1000

n_models <- nrow(distinct(res, estimator))
n_opts <- nrow(distinct(res, beta))
rng <- 1:(n_models*n_opts)
roc1 <- data.frame(ind=rep(1:nlambda, 2), Metric=rep(c("TP","FP"), each=nlambda),
                   t(as.matrix( res[rng, -(1:10)][, 1:(2*nlambda)] )))
colnames(roc1)[-c(1:2)] <- apply(
  distinct(res, beta, estimator),
  MARGIN=1, paste, collapse="_")

roc2 <- roc1 %>% gather(key="estimator", value="Rate", -c(ind, Metric)) %>%
  reshape(idvar=c("ind", "estimator"), timevar="Metric",
          direction="wide")
roc <- data.frame(separate_wider_delim(roc2, cols=estimator, delim="_",
                                       names=colnames(res)[c(9,10)]))

roc$Rate.TP <- roc$Rate.TP / (s_0)
roc$Rate.FP <- roc$Rate.FP / (p-s_0)

p_beta <- ggplot(data=roc, aes(x=Rate.FP, y=Rate.TP, color=beta, linetype=estimator)) +
  geom_line(linewidth=0.8) +
  guides(color = guide_legend(order = 2), 
         shape = guide_legend(order = 1)) +
  theme(panel.background=element_rect(fill="white", color="white"),
        # remove the vertical grid lines
        panel.grid.major.x = element_line( linewidth=.1, color="white"),
        panel.grid.minor.x = element_line( linewidth=.1, color="white"),
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( linewidth=.1, color="black"),
        panel.grid.minor.y = element_line( linewidth=.1, color="black")
  ) +
  labs(x="False Positive Rate", y="True Positive Rate",
       shape="estimator", color=leg)

################################## final figure ################################
# Put subplots together
p_roc <- list(p_beta, p_tnsr, p_rho)
figure <- ggarrange(plotlist=p_roc,
                    labels=c("a", "b", "c"),
                    common.legend=FALSE,
                    ncol = 2, nrow = 2)
figure