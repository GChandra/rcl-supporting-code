library(tidyverse)
library(mice)
library(MASS)

select <- dplyr::select

setwd("~/supporting_code")

dat <- read.csv("package-idata-80.2024-11-01/Urine and Serum Metabolomics Results/IDATA_WAS_NCIA030321_urine_20240827/IDATA_WAS_NCIA030321_urine_rslts.csv")

urn <- dat %>% filter(RSLT_TYPE=="Scaled", QC_FLAG==0)
# remove duplicates
urn <- urn[!( urn %>% select(SUBJECT_ID, VISITTYPE) %>% duplicated() ), ]
# select only relevant variables
urn <- urn %>% select(-c(CLIENT_IDENTIFIER, CLIENT_SAMPLE_ID, STUDY_VAR,
                         QC_FLAG, DUP_FLAG, RSLT_TYPE))

# remove metabolites with any subjects missing all their replicates
metabs <- colnames(urn)[-c(1:2)]
miss_sum <- rep(0, length(metabs))
for (metab in metabs){
  for (subj in unique(urn$SUBJECT_ID)){
    miss_sum[which(metab==metabs)]  <- miss_sum[which(metab==metabs)] +
      all(is.na( urn[urn$SUBJECT_ID==subj, metab] ))
  }
}

metabs_sel <- metabs[miss_sum==0]
urn_sel <- select(urn, SUBJECT_ID, VISITTYPE, all_of(metabs_sel))

# convert visit to categorical
urn_dat <- urn_sel %>% mutate(VISIT=ifelse(VISITTYPE=="Visit 1", 1, 2)) %>%
  select(-VISITTYPE)
urn_dat$VISIT <- as.factor(urn_dat$VISIT)

urn_transf <- urn_dat

# clean up
rm(urn)
gc()

# test error assumptions
summ_dat <- urn_transf %>%
  select(-VISIT) %>%
  group_by(SUBJECT_ID) %>%
  summarise(across( .cols=metabs_sel,
                    .fns=list(mean=mean, sd=sd),
                    .names="{.col}.{.fn}")
  )

summ_dat_grouped <- summ_dat %>% pivot_longer(
  cols=-SUBJECT_ID,
  names_to=c("metabolite", "summary"),
  names_sep="\\.",
  values_to="level"
) %>%
  pivot_wider(names_from=summary, values_from=level) %>%
  arrange(SUBJECT_ID, metabolite)

err_test <- list()
err_slope <- list()
for (metab in metabs_sel){
  fixef <- lm(sd ~ mean, data=summ_dat_grouped, subset=(metabolite==metab))
  
  # sometimes there are no replicates to estimate mean; in this case, set p-value
  # to 0 so it will be removed in the next step.
  if ( nrow(summary(fixef)$coefficients) == 2 ){
    err_test[[metab]] <- summary(fixef)$coefficients[2,4]
    err_slope[[metab]] <- summary(fixef)$coefficients[2,1]
  } else {
    err_test[[metab]] <- 0
    err_slope[[metab]]
  }
}

metabs_rem <- (err_test < 0.01) & (err_slope > 0.5)
metabs_sel <- names( metabs_rem[metabs_rem==FALSE] )

summ_subdat_grouped <- summ_dat_grouped %>%
  filter( metabolite %in% sample(metabs_sel, size=25, replace=FALSE) )

transformed_metab_dat <- urn_transf %>%
  select(SUBJECT_ID, VISIT, all_of(metabs_sel))
# Save results for the next step
saveRDS(transformed_metab_dat, "preprocessed_dat.rds")
