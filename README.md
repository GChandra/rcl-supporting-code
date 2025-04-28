This zip file contains all R code required to reproduce the figures and tables of this article, as well as those in its supporting information. While code to analyze the data used in the data analysis section is provided, the data itself is not provided. Below are descriptions of each of the files in this folder.

## run_sim.R

This file contains a function (run_sim) used to run many of the simulation studies discussed in the article. It also contains a function (get_elbow) used to obtain the elbow plots for the implementation of CS Lasso for logistic regression, discussed in the supporting information. Finally, it contains a helper function (choose_matrix) used to generate diagonal and block-diagonal covariance matrices, which are used frequently in the simulations of this article.

## figure1.R

This script can be used to generate Figure 1.

## main.R

This script should be run to perform the simulation study in Section 5.1 of the article. Ensure that the working directory is set to the location of this folder on your computer. This script will produce a csv file, results.csv, which will be used to generate Figure 2 and Table 1.

## options.csv

This file provides the parameters used for the simulation study in Section 5.1 of the article. It is used in main.R.

## plot_main.R

This script should be run to generate Figure 2. It also produces data frames which can be used to recreate Table 1. Ensure that the working directory is set to the location of this folder on your computer.

## main_paramsearch.R

This script should be run to perform the simulation studies in Section 5.3 of the article. Ensure that the working directory is set to the location of this folder on your computer. This script should be used to produce 3 csv files: param_search_tnsr.csv, param_search_rho.csv, and param_search_beta.csv. To produce each of these, the appropriate lines of the script must be changed to 'tnsr', 'rho', or 'beta', as specified in the comments in the script. These csv files will be used to generate Figure 3.

## run_sim_param_search.R

Analogous to run_sim.R, except that it is used specifically for the simulations of Section 5.3 of the article.

## plot_paramsearch.R

Analogous to plot_main.R, except that it is used to generate Figure 3 of the article. Ensure that the working directory is set to the location of this folder on your computer.

## options_beta.csv

This file provides the parameters used for the simulation study in Section 5.3 of the article, specifically for Figure 3, subplot (a).

## options_tnsr.csv

This file provides the parameters used for the simulation study in Section 5.3 of the article, specifically for Figure 3, subplot (b).

## options_rho.csv

This file provides the parameters used for the simulation study in Section 5.3 of the article, specifically for Figure 3, subplot (c).

## data_analysis_filtering.R

Contains code to preprocess metabolomics and dietary recall data used in the data analysis section of the article. This code produces an RDS file, preprocessed_dat.rds. This will be used to generate the results of Table 2.

## rcl_analysis.R

Runs Naive Lasso, RC Lasso, and CS Lasso on preprocessed metabolomics and dietary recall data to produce the results of Table 2.

## main_mu_Z.R

Analogous to main.R, except that it is used to run the simulations of Section S1 of the Supporting Information.

## run_sim_mu_Z.R

Analogous to run_sim.R, except that it is used to run the simulations of Section S1 of the Supporting Information.

## options_mu_Z.csv

This file provides the parameters used for the simulation study in Section S1 of the Supporting Information.

## plot_mu_Z.csv

Analogous to plot_main.R, except that it is used to generate Figure S1.1 and Table S1.1 of the Supporting Information. Ensure that the working directory is set to the location of this folder on your computer.

## elbow_sim.R

Analogous to main.R, except that it is used to run the simulations to generate the elbow plots in Figure S2.1 of the Supporting Information.

## elbow_options.R

This file provides the parameters used for the simulation study in Section S2 of the Supporting Information.