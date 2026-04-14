This repository contains the source code for the paper "Vine Copula Knockoffs for Variable Selection in Gene Expression Studies"

The computational implementations of the proposed methods were carried out using R, employing versions 3.6.3.

The code for running the simulations and the application to a real dataset uses several packages from R. The packages (with versions) are the following:

ggplot2_3.3.5            survival_3.4-0         dplyr_1.0.9            doParallel_1.0.17     
iterators_1.0.14         knockoff_0.3.5         foreach_1.5.2          glmnet_4.1-4          
Matrix_1.5-3             TSP_1.2-4              VineCopula_2.4.5       rvinecopulib_0.6.2.1.3
seqknockoff_0.0.0.9000   latentcor_1.2.0 
knockofftools 1.0.0

The Simulations folder contains code implementing the four distinct data-generating processes (DGPs) for the predictor matrix 
X considered in the paper: a multivariate normal distribution (Section 3.1), a heavy-tailed Markov chain (Section 3.2), survival regression with a block-correlated AR(1) structure (Section 3.3), and a star-type dependence structure induced by a C-vine copula (Appendix B). It also includes code for the remaining appendices, which examine computational time and the impact of key hyperparameters on the proposed procedures.

The Application folder provides code for applying the proposed methodology to a gene expression survival dataset (TCGA_eset) in the context of ovarian cancer research. The TCGA_eset dataset is available through the CuratedOvarianData R package.
