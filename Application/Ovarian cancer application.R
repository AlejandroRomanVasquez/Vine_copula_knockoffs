#### libraries ####

library(curatedOvarianData)
library(Biobase)
library(rvinecopulib)
library(VineCopula)
library(TSP)
library(parallel)
library(glmnet)
library(foreach)
library(knockoff)
library(doParallel)
library(dplyr)
library(tidyr)
library(survival)
library(ggplot2)
library(ncvreg)
library(knockofftools)
library(clustermq)

##### Auxiliary functions #####

#Custom function to select the most variable genes
#
#method: var stands for variance and cv for coefficient of variation 
#remove_multimapped: Some columns have a label containing two or more genes that must be removed.

select_top_genes <- function(eset,
                             k = 1000,
                             method = c("var", "cv"),
                             remove_multimapped = TRUE,
                             remove_missing = TRUE) {
  
  method <- match.arg(method)
  
  ## Extract expression matrix
  X <- exprs(eset)
  gene_names <- featureNames(eset)
  
  ## Filter probes
  keep <- rep(TRUE, length(gene_names))
  
  if (remove_multimapped) {
    keep <- keep & !grepl("///", gene_names)
  }
  
  if (remove_missing) {
    keep <- keep & !is.na(gene_names) & gene_names != ""
  }
  
  X <- X[keep, , drop = FALSE]
  gene_names <- gene_names[keep]
  rownames(X) <- gene_names
  
  ## Compute selection score
  score <- switch(
    method,
    var = apply(X, 1, var, na.rm = TRUE),
    cv  = {
      mu <- rowMeans(X, na.rm = TRUE)
      sd <- apply(X, 1, sd, na.rm = TRUE)
      sd / pmax(mu, .Machine$double.eps)
    }
  )
  
  #Select top-k genes
  k <- min(k, nrow(X))
  idx <- order(score, decreasing = TRUE)[seq_len(k)]
  
  X[idx, , drop = FALSE]
}


#A function that bins the continuous variables of a matrix X by 
#a specific quantile corresponding to the given probability prob

Continuous_to_binary_by_quantile <- function(X, column_type, prob=0.5){
  
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  
  if (is.null(column_type)) {
    stop("Argument column_type is missing")  
  }
  
  #Number of variables p 
  p <- dim(X)[2]
  
  X_mixed <- X
  for(i in 1:p) {   
    if (column_type[i]=="num"){
      X_mixed[,i] <- X[,i]
    }
    else {
      quantile_i <- as.numeric(quantile(X[,i],prob=prob)) 
      X_mixed[,i] <- ifelse(X[,i]<quantile_i,0,1)
    } 
  }
  return(X_mixed)
}

#Ordinal to uniform transformation (more than 10 levels is considered a numerical variable)
ordinal_to_uniform <- function(column){
  
  if (is.null(column)) {
    stop("Argument column is missing")  
  }
  
  # A discrete variable with more than 10 levels is considered a numerical variable
  if(length(unique(column)) < 11){
    
    column_lower <- column - 1
    
    #Empirical Cumulative Distribution Function
    ecdf_column <- ecdf(column)  
    
    u_upper <-   ecdf_column(column)
    u_lower <-  ecdf_column(column_lower)
    
    u_cond <- runif(n=length(column),min=u_lower, max= u_upper)
    
    return (u_cond)
  }
  
  else{
    return (column)
  }
}


#Function to identify if the variable is numerical(num) or ordinal (ord)
column_type_identification <- function(column){ 
  
  if (is.null(column)) {
    stop("Argument column is missing")  
  }
  
  # A discrete variable with more than 10 levels is considered a numerical variable
  if(length(unique(column)) < 11) 
  {type<-"ord"}
  else
  {type<-"num"} 
  return(type)
}    

#Uniform to original transformation
uniform_to_original <- function(X, u_Xk, column_type){
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  
  if (is.null(u_Xk)) {
    stop("Argument u_Xk is missing")  
  }
  
  if (is.null(column_type)) {
    stop("Argument column_type is missing")  
  }
  
  #Number of variables p 
  p <- dim(X)[2]
  
  Xk <- X
  for(i in 1:p) {   
    if (column_type[i]=="num"){
      Xk[,i] <- as.vector(quantile(X[,i], probs = punif(u_Xk[,i],min=0, max=1), type=8))
    }
    else {
      Xk[,i]<- round(as.vector(quantile(X[,i],  probs = u_Xk[,i], type=1)),0)
    } 
  }
  return(Xk)
  
}


#### Utility functions for the e-values procedure ####

#These functions are taken from
#https://github.com/zhimeir/derandomized_knockoffs_fdr

# The eBH procedure
### Input: 
###   E: e-values
###   alpha: target FDR level
### Output:
###   Variables selected by the e-BH procedure

ebh <- function(E, alpha){
  
  p <- length(E)
  E_ord <- order(E, decreasing = TRUE)
  E <- sort(E, decreasing = TRUE)
  comp <- E >= (p / alpha / (1:p))
  id <- suppressWarnings(max(which(comp>0)))
  if(id > 0){
    rej <- E_ord[1:id]
  }else{
    rej <- NULL
  }
  return(list(rej = rej))
}

## Computing the early stopping time ##
### Input:
###   W: vector of knockoff feature importance statistics 
###   gamma: alpha_kn 
###   offset: value between 0 and 1
### Output: 
###   The modified knockoff stopping time defined in (14)

stop_early <- function(W, gamma, offset){
  
  tau <- alphakn_threshold(W, fdr =  gamma, offset = offset) 
  ord_W <- order(abs(W), decreasing = TRUE)
  sorted_W <- W[ord_W]
  
  if(sum(W>0) >= 1 / gamma){
    pos_ind <- which(sorted_W > 0)
    tau1 <- sorted_W[pos_ind[ceiling(1/gamma)-1]]
  }else{
    tau1 <- 0
  }
  tau <- min(tau,tau1) 
  
  return(tau)
}



## Compute stopping time w/ diff alpha_kn and offset ##
### Input:
###   W: a length p vector of knockoff feature importance statistics
###   fdr: the target FDR level
###   offset: 0 or 1 
### Output: 
###   the knockoff selection threshold

alphakn_threshold <- function(W, fdr, offset) {
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}




#### R functions related to vines #####

#get_dvine_order() -->   This function runs an heuristic procedure to determine the 
#order for the first tree in a D-vine structure using the TSP R package to solve 
#the traveling salesman problem. To solve it, we need to identify the shortest 
#Hamiltonian path by assigning weights based on the pairwise Kendall’s τ

#Arguments:
#X: matrix of predictors

#Value: an integer vector with the new indices

get_dvine_order <- function(X, random_order=TRUE){
  
  if (is.null(X)) {
    stop("Argument X is null")  
  }
  
  #Matrix transformation
  X <- as.matrix(X)
  
  #Matrix of 1 - tau_ij
  M_tau <- 1 - abs(TauMatrix(X))
  
  #Hamiltonian path and solution (functions of package TSP)
  hamilton <- insert_dummy(TSP(M_tau), label="cut")
  
  if (random_order==FALSE){
    sol <- solve_TSP(hamilton,method="identity")
  }
  else{
    sol <- solve_TSP(hamilton,method="repetitive_nn")
  }
  
  #Reordering
  TSP_order <- cut_tour(sol,"cut")
  names(TSP_order) <- NULL
  
  return(TSP_order)
  
}

#get_cvine_order() -->   This function runs an heuristic procedure to determine the 
#order of a C-vine structure using the  Kendall's tau:
#   it orders the variables from highest to lowest score, 
#   producing the C-vine order in which the most influential variable 
#   becomes the root of the first tree, 
#   followed by the next most dependent variable, and so on.

#Arguments:
#X: matrix of predictors

get_cvine_order <- function(X){
  
  if (is.null(X)) {
    stop("Argument X is null")  
  }
  
  #Matrix transformation
  X <- as.matrix(X)
  
  #C-vine order determine by Kendall's tau:
  tau <- cor(X, method = "kendall")
  scores <- rowSums(abs(tau))
  Cvine_order <- order(scores, decreasing = TRUE)
  
  return(Cvine_order)
}


#X_Xk_dvine_distributions()-->Function to fit the dvine distribution for X and X_X matrices
#Arguments:
#X_cont: matrix of continuous predictors (after transformation of ordinal to conditional uniforms)
#vinecop_family : String related to the family of pair copulas used in the dvine fitting. 
#                 Common options are "parametric", "nonparametric", "onepar". More details can be found
#                 in the documentation of R package rvinecopulib
#                  https://cran.r-project.org/web/packages/rvinecopulib/rvinecopulib.pdf
# n_cores: int -> number of cores for parallel processing (default 1)
# threshold: parameter for thresholded vine copulas (default 0)
# psi0: prior probability of a non-independence copula (default 0.95)

#Value: This function returns a list that contains objects of class vinecop_dist for X and X_X
#Note: more information about objects of class vinecop_dist can be found in 
#https://cran.r-project.org/web/packages/rvinecopulib/rvinecopulib.pdf

X_Xk_dvine_distributions <- function(X_cont, vinecop_family="parametric", n_cores=1, 
                                     threshold=0, psi0=0.95){
  
  if (is.null(X_cont)) {
    stop("Argument X_cont is null")  
  }
  
  #Number of variables p and sample size n
  n <- dim(X_cont)[1]
  p <- dim(X_cont)[2]
  
  
  #dstructures for dvines
  X_X_dstructure <- dvine_structure((2*p):1)
  X_dstructure <- dvine_structure(p:1)
  
  #Dataset column binding
  X_X <- cbind(X_cont,X_cont)
  
  #Seudo-Observations
  u_X <- pseudo_obs(X_cont)
  u_X_X <- pseudo_obs(X_X)
  
  #Fitting dvine distribution for X_X
  dvine_fitting_time <- system.time(
    fit_dvine_trunc <- vinecop(u_X_X, family_set=c(vinecop_family), structure= X_X_dstructure, 
                               presel=TRUE, selcrit='mbicv', par_method='mle', 
                               threshold=threshold, psi0=psi0, 
                               show_trace=FALSE, cores=n_cores, trunc_lvl=p-1)
  )
  
  #Printing dvine X_X fitting time
  print("dvine fitting time in seconds:")
  print(dvine_fitting_time)
  
  #Pair-copula list for X_X
  X_X_dvine_pclist <- fit_dvine_trunc$pair_copulas
  
  #dvine distribution for X_Xk 
  X_X_dvine_dist <- vinecop_dist(X_X_dvine_pclist, X_X_dstructure)
  
  #Pair-copula list for X
  X_dvine_pclist <- list(rep(list(""),p-1))
  
  #Iniziating with Independent copula
  for (i in 1:(p-1)){
    bicop <- bicop_dist("indep",)
    X_dvine_pclist[i] <- list(rep(list(bicop),p-i))
  }
  
  #Pair copula list just for X dependencies
  for (i in 1:(p-1)){
    J <- p-i
    
    for (j in 1:J){
      X_dvine_pclist[[i]][j] <- X_X_dvine_pclist[[i]][j] 
      
    } 
  }
  
  # dvine distribution for X
  X_dvine_dist <- vinecop_dist(X_dvine_pclist, X_dstructure)
  
  #List with dvine distributions
  dvine_distributions <- list(X_dvine_dist=X_dvine_dist, 
                              X_X_dvine_dist=X_X_dvine_dist,
                              u_X=u_X)
  
  return(dvine_distributions)
}

#X_Xk_Rvine_distributions()-->Function to fit a C-vine distribution for X and 
# then fit an R-vine for the vector (X,Xk)

#Arguments:
#X_cont: matrix of continuous predictors (after transforming ordinals variables)
#vinecop_family : String related to the family of pair copulas used in the vine fitting. 
#                 Common options are "parametric", "nonparametric", "onepar". More details can be found
#                 in the documentation of R package rvinecopulib
#                  https://cran.r-project.org/web/packages/rvinecopulib/rvinecopulib.pdf
# n_cores: int -> number of cores for parallel processing (default 1)
# threshold: parameter for thresholded vine copulas (default 0)
# psi0: prior probability of a non-independence copula (default 0.95)

#Value: This function returns a list that contains objects of class vinecop_dist for X and X_Xk
#Note: more information about objects of class vinecop_dist can be found in 
#https://cran.r-project.org/web/packages/rvinecopulib/rvinecopulib.pdf

X_Xk_rvine_distributions <- function(X_cont, vinecop_family="parametric", n_cores=1, 
                                     threshold=0, psi0=0.95){
  
  if (is.null(X_cont)) {
    stop("Argument X_cont is null")  
  }
  
  #Number of variables p and sample size n
  n <- dim(X_cont)[1]
  p <- dim(X_cont)[2]
  
  #A C-vine structure with the root of the first tree in the first column, 
  #the next root for the next tree involving the first and second columns, and so on.
  X_Cvine_structure <- cvine_structure(seq(p,1))
  X_Cvine_matrix <- as_rvine_matrix(X_Cvine_structure)
  
  #Seudo-observations for X_cont
  u_X_cont <- pseudo_obs(X_cont)
  
  #Fitting C-vine distribution for u_X_cont
  Cvine_fitting_time <- system.time(
    X_Cvine_fit <- vinecop(data = u_X_cont, family_set=c(vinecop_family), structure= X_Cvine_structure, 
                           presel=TRUE, selcrit='mbicv', par_method='mle', 
                           threshold=threshold, psi0=psi0, show_trace=FALSE, cores=n_cores)
  )
  
  #Printing C-vine fitting time
  print("C-vine fitting time in seconds:")
  print(Cvine_fitting_time)
  
  #Creation of the Rmatrix for the R-vine structure
  Rmatrix <- matrix(rep(0,2*p*2*p), 2*p, 2*p)
  
  #The upper right triangular matrix is a copy of the C-vine matrix of X
  Rmatrix[1:p, 1:(p-1)] <- X_Cvine_matrix[1:p,2:p] + p
  
  #For loops to fill some cells of the matrix
  for (i in 1:p){
    for (j in 0:(p-1)) {
      Rmatrix[i+j, p-j] <- i 
    }
  }
  
  for (i in 1:p){
    Rmatrix[p + i, p + 1 - i] <- p + i
  }
  
  #The upper left triangular matrix is a copy of a part of the C-vine matrix of X
  Rmatrix[1:p, (p+1):(2*p)] <- X_Cvine_matrix[1:p,1:p]
  
  #Rvine_matrix
  Rvine_matrix <- rvine_matrix(Rmatrix)
  
  #R-vine structure for the (X,Xk) vector
  X_Xk_Rvine_structure <- as_rvine_structure(Rvine_matrix)
  
  #Pair-copula list of X from the C-vine fit
  X_Cvine_pclist <- X_Cvine_fit$pair_copulas
  
  #Initiating the construction of pair-copula list for (X,Xk)
  X_Xk_Rvine_pclist <- list(rep(list(""),2*p-1))
  
  #Independent copula
  bicop <- bicop_dist("indep",)
  
  #Filling with independent copula
  for (i in 1:(2*p-1)){
    X_Xk_Rvine_pclist[i] <- list(rep(list(bicop),2*p-i))
  }
  
  #List of pair copulas for the C-vine structure of X
  for (i in 1:(p-1)){
    J <- p-i
    
    for (j in 1:J){
      X_Xk_Rvine_pclist[[i]][j] <- X_Cvine_pclist[[i]][j] 
      
    } 
  }
  
  #List of pair copulas for the C-vine structure of Xk (equal to X)
  for (i in 1:(p-1)){
    J <- 2*p - i
    
    for (j in (p+1):J){
      X_Xk_Rvine_pclist[[i]][j] <- X_Cvine_pclist[[i]][j-p] 
      
    } 
  }
  
  #C-vine distribution of X
  X_Cvine_distribution <- vinecop_dist(X_Cvine_pclist, X_Cvine_structure)
  
  #R-vine distribution of (X,Xk)
  X_Xk_Rvine_distribution <- vinecop_dist(X_Xk_Rvine_pclist, X_Xk_Rvine_structure)
  
  
  #List with some vine distributions
  Rvine_distributions <- list(X_Cvine_distribution=X_Cvine_distribution, 
                              X_Xk_Rvine_distribution=X_Xk_Rvine_distribution,
                              u_X_cont=u_X_cont)
  
  return(Rvine_distributions)
}



# create_dvine_Knockoffs()--> Function to sample dvine knockoffs

#Arguments:
#X: matrix of original predictors
#X_cont: matrix of continuous predictors (after transformation of ordinal to conditional uniforms)
#column_type: vector with the type of each variable of X ("num" numeric, "ord" ordinal)
#dvine_distributions: a list containing: 
#   1) object of class vinecop_dist for X_X,
#   2) Object of class vinecop_dist for X_cont,
#   3) Pseudo observations of X_cont saved in u_X
#n_cores: int -> number of cores for parallel processing
#seed_Xk: integer value to reproduce the sampled Knockoffs. The default is NULL to ensure a random knockoffs simulation 
#Note: more information about objects of class vinecop_dist can be found in 
#https://cran.r-project.org/web/packages/rvinecopulib/rvinecopulib.pdf

#Value: This function returns a matrix Xk of knockoffs 

create_dvine_Knockoffs <- function(X, X_cont, column_type, dvine_distributions, n_cores=1, seed_Xk=NULL){
  
  if (is.null(X)) {
    stop("Argument X is null")  
  }
  if (is.null(X_cont)) {
    stop("Argument X_cont is null")  
  }
  if (is.null(column_type)) {
    stop("Argument colum_type is null")  
  }
  if (is.null( dvine_distributions)) {
    stop("Argument dvine_distributions is null")  
  }
  
  set.seed(seed_Xk) #For reproducibility
  
  #Number of variables p and sample size n
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Pseudo observations
  u_X <- dvine_distributions$u_X
  
  #Independent uniforms w
  w_X <- rosenblatt(x=u_X, model=dvine_distributions$X_dvine_dist, cores = n_cores)
  w_Xk <- matrix(runif(n=p*n,min=0,max=1),nrow=n,ncol=p)
  w_X_Xk <- cbind(w_X,w_Xk)
  
  #Knockoff sampling Xk
  u_X_Xk <- inverse_rosenblatt(u=w_X_Xk, model= dvine_distributions$X_X_dvine_dist, cores = n_cores)
  u_Xk <- u_X_Xk[,(p+1):(2*p)]
  
  #Transformation to original predictors
  Xk <- uniform_to_original(X, u_Xk, column_type ) 
  
  return(Xk)
}

# create_cvine_Knockoffs()--> Function to sample C-vine knockoffs

#Arguments:
#X: matrix of original predictors
#X_cont: matrix of continuous predictors (after transforming ordinals predictors)
#column_type: vector with the type of each variable of X ("num" numeric, "ord" ordinal)
#rvine_distributions: a list containing: 
#   1) Object of class vinecop_dist (R-vine) for X_Xk,
#   2) Object of class vinecop_dist (C-vine) for X_cont,
#   3) Pseudo observations of X_cont saved in u_reordered_X_cont
#n_cores: int -> number of cores for parallel processing
#seed_Xk: integer value to reproduce the sampled Knockoffs. The default is NULL to ensure a random knockoffs simulation 
#Note: more information about objects of class vinecop_dist can be found in 
#https://cran.r-project.org/web/packages/rvinecopulib/rvinecopulib.pdf

#Value: This function returns a matrix Xk of knockoffs 

create_cvine_Knockoffs <- function(X, X_cont, column_type,
                                   rvine_distributions, n_cores=1, 
                                   seed_Xk=NULL){
  
  if (is.null(X)) {
    stop("Argument X is null")  
  }
  if (is.null(X_cont)) {
    stop("Argument X_cont is null")  
  }
  if (is.null(column_type)) {
    stop("Argument colum_type is null")  
  }
  if (is.null(rvine_distributions)) {
    stop("Argument rvine_distributions is null")  
  }
  
  set.seed(seed_Xk) #For reproducibility
  
  #Number of variables p and sample size n
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Pseudo observations
  u_X_cont <- rvine_distributions$u_X_cont
  
  #Independent uniforms w
  w_X <- rosenblatt(x=u_X_cont, model=rvine_distributions$X_Cvine_distribution, cores = n_cores)
  w_Xk <- matrix(runif(n=p*n,min=0,max=1),nrow=n,ncol=p)
  w_X_Xk <- cbind(w_X,w_Xk)
  
  #Knockoff sampling Xk
  u_X_Xk <- inverse_rosenblatt(u=w_X_Xk, model= rvine_distributions$X_Xk_Rvine_distribution, cores = n_cores)
  u_Xk <- u_X_Xk[,(p+1):(2*p)]
  
  #Transformation to original predictors
  Xk <- uniform_to_original(X, u_Xk, column_type ) 
  
  #Xk column names
  Xk_column_names <- sapply(1:p, function(number) paste0("Xk", number))
  colnames(Xk) <- Xk_column_names
  
  return(Xk)
}


# stable_lasso_glmnet_parallel()--> Function to fit a regularized lasso regresion model using 
#some functions of the R package glmnet.
#It implements the stabilizing procedure of Roberts and Nowak (2014) to diminish sensitivity
#to the fold assignment used in cross-validation to select the hyperparameter lambda

#Arguments:
#X: matrix of predictors
#y: vector or matrix of response
#lasso_family: a string to select linear regression "gaussian" or logistic regression "binomial"
#M_lasso: integer related to the number of runs for the stabilization against CV (Roberts and Nowak, 2014)
#n_folds: integer indicating the number of cross validations
#random_cv: logic value for a random assignation of the folds in the Cross-validation
#        TRUE is the default. FALSE is for reproducibility purposes. 
#Note: more information about the R package glmnet can be found in 
#https://cran.r-project.org/web/packages/glmnet/glmnet.pdf
#Note 2: this function runs in parallel for the stabilization against CV (Roberts and Nowak, 2014)

#Value: This function returns a vector of the estimated coeficientes (without the intercept)


stable_lasso_glmnet_parallel <- function(X, y, lasso_family, M_lasso = 10, n_folds = 5, random_cv=TRUE){
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  if (is.null(y)) {
    stop("Argument y is missing")  
  }
  if (is.null( lasso_family )) {
    stop("Argument lasso_family is missing")  
  }
  
  y_vec <- as.vector(y)
  X_matrix <- as.matrix(X)
  
  #Stabilizing the lasso against CV (Roberts and Nowak, 2014)
  lambdas <- rep(0,M_lasso)
  
  #Random assignation of the folds in the Cross-validation
  if(random_cv==TRUE){
    time_cv <- system.time(  
      lambdas <- foreach(i = 1:M_lasso, .combine=c,.packages=c("glmnet")) %dopar% {
        set.seed(NULL)
        cvfit <- cv.glmnet(X_matrix, y_vec, alpha=1, family = lasso_family, nfolds = n_folds, standardize = TRUE)
        cvfit$lambda.min
      }
    )
  }
  #Non-random assignation of the folds in the Cross-validation
  else{
    time_cv <- system.time(  
      lambdas <- foreach(i = 1:M_lasso, .combine=c,.packages=c("glmnet")) %dopar% {
        set.seed(i)
        cvfit <- cv.glmnet(X_matrix, y_vec, alpha=1, family = lasso_family, nfolds = n_folds, standardize = TRUE)
        cvfit$lambda.min
      }
    )
  }
  
  #Selecting the median of the lambdas distribution
  lambda50 <- as.numeric(quantile(lambdas,probs=0.5))
  fit_coef <- coef(glmnet(X_matrix, y_vec, alpha = 1, lambda = lambda50, family = lasso_family, standardize = TRUE))
  
  fit_coef_vec <- as.vector(fit_coef)
  fit_coef_vec <- fit_coef_vec[-1] 
  
  return(fit_coef_vec)
}

# stable_lasso_glmnet()--> Function to fit a regularized lasso regresion model using 
#some functions of the R package glmnet.
#It implements the stabilizing procedure of Roberts and Nowak (2014) to diminish sensitivity
#to the fold assignment used in cross-validation to select the hyperparameter lambda

#Arguments:
#X: matrix of predictors
#y: vector of response (In the context of survival regression, y needs to be a survival object from the Survival package)
#lasso_family: a string to select linear regression "gaussian", logistic regression "binomial" 
#               or "cox" for survival cox regression model
#M_lasso: integer related to the number of runs for the stabilization against CV (Roberts and Nowak, 2014)
#n_folds: integer indicating the number of cross validations
#random_cv: logic value for a random assignation of the folds in the Cross-validation
#        TRUE is the default. FALSE is for reproducibility purposes. 
#Note: more information about the R package glmnet can be found in 
#https://cran.r-project.org/web/packages/glmnet/glmnet.pdf
#Note 2: this function runs in parallel for the stabilization against CV (Roberts and Nowak, 2014)

#Value: This function returns a vector of the estimated coeficientes (without the intercept)


stable_lasso_glmnet <- function(X, y, lasso_family, M_lasso = 10, n_folds = 5, random_cv=TRUE){
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  if (is.null(y)) {
    stop("Argument y is missing")  
  }
  if (is.null( lasso_family )) {
    stop("Argument lasso_family is missing")  
  }
  
  if(lasso_family=="gaussian"){
    y_glmnet <- as.vector(y)}
  else if(lasso_family=="binomial"){
    y_glmnet  <- as.factor(y)}
  else { #This is the case for the survival regression model
    # y should be a survival object
    y_glmnet <- y}
  
  X_matrix <- as.matrix(X)
  
  #Stabilizing the lasso against CV (Roberts and Nowak, 2014)
  lambdas <- rep(0,M_lasso)
  
  #Random assignation of the folds in the Cross-validation
  if(random_cv==TRUE){
    start_time <- Sys.time()  
    for (i in 1:M_lasso){
      set.seed(NULL)
      cvfit <- cv.glmnet(X_matrix, y_glmnet, alpha=1, family = lasso_family, nfolds = n_folds, standardize = TRUE)
      lambdas[i] <- cvfit$lambda.min
    }
    end_time <- Sys.time()
    print("Time for running Stabilizing the lasso:")
    print(end_time-start_time)
  }
  #Non-random assignation of the folds in the Cross-validation
  else{
    start_time <- Sys.time() 
    for( i in 1:M_lasso){
      set.seed(i)
      cvfit <- cv.glmnet(X_matrix, y_glmnet, alpha=1, family = lasso_family, nfolds = n_folds, standardize = TRUE)
      lambdas[i] <- cvfit$lambda.min
    }
    end_time <- Sys.time()
    print("Time for running Stabilizing the lasso:")
    print(end_time-start_time)
  }
  
  #Selecting the median of the lambdas distribution
  lambda50 <- as.numeric(quantile(lambdas,probs=0.5))
  fit_coef <- coef(glmnet(X_matrix, y_glmnet, alpha = 1, lambda = lambda50, family = lasso_family, standardize = TRUE))
  
  fit_coef_vec <- as.vector(fit_coef)
  
  if(lasso_family=="cox")
    fit_coef_vec <- fit_coef_vec
  else {
    fit_coef_vec <- fit_coef_vec[-1]   
  }
  return(fit_coef_vec)
}


# ekn_dvines()--> Function to derandomized knockoffs using e-values for FDR control. This function
# considers the dvine knockoff procedure.
#The code to implement this function is adapted from 
#https://github.com/zhimeir/derandomized_knockoffs_fdr

#Arguments:
#X: matrix of predictors
#X_cont: matrix of continuous predictors (after transformation of ordinal to conditional uniforms)
#column_type: vector with the type of each variable of X ("num" numeric, "ord" ordinal)
#y: vector or matrix of response
#dvine_distributions: a list containing: 
#   1) object of class vinecop_dist for X_X,
#   2) Object of class vinecop_dist for X_cont,
#   3) Pseudo observations of X_cont saved in u_X
#M: integer denoting the number of generated copies of the knockff matrix Xk.
#M_lasso: integer related to the number of runs for the stabilzation against CV
#alpha: integer indicating FDR target level
#gamma: integer denoting target level for the knockoff threshold. According to Ren & Barber (2023),
#       experimentally, gamma=alpha/2 works well.           
#lasso_family: a string to select linear regression "gaussian" or logistic regression "binomial" 
#n_cores: int -> number of cores for parallel processing
#random_cv: logic value for a random assignation of the folds in the Cross-validation used in LASSO regression
#        TRUE is the default. FALSE is for reproducibility purposes. 
#random_evalues: logic value for a random sampling of the knockoffs in the derandomized knockoff procedure
#        TRUE is the default. FALSE is for reproducibility purposes. 


#Note: the knockoff.threshold() function from the R knockoff package is used for 
#setting the Knockoff rejection threshold (https://cran.r-project.org/web/packages/knockoff/knockoff.pdf)

#Value: This function returns a list with the selected variables of the procedure

ekn_dvines <- function(X, X_cont, column_type, y, dvine_distributions, M=50, M_lasso=10, alpha=0.2, 
                       gamma=0.1, lasso_family, n_cores=1, n_folds = 5, random_cv=TRUE, random_evalues=TRUE){
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  if (is.null( X_cont )) {
    stop("Argument X_cont is missing")  
  }
  if (is.null( column_type )) {
    stop("Argument column_type is missing")  
  }
  if (is.null(y)) {
    stop("Argument y is missing")  
  }
  
  if (is.null( dvine_distributions )) {
    stop("Argument dvine_distributions is missing")  
  }
  if (is.null( lasso_family )) {
    stop("Argument lasso_family is missing")  
  }
  
  #Number of variables p and sample size n  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Initial matrix of E-values 
  E <- matrix(0, M, p)
  
  # Creating empty lists to store the X_Xk matrices
  ls_X_Xk <- list()
  
  for(m in 1:M){
    
    if(random_evalues==TRUE){
      #dvine Knockoffs sampling
      Xk <- create_dvine_Knockoffs(X=X, X_cont=X_cont, column_type=column_type, 
                                   dvine_distributions=dvine_distributions, 
                                   n_cores=n_cores,seed_Xk=NULL)
    }
    else {
      #dvine Knockoffs sampling
      Xk <- create_dvine_Knockoffs(X=X, X_cont=X_cont, column_type=column_type, 
                                   dvine_distributions=dvine_distributions, 
                                   n_cores=n_cores, seed_Xk=m)    
    }
    #X and Xk column binding
    X_Xk <- cbind(X, Xk)
    
    # Storing the data in lists
    ls_X_Xk <- append(ls_X_Xk, list(X_Xk))
    
  }
  
  start_time <- Sys.time()  
  Z <- foreach(m = 1:M,.packages=c("glmnet"),.export = "stable_lasso_glmnet") %dopar% {
    
    stable_lasso_glmnet(X=ls_X_Xk[[m]], y=y, lasso_family=lasso_family,  
                        M_lasso=M_lasso, n_folds=n_folds, 
                        random_cv=random_cv)
  }
  end_time <- Sys.time()
  print("Time for getting Z:")
  print(end_time-start_time)
  
  for(m in 1:M){
    #Importance statistic
    W <- abs(Z[[m]][1:p]) - abs(Z[[m]][(p+1):length(Z[[m]])])
    
    #Knockoff rejection threshold - conservative procedure ("knockoffs+" offset = 1)
    tau <- stop_early(W, gamma, offset=1) 
    
    #E-vales for all the variables (columns) for m run
    E[m,] <- (W >= tau) / (1 + sum(W <= -tau))
    
  }
  
  #Averaging the e-values to select set of discoveries
  E <- p*colMeans(E)
  rej <- ebh(E, alpha)$rej
  
  return(list(rej = rej, E = E)) 
  
}
# ekn_cvines()--> Function to derandomized knockoffs using e-values for FDR control. 
# This function considers the C-vine knockoff procedure.
# Part of the  code to implement this function is adapted from 
# https://github.com/zhimeir/derandomized_knockoffs_fdr

#Arguments:
#X: matrix of predictors
#X_cont: matrix of continuous predictors (after transforming ordinal to continuous)
#column_type: vector with the type of each variable of X ("num" numeric, "ord" ordinal)
#y: vector or matrix of response
#rvine_distributions: a list containing: 
#   1) Object of class vinecop_dist (R-vine) for X_Xk,
#   2) Object of class vinecop_dist (C-vine) for X,
#   3) Pseudo observations of X_cont saved in u_X_cont
#M: integer denoting the number of generated copies of the knockff matrix Xk.
#M_lasso: integer related to the number of runs for the stabilzation against CV
#alpha: integer indicating FDR target level
#gamma: integer denoting target level for the knockoff threshold. According to Ren & Barber (2023),
#       experimentally, gamma=alpha/2 works well.           
#lasso_family: a string to select linear regression "gaussian" or logistic regression "binomial" 
#n_cores: int -> number of cores for parallel processing
#random_cv: logic value for a random assignation of the folds in the Cross-validation used in LASSO regression
#        TRUE is the default. FALSE is for reproducibility purposes. 
#random_evalues: logic value for a random sampling of the knockoffs in the derandomized knockoff procedure
#        TRUE is the default. FALSE is for reproducibility purposes. 


#Note: the knockoff.threshold() function from the R knockoff package is used for 
#setting the Knockoff rejection threshold (https://cran.r-project.org/web/packages/knockoff/knockoff.pdf)

#Value: This function returns a list with the selected variables of the procedure

ekn_cvines <- function(X, X_cont, column_type, y, rvine_distributions, M=50, M_lasso=10, alpha=0.2, 
                       gamma=0.1, lasso_family, n_cores=1, n_folds = 5, random_cv=TRUE, random_evalues=TRUE){
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  if (is.null( X_cont )) {
    stop("Argument X_cont is missing")  
  }
  if (is.null( column_type )) {
    stop("Argument column_type is missing")  
  }
  if (is.null(y)) {
    stop("Argument y is missing")  
  }
  
  if (is.null( rvine_distributions )) {
    stop("Argument rvine_distributions is missing")  
  }
  if (is.null( lasso_family )) {
    stop("Argument lasso_family is missing")  
  }
  
  #Number of variables p and sample size n  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Initial matrix of E-values 
  E <- matrix(0, M, p)
  
  # Creating empty lists to store the X_Xk matrices
  ls_X_Xk <- list()
  
  for(m in 1:M){
    
    if(random_evalues==TRUE){
      #cvine Knockoffs sampling
      Xk <- create_cvine_Knockoffs(X=X, X_cont=X_cont, column_type=column_type, 
                                   rvine_distributions=rvine_distributions, 
                                   n_cores=n_cores,seed_Xk=NULL)
    }
    else {
      #dvine Knockoffs sampling
      Xk <- create_cvine_Knockoffs(X=X, X_cont=X_cont, column_type=column_type, 
                                   rvine_distributions=rvine_distributions, 
                                   n_cores=n_cores, seed_Xk=m)    
    }
    #X and Xk column binding
    X_Xk <- cbind(X, Xk)
    
    # Storing the data in lists
    ls_X_Xk <- append(ls_X_Xk, list(X_Xk))
    
  }
  
  start_time <- Sys.time()  
  Z <- foreach(m = 1:M,.packages=c("glmnet"),.export = "stable_lasso_glmnet") %dopar% {
    
    stable_lasso_glmnet(X=ls_X_Xk[[m]], y=y, lasso_family=lasso_family,  
                        M_lasso=M_lasso, n_folds=n_folds, 
                        random_cv=random_cv)
  }
  end_time <- Sys.time()
  print("Time for getting Z:")
  print(end_time-start_time)
  
  for(m in 1:M){
    #Importance statistic
    W <- abs(Z[[m]][1:p]) - abs(Z[[m]][(p+1):length(Z[[m]])])
    
    #Knockoff rejection threshold - conservative procedure ("knockoffs+" offset = 1)
    tau <- stop_early(W, gamma, offset=1) 
    
    #E-vales for all the variables (columns) for m run
    E[m,] <- (W >= tau) / (1 + sum(W <= -tau))
    
  }
  
  #Averaging the e-values to select set of discoveries
  E <- p*colMeans(E)
  rej <- ebh(E, alpha)$rej
  
  return(list(rej = rej, E = E)) 
  
}


# ekn_second_order()--> Function to derandomized knockoffs using e-values for FDR control. This function
# considers the second-order approximation of the knockoff R package
#The code to implement this function is adapted from 
#https://github.com/zhimeir/derandomized_knockoffs_fdr

#Arguments:
#X: matrix of predictors
#y: vector or matrix of response
#M: integer denoting the number of generated copies of the knockff matrix Xk.
#M_lasso: integer related to the number of runs for the stabilization against Cross Validation  (Roberts and Nowak, 2014)
#alpha: integer indicating FDR target level
#gamma: integer denoting target level for the knockoff threshold. According to Ren & Barber (2023),
#       experimentally, gamma=alpha/2 works well.           
#lasso_family: a string to select linear regression "gaussian" or logistic regression "binomial" 
#n_cores: int -> number of cores for parallel processing
#random_cv: logic value for a random assignation of the folds in the Cross-validation used in LASSO regression
#        TRUE is the default. FALSE is for reproducibility purposes. 
#random_evalues: logic value for a random sampling of the knockoffs in the derandomized knockoff procedure
#        TRUE is the default. FALSE is for reproducibility purposes. 


#Note: the knockoff.threshold() function from the R knockoff package is used for 
#setting the Knockoff rejection threshold (https://cran.r-project.org/web/packages/knockoff/knockoff.pdf)

#Value: This function returns a list with the selected variables of the procedure

ekn_second_order <- function(X, y, M=50, M_lasso=10, alpha=0.2, gamma=0.1, 
                             lasso_family, n_cores=1,n_folds = 5, random_cv=TRUE, random_evalues=TRUE){
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  if (is.null(y)) {
    stop("Argument y is missing")  
  }
  if (is.null( lasso_family )) {
    stop("Argument lasso_family is missing")  
  }
  
  
  #Number of variables p and sample size n  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Matrix of E-values 
  E <- matrix(0, M, p)
  
  # Creating empty lists to store the X_Xk matrices
  ls_X_Xk <- list()
  
  for(m in 1:M){
    
    if(random_evalues==TRUE){
      set.seed(NULL) 
    }
    else {  
      set.seed(m) #The seed is set for reproducibility.
    }
    #Gaussian Knockoffs copy selection from the list object
    Xk <- create.second_order(X)
    
    #X and Xk column binding
    X_Xk <- cbind(X, Xk)
    
    # Storing the data in lists
    ls_X_Xk <- append(ls_X_Xk, list(X_Xk))
    
  }
  
  start_time <- Sys.time()  
  Z <- foreach(m = 1:M,.packages=c("glmnet"),.export = "stable_lasso_glmnet") %dopar% {
    
    stable_lasso_glmnet(X=ls_X_Xk[[m]], y=y, lasso_family=lasso_family,  
                        M_lasso=M_lasso, n_folds=n_folds, 
                        random_cv=random_cv)
  }
  end_time <- Sys.time()
  print("Time for getting Z:")
  print(end_time-start_time)
  
  for(m in 1:M){
    #Importance statistic
    W <- abs(Z[[m]][1:p]) - abs(Z[[m]][(p+1):length(Z[[m]])])
    
    #Knockoff rejection threshold - conservative procedure ("knockoffs+" offset = 1)
    tau <- stop_early(W, gamma, offset=1) 
    
    #E-vales for all the variables (columns) for m run
    E[m,] <- (W >= tau) / (1 + sum(W <= -tau))
    
  }
  
  #Averaging the e-values to select set of discoveries
  E <- p*colMeans(E)
  rej <- ebh(E, alpha)$rej
  
  return(list(rej = rej, E = E)) 
  
}






#### Loading and processing the dataset ####

# ---- Load dataset ---
data(TCGA_eset)

# ---- Extract phenotype (clinical) data ---
pheno <- pData(TCGA_eset)

# ---- Clinical variables to be included ---
clin_var <- c(
  "sample_type",
  "primarysite",
  "summarygrade", 
  "summarystage",
  "tumorstage",
  "grade",
  "age_at_initial_pathologic_diagnosis",
  "days_to_death",
  "vital_status",
  "percent_normal_cells",
  "percent_stromal_cells",
  "percent_tumor_cells"
)

# ---- Subset clinical variables ---
clinical <- pheno %>% dplyr::select(all_of(clin_var))

# ---- Selection of the top-k most variable genes (Coefficient of variation cv) ---
k <- 100
expr_top <- select_top_genes(
  eset = TCGA_eset,
  k = k,
  method = "cv"
)

# ---- Transpose expression data and convert to data frame ---
expr_df <- as.data.frame(t(expr_top))

# ---- Verify sample alignment between clinical and expression data ---
all(rownames(expr_df) == rownames(clinical))

# ---- Combine clinical and gene expression data ---
final_data <- cbind(clinical, expr_df)

# ---- Percentage of missing values in clinical variables ---
round(100 * colSums(is.na(clinical)) / nrow(clinical), 2)

# ---- Imputation of missing values ---
# Categorical variables: impute with the mode
final_data$vital_status <- replace_na(
  final_data$vital_status, 
  names(sort(table(final_data$vital_status), decreasing = TRUE))[1]
)

final_data$primarysite <- replace_na(
  final_data$primarysite, 
  names(sort(table(final_data$primarysite), decreasing = TRUE))[1]
)

final_data$summarygrade <- replace_na(
  final_data$summarygrade, 
  names(sort(table(final_data$summarygrade), decreasing = TRUE))[1]
)

final_data$summarystage <- replace_na(
  final_data$summarystage, 
  names(sort(table(final_data$summarystage), decreasing = TRUE))[1]
)

# Ordinal variables: impute with the most frequent level
final_data$tumorstage <- replace_na(
  final_data$tumorstage, 
  as.integer(names(sort(table(final_data$tumorstage), decreasing = TRUE))[1])
)

final_data$grade <- replace_na(
  final_data$grade, 
  as.integer(names(sort(table(final_data$grade), decreasing = TRUE))[1])
)

# Discrete / continuous variables: impute with the mean
final_data <- final_data %>%
  mutate(days_to_death = dplyr::if_else(
    is.na(days_to_death),
    round(mean(days_to_death, na.rm = TRUE)),
    as.numeric(days_to_death)
  ))

final_data <- final_data %>%
  mutate(age_at_initial_pathologic_diagnosis = dplyr::if_else(
    is.na(age_at_initial_pathologic_diagnosis),
    round(mean(age_at_initial_pathologic_diagnosis, na.rm = TRUE)),
    as.numeric(age_at_initial_pathologic_diagnosis)
  ))

final_data <- final_data %>%
  mutate(percent_normal_cells = dplyr::if_else(
    is.na(percent_normal_cells),
    round(mean(percent_normal_cells, na.rm = TRUE)),
    as.numeric(percent_normal_cells)
  ))

final_data <- final_data %>%
  mutate(percent_stromal_cells = dplyr::if_else(
    is.na(percent_stromal_cells),
    round(mean(percent_stromal_cells, na.rm = TRUE)),
    as.numeric(percent_stromal_cells)
  ))

final_data <- final_data %>%
  mutate(percent_tumor_cells = dplyr::if_else(
    is.na(percent_tumor_cells),
    round(mean(percent_tumor_cells, na.rm = TRUE)),
    as.numeric(percent_tumor_cells)
  ))

# ---- Percentage of missing values after imputation ---
round(100 * colSums(is.na(final_data)) / nrow(final_data), 2)

# ---- Survival outcome encoding ---
# Convert survival time to numeric and event indicator to binary
final_data$days_to_death <- as.numeric(final_data$days_to_death)
final_data$vital_status <- as.numeric(final_data$vital_status == "deceased")

# ---- Binary encoding of selected categorical variables ---
final_data$sample_type   <- ifelse(final_data$sample_type == "healthy", 0, 1)
final_data$primarysite   <- ifelse(final_data$primarysite == "other", 0, 1)
final_data$summarygrade  <- ifelse(final_data$summarygrade == "low", 0, 1)
final_data$summarystage  <- ifelse(final_data$summarystage == "early", 0, 1)

# ---- Censoring rate ---
Censoring_survival <- round(
  100 * as.numeric(table(final_data$vital_status)[1] / length(final_data$vital_status)),
  2
)
print(Censoring_survival)

# ---- Design matrix and survival response ---
X <- final_data %>%
  dplyr::select(-c("days_to_death", "vital_status")) %>%
  as.matrix()

# Overall survival object
y_survival <- survival::Surv(
  final_data$days_to_death,
  final_data$vital_status
)

# ---- Number of predictors ---
p <- dim(X)[2]

#### Univariete Cox Regression model for clinical variables ####

#percent_normal_cells
cox_percent_normal_cells <- coxph(
  Surv(time = days_to_death, event = vital_status) ~ 
    percent_normal_cells,
  data = final_data
)

summary(cox_percent_normal_cells)

#percent_stromal_cells
cox_percent_stromal_cells <- coxph(
  Surv(time = days_to_death, event = vital_status) ~ 
    percent_stromal_cells,
  data = final_data
)

summary(cox_percent_stromal_cells)

#sample_type
cox_sample_type <- coxph(
  Surv(time = days_to_death, event = vital_status) ~ 
    sample_type,
  data = final_data
)

summary(cox_sample_type)

#age_at_initial_pathologic_diagnosis
cox_age_at_initial_pathologic_diagnosis <- coxph(
  Surv(time = days_to_death, event = vital_status) ~ 
    age_at_initial_pathologic_diagnosis,
  data = final_data
)

summary(cox_age_at_initial_pathologic_diagnosis)


#summarygrade
cox_summarygrade <- coxph(
  Surv(time = days_to_death, event = vital_status) ~ 
    summarygrade,
  data = final_data
)

summary(cox_summarygrade)


#summarystage
cox_summarystage <- coxph(
  Surv(time = days_to_death, event = vital_status) ~ 
    summarystage,
  data = final_data
)

summary(cox_summarystage)

#tumorstage
cox_tumorstage <- coxph(
  Surv(time = days_to_death, event = vital_status) ~ 
    tumorstage,
  data = final_data
)

summary(cox_tumorstage)


#grade
cox_grade <- coxph(
  Surv(time = days_to_death, event = vital_status) ~ 
    grade,
  data = final_data
)

summary(cox_grade)


#### Non-parametric DTLDCKe (D-vines) #####

# ---- Parallel computing setup ---
n_cores <- detectCores() - 1           # Number of cores used for parallel execution
cl <- makeCluster(n_cores)
registerDoParallel(cl)                # Register parallel backend for foreach

# ---- Method parameters ---
alpha_dvine <- 0.3                    # Target false discovery rate (FDR)
vinecop_family <- 'nonparametric'     # Type of bivariate copula used in D-vine fitting
M <- 50                               # Number of e-value repetitions
M_lasso <- 10                         # Number of Lasso repetitions per e-value
lasso_family <- 'cox'                 # Lasso model family
n_folds <- 5                          # Number of folds for cross-validated Lasso

start_time <- Sys.time()

# ---- D-vine ordering ---
# The variable order strongly affects D-vine modeling; 
# here it is obtained via a TSP-based heuristic
TSP_order <- get_dvine_order(X, random_order = TRUE)
X_dvine_order <- X[, TSP_order]

# ---- Column type identification ---
# Identify variable types (e.g., ordinal vs continuous) for later transformations
column_type <- as.vector(
  apply(X_dvine_order, MARGIN = 2, FUN = column_type_identification)
)

# ---- Ordinal-to-uniform transformation ---
# Ordinal variables are mapped to the copula scale (uniform margins)
X_cont_dvine_order <- apply(
  X = X_dvine_order, MARGIN = 2, FUN = ordinal_to_uniform
)

# ---- D-vine fitting ---
# Fit pair-copula distributions for all edges in the D-vine
dvine_distributions <- X_Xk_dvine_distributions(
  X_cont = X_cont_dvine_order, 
  vinecop_family, 
  n_cores
)

# ---- Matrix conversion ---
# Convert to matrix format required by the knockoff procedure
X_dvine_order_matrix <- as.matrix(X_dvine_order)

# ---- Derandomized knockoff procedure using e-values ---
res_dvines <- ekn_dvines(
  X = X_dvine_order_matrix, 
  X_cont = X_cont_dvine_order, 
  column_type = column_type, 
  y = y_survival, 
  dvine_distributions, 
  M = M, 
  M_lasso = M_lasso,
  alpha = alpha_dvine, 
  gamma = alpha_dvine / 2, 
  lasso_family = lasso_family, 
  n_cores = n_cores,
  n_folds = n_folds, 
  random_cv = TRUE, 
  random_evalues = TRUE
)

end_time <- Sys.time()
nonDTLDCKe_time <- end_time - start_time
print("Time for running the non-DTLDCKe procedure:")
print(nonDTLDCKe_time)

stopCluster(cl)                       # Shut down the parallel cluster

# ---- Extract selected variables ---
# Indices of variables rejected as null (i.e., selected)
just_rejections_nonpar_dvines <- sort(res_dvines$rej)
print(paste0(
  "The number of selected variables is: ", 
  length(just_rejections_nonpar_dvines)
))

# Binary rejection indicator for all variables (0 = null, 1 = selected)
rejections_nonpar_dvines <- rep(0, p)
rejections_nonpar_dvines[just_rejections_nonpar_dvines] <- 1

# Names of selected variables
selected_columns_nonpar_dvines <- 
  colnames(X_dvine_order)[rejections_nonpar_dvines == 1]
print(selected_columns_nonpar_dvines)


#### Parametric DTLDCKe (D-vines) #####

# ---- Parallel computing setup ---
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# ---- Method parameters ---
alpha_dvine <- 0.35                   # Target false discovery rate (FDR)
vinecop_family <- 'parametric'        # Parametric pair-copula family
M <- 50
M_lasso <- 10
lasso_family <- 'cox'
n_folds <- 5

start_time <- Sys.time()

# ---- D-vine ordering ---
# Variable ordering is obtained via a randomized TSP-based procedure
TSP_order <- get_dvine_order(X, random_order = TRUE)
X_dvine_order <- X[, TSP_order]

# ---- Column type identification ---
column_type <- as.vector(
  apply(X_dvine_order, MARGIN = 2, FUN = column_type_identification)
)

# ---- Ordinal-to-uniform transformation ---
X_cont_dvine_order <- apply(
  X = X_dvine_order, MARGIN = 2, FUN = ordinal_to_uniform
)

# ---- D-vine fitting ---
# Fit parametric pair-copulas with a minimum dependence threshold (psi0)
dvine_distributions <- X_Xk_dvine_distributions(
  X_cont = X_cont_dvine_order, 
  vinecop_family, 
  n_cores, 
  psi0 = 0.95
)

# ---- Matrix conversion ---
X_dvine_order_matrix <- as.matrix(X_dvine_order)

# ---- Derandomized knockoff procedure using e-values ---
res_dvines <- ekn_dvines(
  X = X_dvine_order_matrix, 
  X_cont = X_cont_dvine_order, 
  column_type = column_type, 
  y = y_survival, 
  dvine_distributions, 
  M = M, 
  M_lasso = M_lasso,
  alpha = alpha_dvine, 
  gamma = alpha_dvine / 2, 
  lasso_family = lasso_family, 
  n_cores = n_cores,
  n_folds = n_folds, 
  random_cv = TRUE, 
  random_evalues = TRUE
)

end_time <- Sys.time()
DTLDCKe_time <- end_time - start_time
print("Time for running the parametric DTLDCKe procedure:")
print(DTLDCKe_time)

stopCluster(cl)

# ---- Extract selected variables ---
just_rejections_par_dvines <- sort(res_dvines$rej)
print(paste0(
  "The number of selected variables is: ", 
  length(just_rejections_par_dvines)
))

rejections_par_dvines <- rep(0, p)
rejections_par_dvines[just_rejections_par_dvines] <- 1

selected_columns_par_dvines <- 
  colnames(X_dvine_order)[rejections_par_dvines == 1]
print(selected_columns_par_dvines)

#### Non-parametric DTLCCKe (C-vines) #####

# ---- Parallel computing setup --
n_cores <- detectCores() - 1           # Number of cores used for parallel execution
cl <- makeCluster(n_cores)
registerDoParallel(cl)                # Register parallel backend for foreach

# ---- Method parameters ---
alpha_cvine <- 0.3                    # Target false discovery rate (FDR)
vinecop_family <- 'nonparametric'     # Type of bivariate copula used in C-vine fitting
M <- 50                               # Number of e-value repetitions
M_lasso <- 10                         # Number of Lasso repetitions per e-value
lasso_family <- 'cox'                 # Lasso model family
n_folds <- 5                          # Number of folds for cross-validated Lasso

start_time <- Sys.time()

# ---- C-vine ordering ---
# Root node selection and variable ordering defining the C-vine structure
cvine_order <- get_cvine_order(X)
X_cvine_order <- X[, cvine_order]

# ---- Column type identification ---
# Identify variable types for mixed-data handling
column_type <- as.vector(
  apply(X_cvine_order, MARGIN = 2, FUN = column_type_identification)
)

# ---- Ordinal-to-uniform transformation ---
# Map ordinal variables to uniform margins for copula modeling
X_cont_cvine_order <- apply(
  X = X_cvine_order, MARGIN = 2, FUN = ordinal_to_uniform
)

# ---- C-vine fitting ---
# Fit pair-copula distributions for all edges in the C-vine
rvine_distributions <- X_Xk_rvine_distributions(
  X_cont = X_cont_cvine_order, 
  vinecop_family, 
  n_cores
)

# ---- Matrix conversion ---
# Convert data to matrix format required by the knockoff procedure
X_cvine_order_matrix <- as.matrix(X_cvine_order)

# ---- Derandomized knockoff procedure using e-values ---
res_cvines <- ekn_cvines(
  X = X_cvine_order_matrix, 
  X_cont = X_cont_cvine_order, 
  column_type = column_type, 
  y = y_survival, 
  rvine_distributions, 
  M = M, 
  M_lasso = M_lasso,
  alpha = alpha_cvine, 
  gamma = alpha_cvine / 2, 
  lasso_family = lasso_family, 
  n_cores = n_cores, 
  random_cv = FALSE, 
  random_evalues = FALSE
)

end_time <- Sys.time()
nonDTLCCKe_time <- end_time - start_time
print("Time for running the nonparametric DTLCCKe procedure:")
print(nonDTLCCKe_time)

stopCluster(cl)                       # Shut down the parallel cluster

# ---- Extract selected variables ---
# Indices of variables rejected as null (i.e., selected)
just_rejections_nonpar_cvines <- sort(res_cvines$rej)
print(paste0(
  "The number of selected variables is: ", 
  length(just_rejections_nonpar_cvines)
))

# Binary rejection indicator for all variables (0 = null, 1 = selected)
rejections_nonpar_cvines <- rep(0, p)
rejections_nonpar_cvines[just_rejections_nonpar_cvines] <- 1

# Names of selected variables
selected_columns_nonpar_cvines <- 
  colnames(X_cvine_order)[rejections_nonpar_cvines == 1]
print(selected_columns_nonpar_cvines)


#### Parametric DTLCCKe (C-vines) #####

# ---- Parallel computing setup ---
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# ---- Method parameters ---
alpha_cvine <- 0.3                    # Target false discovery rate (FDR)
vinecop_family <- 'parametric'        # Parametric pair-copula family
M <- 50
M_lasso <- 10
lasso_family <- 'cox'
n_folds <- 5

start_time <- Sys.time()

# ---- C-vine ordering ---
# Root node selection and ordering defining the C-vine structure
cvine_order <- get_cvine_order(X)
X_cvine_order <- X[, cvine_order]

# ---- Column type identification ---
column_type <- as.vector(
  apply(X_cvine_order, MARGIN = 2, FUN = column_type_identification)
)

# ---- Ordinal-to-uniform transformation ---
X_cont_cvine_order <- apply(
  X = X_cvine_order, MARGIN = 2, FUN = ordinal_to_uniform
)

# ---- C-vine fitting ---
# Fit parametric pair-copulas for the C-vine structure
rvine_distributions <- X_Xk_rvine_distributions(
  X_cont = X_cont_cvine_order, 
  vinecop_family, 
  n_cores
)

# ---- Matrix conversion ---
X_cvine_order_matrix <- as.matrix(X_cvine_order)

# ---- Derandomized knockoff procedure using e-values ---
res_cvines <- ekn_cvines(
  X = X_cvine_order_matrix, 
  X_cont = X_cont_cvine_order, 
  column_type = column_type, 
  y = y_survival, 
  rvine_distributions, 
  M = M, 
  M_lasso = M_lasso,
  alpha = alpha_cvine, 
  gamma = alpha_cvine / 2, 
  lasso_family = lasso_family, 
  n_cores = n_cores, 
  random_cv = FALSE, 
  random_evalues = FALSE
)

end_time <- Sys.time()
DTLCCKe_time <- end_time - start_time
print("Time for running the parametric DTLCCKe procedure:")
print(DTLCCKe_time)

stopCluster(cl)

# ---- Extract selected variables ---
just_rejections_par_cvines <- sort(res_cvines$rej)
print(paste0(
  "The number of selected variables is: ", 
  length(just_rejections_par_cvines)
))

rejections_par_cvines <- rep(0, p)
rejections_par_cvines[just_rejections_par_cvines] <- 1

selected_columns_par_cvines <- 
  colnames(X_cvine_order)[rejections_par_cvines == 1]
print(selected_columns_par_cvines)


#### Model-X knockoff (second-order approximation) ####

# ---- Parallel computing setup ---
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# ---- Method parameters ---
alpha_second_order <- 0.9              # Target false discovery rate (FDR)
M <- 50                                # Number of e-value repetitions
M_lasso <- 10                          # Number of Lasso repetitions per e-value
lasso_family <- 'cox'                  # Lasso model family
n_folds <- 5                           # Number of folds for cross-validated Lasso

start_time <- Sys.time()

# ---- Second-order Model-X knockoff procedure ---
# Gaussian second-order approximation for knockoff generation
res_second_order <- ekn_second_order(
  X = X, 
  y = y_survival, 
  M = M, 
  M_lasso = M_lasso, 
  alpha = alpha_second_order, 
  gamma = alpha_second_order / 2, 
  lasso_family = lasso_family,
  n_folds = n_folds,
  n_cores = n_cores, 
  random_cv = FALSE, 
  random_evalues = FALSE
)

end_time <- Sys.time()
second_time <- end_time - start_time
print("Time for running the second-order Model-X:")
print(second_time)

stopCluster(cl)

# ---- Extract selected variables ---
just_rejections_second_order <- sort(res_second_order$rej)
print(paste0(
  "The number of selected variables is: ", 
  length(just_rejections_second_order)
))

# Binary rejection indicator for all variables
rejections_second_order <- rep(0, p)
rejections_second_order[just_rejections_second_order] <- 1

# Names of selected variables
selected_columns_second_order <- 
  colnames(X)[rejections_second_order == 1]
print(selected_columns_second_order)

#### Applying the Sequential Knockoffs ####

# ---- Parallel computing configuration ---
options(clustermq.scheduler = "multiprocess")
n_cores <- detectCores() - 2           # Number of cores reserved for computation

# ---- Method parameters ---
alpha <- 0.35                          # Target false discovery rate (FDR)
M <- 50                                # Number of knockoff repetitions

# ---- Data preparation ---
# Convert design matrix to a tibble for formula-based modeling
df_X <- as_tibble(X)

# Encode categorical binary variables as factors
df_X$sample_type   <- as.factor(df_X$sample_type)
df_X$primarysite   <- as.factor(df_X$primarysite)
df_X$summarygrade  <- as.factor(df_X$summarygrade)
df_X$summarystage  <- as.factor(df_X$summarystage)
df_X$tumorstage    <- as.factor(df_X$tumorstage)
df_X$grade         <- as.factor(df_X$grade)

start_time <- Sys.time()

# ---- Knockoff statistics computation ---
# Compute knockoff statistics for survival outcomes
W <- knockoff.statistics(
  y = y_survival, 
  X = df_X, 
  type = "survival", 
  M = M
)

# ---- Variable selection ---
# Select variables controlling the FDR at the specified level
S <- variable.selections(W, level = alpha, error.type = "fdr")
S$stable.variables

end_time <- Sys.time()
sequential_time <- end_time - start_time
print("Running time of the Sequential knockoffs procedure for the simulations:")
print(sequential_time)


#### Fitting the Cox proportional hazards model with Lasso penalization ####

# ---- Parallel computing setup ---
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)                # Register parallel backend for foreach

# ---- Method parameters ---
M_lasso <- 10                          # Number of stability-selection repetitions
lasso_family <- 'cox'                 # Cox proportional hazards model
n_folds <- 10                          # Number of folds for cross-validated Lasso

start_time <- Sys.time()

# ---- Stable Lasso fitting ---
coef <- stable_lasso_glmnet_parallel(
  X, 
  y_survival, 
  lasso_family, 
  M_lasso, 
  n_folds, 
  random_cv = TRUE
)

end_time <- Sys.time()
lasso_time <- end_time - start_time
print("Time for running the Cox Lasso Model:")
print(lasso_time)

stopCluster(cl)                       # Shut down the parallel cluster

# ---- Extract selected variables ---
# Indices of variables with nonzero estimated coefficients
filtered_by_lasso <- which(coef != 0)
print(paste0(
  "The number of selected variables is: ", 
  length(filtered_by_lasso)
))

# Binary rejection indicator for all variables
rejections_lasso <- rep(0, p)
rejections_lasso[filtered_by_lasso] <- 1

# Names of selected variables
selected_columns_lasso <- 
  colnames(X)[rejections_lasso == 1]
print(selected_columns_lasso)


#### Fitting the Cox proportional hazards model with SCAD penalization ####

# ---- Parallel computing setup ---
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# ---- Method parameters ---
penalty_type <- "SCAD"                 # Nonconvex SCAD penalty
n_folds <- 10                          # Number of folds for cross-validation

start_time <- Sys.time()

# ---- Cross-validated SCAD fitting ---
cvfit <- cv.ncvsurv(
  X = X, 
  y = y_survival, 
  penalty = penalty_type, 
  cluster = cl, 
  nfolds = n_folds
)

fit <- ncvsurv(
  X = X, 
  y = y_survival, 
  penalty = penalty_type, 
  cluster = cl, 
  lambda = cvfit$lambda.min
)

Z <- coef(fit) # No intercept term to exclude

end_time <- Sys.time()
scad_time <- end_time - start_time
print("Time for running the SCAD penalty regression model:")
print(scad_time)

stopCluster(cl)

# ---- Extract selected variables ---
# Indices of variables with nonzero SCAD coefficients
just_rejections_SCAD <- which(abs(Z) > 0)
print(paste0(
  "The number of selected variables is: ", 
  length(just_rejections_SCAD)
))

# Binary rejection indicator for all variables
rejections_SCAD <- rep(0, p)
rejections_SCAD[just_rejections_SCAD] <- 1

# Names of selected variables
selected_columns_SCAD <- 
  colnames(X)[rejections_SCAD == 1]
print(selected_columns_SCAD)