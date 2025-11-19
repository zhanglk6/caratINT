random.schemes <- function(N, pi, S, profiles) {
  SRS_A = SRS(N, pi)
  SBR_A = SBR(S, pi, bsize = 6)
  SBCD_A = SBCD(S, pi)
  PS_A = PS(profiles, rep(1/ncol(profiles), ncol(profiles)), pi, lambda = 0.75)
  return(list(SRS_A = SRS_A, SBR_A = SBR_A, SBCD_A = SBCD_A, PS_A = PS_A))
}

all_outputs <- function(pi, sample_data) {
  S <- sample_data[["S"]]
  X <- sample_data[["X"]]
  
  SRS_Y = sample_data[["SRS_Y"]]
  SRS_A = sample_data[["SRS_A"]]
  SRS_Y1_fit = sample_data[["SRS_Y1_fit"]]
  SRS_Y0_fit = sample_data[["SRS_Y0_fit"]]
  
  SBR_Y = sample_data[["SBR_Y"]]
  SBR_A = sample_data[["SBR_A"]]
  SBR_Y1_fit = sample_data[["SBR_Y1_fit"]]
  SBR_Y0_fit = sample_data[["SBR_Y0_fit"]]
  
  SBCD_Y = sample_data[["SBCD_Y"]]
  SBCD_A = sample_data[["SBCD_A"]]
  SBCD_Y1_fit = sample_data[["SBCD_Y1_fit"]]
  SBCD_Y0_fit = sample_data[["SBCD_Y0_fit"]]
  
  PS_Y = sample_data[["PS_Y"]]
  PS_A = sample_data[["PS_A"]]
  PS_Y1_fit = sample_data[["PS_Y1_fit"]]
  PS_Y0_fit = sample_data[["PS_Y0_fit"]]
  
  outputs = c(
  output(SRS_Y, SRS_A, S, X, pi, pi*(1-pi), SRS_Y1_fit, SRS_Y0_fit),
  output(SBR_Y, SBR_A, S, X, pi, 0, SBR_Y1_fit, SBR_Y0_fit),
  output(SBCD_Y, SBCD_A, S, X, pi, 0, SBCD_Y1_fit, SBCD_Y0_fit),
  output(PS_Y, PS_A, S, X, pi, 0, PS_Y1_fit, PS_Y0_fit)
  )
  
  return(outputs)
}

#' p-values for all methods (ols, HW, modifed, and efficient) produced by one copy of data
#'   suppose that Y1_fit and Y0_fit are fitted conditional on X and Z
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects;
#'@param X a numeric vector of covariate values, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'
#'@param Y1_fit a numeric vector of fitted values for Y1 conditioned on covariates,
#'  expected to be derived by cross-fitting or using true Y1.
#'@param Y0_fit a numeric vector of fitted values for Y0 conditioned on covariates,
#'  expected to be derived by cross-fitting or using true Y0.
#'
#'@return a vector of p-values produced by different methods.
output <- function(Y, A, S, X, pi, q, Y1_fit, Y0_fit) {
  c(
    ols.test.cont(Y, A, X),
    usual.test.cont(Y, A, X), 
    mod.test.cont(Y, A, S, X, pi, q), 
    eff.test.cont(Y, A, S, X, pi, Y1_fit, Y0_fit)
  )
}

#' Generate fold indexes of size n used for cross-fitting;
#'   first M-1 folds are of equal size
#'@param n an Integer for sample size
#'@param M an Integer for number of folds for M-fold cross-fitting 
#'
#'@return a vector of p-values produced by different methods.
create.fold <- function(n, M = 2) {
  n_fold = numeric(M)
  fold_seq = c()
  for (m in 1:(M-1)) {
    n_fold[m] = floor(n / M)
    fold_seq = c(fold_seq, rep(m, n_fold[m]))
  }
  n_fold[M] = n - sum(n_fold[1:(M-1)])
  fold_seq = c(fold_seq, rep(M, n_fold[M]))
  fold = sample(fold_seq, n)
  return(fold)
}