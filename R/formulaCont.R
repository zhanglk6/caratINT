#' The modified variance estimator for interaction effect for continuous X
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects;
#'@param X a numeric vector of covariate values, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return The modified variance estimator for interaction effect.
var.modified <- function(Y, A, S, X, pi, q) {
  n = length(Y)
  n1 = sum(A)
  n0 = sum(1-A)
  
  meanX = mean(X)
  sigmaXSq = var(X)
  meanY1 = sum(A * Y) / n1
  meanY0 = sum((1-A) * Y) / n0
  
  r = transform.outcome(Y, A, S, X, pi, q)
  meanr1 = sum(A * r) / n1
  meanr0 = sum((1-A) * r) / n0
  
  zeta.tilde.r = zeta.tilde.r(r, A, S, X, pi, q)
  zeta.H.r = zeta.H.r(r, A, S, X, pi, q)
  zeta.A.r = zeta.A.r(r, A, S, X, pi, q)
 
  ## strong balance without the third term
  var.est = 1/sigmaXSq^2 * (zeta.tilde.r + zeta.H.r + zeta.A.r) 
  return(var.est)
}

#' The transformed outcomes used for interaction tests for continuous X
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects;
#'@param X a numeric vector of covariate values, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return A vector of the transformed outcomes r_i = r_i(1)A_i + r_i(0)(1-A_i).
transform.outcome <- function(Y, A, S, X, pi, q) {
  n = length(Y)
  n1 = sum(A)
  n0 = sum(1-A)
  
  meanX = mean(X)
  meanY1 = sum(A * Y) / n1
  meanY0 = sum((1-A) * Y) / n0
  
  index1 = which(A==1)
  index0 = which(A==0)
  
  total_data = data.frame(Y = Y, X = X)
  treated_data = total_data[index1,]
  control_data = total_data[index0,]
  
  lm.result.1 = lm(Y~X, data = treated_data)
  lm.result.0 = lm(Y~X, data = control_data)
  #delta.est = summary(lm.result)$coefficients[4,1]
  gamma1.est = summary(lm.result.1)$coefficients[2,1]
  gamma0.est = summary(lm.result.0)$coefficients[2,1]
  
  r = (X - meanX) * (Y - ( A*meanY1 + (1-A)*meanY0)) - (X - meanX)^2 * ( A*gamma1.est + (1-A)*gamma0.est)
  return(r)
}

#' The stratified mean under treated
#'@param r a numeric vector of transformed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'
#'@return A vector of the stratified means of r under treated.
mu1.cont <- function(r, A, S) {
  nStrata = max(S)
  mu1s = numeric(nStrata)
  for (i in 1:nStrata) {
    mu1s[i] = sum(r*A*ifelse(S == i, 1, 0)) / sum(A*ifelse(S == i, 1, 0))
  }
  return(mu1s)
}

#' The stratified mean under control
#'@param r a numeric vector of transformed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'
#'@return A vector of the stratified means of r under control.
mu0.cont <- function(r, A, S) {
  nStrata = max(S)
  mu0s = numeric(nStrata)
  for (i in 1:nStrata) {
    mu0s[i] = sum(r*(1-A)*ifelse(S == i, 1, 0)) / sum((1-A)*ifelse(S == i, 1, 0))
  }
  return(mu0s)
}

#' The first component in the asymptotic variance of the modified method
#'@param r a numeric vector of transformed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a numeric vector of covariate values, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return A numeric value for the first component.
zeta.tilde.r <- function(r, A, S, X, pi, q) {
  n = length(r)
  n1 = sum(A)
  n0 = sum(1-A)
  
  meanr1 = sum(A * r) / n1
  meanr0 = sum((1-A) * r) / n0
  
  nStrata = max(S)
  ns = numeric(nStrata)
  indicate.s = matrix(0, nrow = n, ncol = nStrata)
  for (i in 1:nStrata) {
    indicate.s[,i] = ifelse(S == i, 1, 0)
    ns[i] = sum(indicate.s[,i])
  }
  nsN = ns / n
  
  mu1s = mu1.cont(r, A, S)
  mu0s = mu0.cont(r, A, S)
  
  summ = 1 / pi * (1/n1 * sum(r^2*A) - sum(nsN * mu1s * mu1s)) + 1 / (1-pi) * (1/n0 * sum(r^2*(1-A)) - sum(nsN * mu0s * mu0s))
  return(summ)
}

#' The second component in the asymptotic variance of the modified method
#'@param r a numeric vector of transformed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a numeric vector of covariate values, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return A numeric value for the second component.
zeta.H.r <- function(r, A, S, X, pi, q) {
  n = length(r)
  n1 = sum(A)
  n0 = sum(1-A)
  
  meanr1 = sum(A * r) / n1
  meanr0 = sum((1-A) * r) / n0
  
  nStrata = max(S)
  ns = numeric(nStrata)
  indicate.s = matrix(0, nrow = n, ncol = nStrata)
  for (i in 1:nStrata) {
    indicate.s[,i] = ifelse(S == i, 1, 0)
    ns[i] = sum(indicate.s[,i])
  }
  nsN = ns / n
  
  mu1s = mu1.cont(r, A, S)
  mu0s = mu0.cont(r, A, S)
  
  summ = sum(nsN * (mu1s - meanr1 - mu0s + meanr0) * (mu1s - meanr1 - mu0s + meanr0))
  return(summ)
}

#' The third component in the asymptotic variance of the modified method
#'@param r a numeric vector of transformed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a numeric vector of covariate values, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return A numeric value for the third component.
zeta.A.r <- function(r, A, S, X, pi, q) {
  n = length(r)
  n1 = sum(A)
  n0 = sum(1-A)
  
  meanr1 = sum(A * r) / n1
  meanr0 = sum((1-A) * r) / n0
  
  nStrata = max(S)
  ns = numeric(nStrata)
  indicate.s = matrix(0, nrow = n, ncol = nStrata)
  for (i in 1:nStrata) {
    indicate.s[,i] = ifelse(S == i, 1, 0)
    ns[i] = sum(indicate.s[,i])
  }
  nsN = ns / n
  
  mu1s = mu1.cont(r, A, S)
  mu0s = mu0.cont(r, A, S)
  
  summ = sum(q* nsN * (1/pi * (mu1s - meanr1) + 1/(1-pi) * (mu0s - meanr0))^2 )
  
  return(summ)
}

#' The semiparametric efficient estimator for interaction effect for continuous X,
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
#'@return The semiparametric efficient estimator for interaction effect.
delta.eff <- function(Y, A, S, X, pi, Y1_fit, Y0_fit) {
  n = length(Y)
  meanX = mean(X)
  
  nStrata = max(S)
  ns = numeric(nStrata)
  n1s = numeric(nStrata)
  n0s = numeric(nStrata)
  indicate.s = matrix(0, nrow = n, ncol = nStrata)
  for (i in 1:nStrata) {
    indicate.s[,i] = ifelse(S == i, 1, 0)
    ns[i] = sum(indicate.s[,i])
    n1s[i] = sum(A*indicate.s[,i])
    n0s[i] = sum((1-A)*indicate.s[,i])
  }
  
  
  pi = n1s / ns
  nsN = ns / n
  
  mu1s = mu1.cont(Y, A, S)
  mu0s = mu0.cont(Y, A, S)
  mu1_mean = sum(nsN * mu1s)
  mu0_mean = sum(nsN * mu0s)
  
  #fitted_model = Ya.fit(Y, A, S, X, Z)
  Y_fit = A*Y1_fit + (1-A)*Y0_fit 
  res = Y - Y_fit
  #Y1_fit = fitted_model$Y1_fit
  #Y0_fit = fitted_model$Y0_fit
  Tns = 1/n * sum(A*(X-meanX)*res / pi) - 1/n * sum((1-A)*(X-meanX)*res / (1-pi)) +
    1/n * sum((X-meanX) * ((Y1_fit - Y0_fit) - (mu1_mean - mu0_mean) ))
  return(Tns / var(X))
}

#' The semiparametric efficieny bound estimator for interaction effect for continuous X,
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
#'@return The semiparametric efficiency bound estimator for interaction effect.
var.eff <- function(Y, A, S, X, pi, Y1_fit, Y0_fit) {
  n = length(Y)
  meanX = mean(X)
  
  #fitted_model = Ya.fit(Y, A,S,  X, Z)
  #res = fitted_model$res
  #Y1_fit = fitted_model$Y1_fit
  #Y0_fit = fitted_model$Y0_fit
  Y_fit = A*Y1_fit + (1-A)*Y0_fit 
  res = Y - Y_fit
  
  delta0 = delta.eff(Y, A, S, X, pi, Y1_fit, Y0_fit)
  
  nStrata = max(S)
  ns = numeric(nStrata)
  n1s = numeric(nStrata)
  n0s = numeric(nStrata)
  indicate.s = matrix(0, nrow = n, ncol = nStrata)
  for (i in 1:nStrata) {
    indicate.s[,i] = ifelse(S == i, 1, 0)
    ns[i] = sum(indicate.s[,i])
    n1s[i] = sum(A*indicate.s[,i])
    n0s[i] = sum((1-A)*indicate.s[,i])
  }
  nsN = ns / n
  mu1s = mu1.cont(Y, A, S)
  mu0s = mu0.cont(Y, A, S)
  mu1_mean = sum(nsN * mu1s)
  mu0_mean = sum(nsN * mu0s)
  
  pi = n1s / ns
  
  res_1s = numeric(nStrata)
  res_0s = numeric(nStrata)
  for (i in 1:nStrata) {
    res_1s[i] = sum((X-meanX)^2*res^2*A*indicate.s[,i]) / n1s[i]
    res_0s[i] = sum((X-meanX)^2*res^2*(1-A)*indicate.s[,i]) / n0s[i]
  }
  
  var_part1 = sum(nsN/pi * res_1s)
  var_part2 = sum(nsN/(1-pi) * res_0s)
  var_part3 = mean(((X-meanX) * ((Y1_fit - Y0_fit) - (mu1_mean - mu0_mean)) - (X-meanX)^2*delta0)^2)
  sum_of_var = var_part1 + var_part2 + var_part3
  #sum_of_var = var_part3
  return(sum_of_var  / ( var(X)^2) )
  
}