# new for cont 25.10.29

#' OLS test for no interaction for continuous X 
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param X a numeric vector of covariate values, whose
#'  treatment-covariate interaction is of interest.
#'
#'@return The two-sieded p-value for the test. 
#'  Can be either conservative or anti-conservariave
#'
#'@export
ols.test.cont <- function(Y, A, X) {
  ## OLS test for no interaction
  lm.result = lm(Y~X + A + X:A)
  summary(lm.result)$coefficients[4,4]
}

#'The usual interaction test
#'
#'Testing the interaction effect based on the simple linear model
#'  and heteroscedasticity-robust variance estimator (Huber--White)
#'
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param X a numeric vector of covariate values, whose
#'  treatment-covariate interaction is of interest.
#'
#'@return The two-sieded p-value for the test.
#'
#'@export
usual.test.cont <- function(Y, A, X) {
  ## product term t-test using Huber White variance estimators
  lm.result = lm(Y~X + A + X:A)
  vcov = vcovHC(lm.result, type = "HC")
  delta.est = summary(lm.result)$coefficients[4,1]
  var.est = vcov[4,4]
  
  t.value =  delta.est / sqrt(var.est)
  return((2 * (1- pnorm(abs(t.value))) ))
}

#'The modified interaction test
#'
#'Testing the interaction effect based on the simple linear model
#'  and modified variance estimator
#'
#'@param Y a numeric vector of observed outcomes. Its length should be the same
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
#'@return The two-sieded p-value for the test.
#'
#'@export
mod.test.cont <- function(Y, A, S, X, pi, q) {
  ## product term t-test using modified variance estimators
  n = length(Y)
  lm.result = lm(Y~X + A + X:A)
  delta.est = summary(lm.result)$coefficients[4,1]
  var.est = var.modified(Y, A, S, X, pi, q)
  
  t.value = sqrt(n) * delta.est / sqrt(var.est)
  return((2 * (1- pnorm(abs(t.value))) ))
}

#'The semiparametric efficient interaction test
#'
#'Testing the interaction effect based on the semiparametric efficient method
#'
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a numeric vector of covariate values, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param Y1_fit a numeric vector of fitted values for Y1 conditioned on covariates,
#'  expected to be derived by cross-fitting or using true Y1.
#'@param Y0_fit a numeric vector of fitted values for Y0 conditioned on covariates,
#'  expected to be derived by cross-fitting or using true Y0.
#'
#'@return The two-sieded p-value for the test.
#'
#'@export
eff.test.cont <- function(Y, A, S, X, pi, Y1_fit, Y0_fit) {
  ## semiparametric efficient interaction test with cross-fitting
  n = length(Y)
  #cross.fit.result <- Ya.fit(Y, A, X, Z)
  delta.est = delta.eff(Y, A, S, X, pi, Y1_fit, Y0_fit)
  var.est = var.eff(Y, A, S, X, pi, Y1_fit, Y0_fit)
  # return(delta.est)
  
  t.value = sqrt(n) * delta.est / sqrt(var.est)
  return((2 * (1- pnorm(abs(t.value))) ))
}
