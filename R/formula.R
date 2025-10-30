#'Transform categorical variables to Integers 0, 1, 2, ...
#'
#'@param categorical_var a categorical variable, e.g, covariate X and stratum S.
#'
#'@return The integer froms of the variable.
cat.to.int <- function(categorical_var) {
  # Ensure the input is a factor (if not already)
  categorical_var <- as.factor(categorical_var)

  # Convert the levels to integers starting from 0
  encoded_var <- as.integer(factor(categorical_var, levels = unique(categorical_var))) - 1

  # Return the new encoded variable
  return(encoded_var)
}

#'Generate contrast matrices for interaction tests 
#'for no heterogeneity for over all levels of X.
#'
#'@param K equals the levels of categorical X minus 1, that is, K = max(cat.to.int(X)).
#'
#'@return The contrast matrix R.
generate.R.matrix <- function(K) {
  ## contrast matrix R
  matrix(cbind(-rep(1,K), diag(K)), nrow = K)
}

#'Generate a batch of frequently used quantities to avoid redundancy.
#'
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations. not used here.
#'
#'@return some frequently used quantities
base.est.multi <- function(Y, A, S, X, pi, q) {
  ## helper for multilevel calculation
  n = length(Y)
  K = max(X)
  indicate.x = matrix(0, nrow = n, ncol = K+1)
  indicate1x = matrix(0, nrow = n, ncol = K+1)
  indicate0x = matrix(0, nrow = n, ncol = K+1)
  nx = numeric(K+1)
  n1x = numeric(K+1)
  n0x = numeric(K+1)
  meanY1x = numeric(K+1)
  meanY0x = numeric(K+1)
  varY1x = numeric(K+1)
  varY0x = numeric(K+1)
  for (i in 0:K) {
    indicate.x[,i+1] = ifelse(X == i, 1, 0)
    indicate1x[,i+1] = indicate.x[,i+1] * A
    indicate0x[,i+1] = indicate.x[,i+1] * (1-A)
    nx[i+1] = sum(indicate.x[,i+1])
    n1x[i+1] = sum(indicate1x[,i+1])
    n0x[i+1] = sum(indicate0x[,i+1])
    meanY1x[i+1] = sum(indicate1x[,i+1] * Y) / n1x[i+1]
    meanY0x[i+1] = sum(indicate0x[,i+1] * Y) / n0x[i+1]
    varY1x[i+1] = sum((Y - meanY1x[i+1])^2 * indicate1x[,i+1]) / n1x[i+1]
    varY0x[i+1] = sum((Y - meanY0x[i+1])^2 * indicate0x[,i+1]) / n0x[i+1]
  }
  px = nx / n
  return(list(n = n, K = K, indicate.x =  indicate.x, indicate1x = indicate1x,
              indicate0x = indicate0x, nx = nx, n1x = n1x, n0x = n0x,
              meanY1x = meanY1x, meanY0x = meanY0x, varY1x =varY1x,
              varY0x = varY0x, px = px))
}

#' OLS test for no interaction for binary X 
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects, not used here.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum, not used here.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations. not used here.
#'
#'@return The two-sieded p-value for the test. 
#'  Can be either conservative or anti-conservariave
ols.test <- function(Y, A, S, X, pi, q) {
  ## OLS test for no interaction
  lm.result = lm(Y~X + A + X:A)
  summary(lm.result)$coefficients[4,4]
}

#' Usual test using Huber White sandwitch variance estimator for no interaction for binary X 
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return The two-sieded p-value for the test. 
#'  Can be conservative
old.test <- function(Y, A, S, X, pi, q) {
  ## product term t-test using Huber White variance estimators
  n = length(Y)
  tau.est = tau.t(Y, A, S, X, pi, q)
  var.est = var.heter(Y, A, S, X, pi, q)

  t.value = sqrt(n) * tau.est / sqrt(var.est)
  return((2 * (1- pnorm(abs(t.value))) ))
}

#' Usual test using Huber White sandwitch variance estimator for no interaction for categorical X 
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return The two-sieded p-value for the test. 
#'  Can be conservative
old.test.multi <- function(Y, A, S, X, pi, q) {
  ## product term Wald test using Huber White variance estimators
  K = max(X)
  n = length(Y)
  tau.est = tau.t.multi(Y, A, S, X, pi, q)
  var.inv.est = var.heter.multi.inv(Y, A, S, X, pi, q)

  wald.value = n * t(tau.est) %*% var.inv.est %*% tau.est
  return(1- pchisq(q = wald.value, df = K) )
}

#' Modified test for no interaction for binary X 
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return The two-sieded p-value for the test. 
#'  Asymptotically exact
new.test <- function(Y, A, S, X, pi, q) {
  ## product term t-test using modified variance estimators
  n = length(Y)
  tau.est = tau.t(Y, A, S, X, pi, q)
  var.est = var.modified(Y, A, S, X, pi, q)

  t.value = sqrt(n) * tau.est / sqrt(var.est)
  return((2 * (1- pnorm(abs(t.value))) ))
}

#' Modified test for no interaction for categorical X 
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return The two-sieded p-value for the test. 
#'  Asymptotically exact
new.test.multi <- function(Y, A, S, X, pi, q) {
  ## product term Wald test using Huber White variance estimators
  K = max(X)
  n = length(Y)
  tau.est = tau.t.multi(Y, A, S, X, pi, q)
  var.inv.est = var.modified.multi.inv(Y, A, S, X, pi, q)

  wald.value = n * t(tau.est) %*% var.inv.est %*% tau.est
  return(1- pchisq(q = wald.value, df = K) )
}

#' Stratified-adjusted test for no interaction for binary X 
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations; actually not used here.
#'
#'@return The two-sieded p-value for the test. 
#'  Asymptotically exact and more powerful
strata.test <- function(Y, A, S, X, pi, q) {
  ## strata combination test
  n = length(Y)
  tau.est = tau.strata(Y, A, S, X, pi, q)
  var.est = var.strata(Y, A, S, X, pi, q)

  t.value = sqrt(n) * tau.est / sqrt(var.est)
  return((2 * (1- pnorm(abs(t.value))) ))
}

#' Stratified-adjusted test for no interaction for categorical X 
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations; actually not used here.
#'
#'@return The two-sieded p-value for the test. 
#'  Asymptotically exact and more powerful
strata.test.multi <- function(Y, A, S, X, pi, q) {
  ## product term Wald test using Huber White variance estimators
  K = max(X)
  n = length(Y)
  tau.est = tau.strata.multi(Y, A, S, X, pi, q)
  var.inv.est = var.strata.multi.inv(Y, A, S, X, pi, q)

  wald.value = n * t(tau.est) %*% var.inv.est %*% tau.est
  return(1- pchisq(q = wald.value, df = K) )
}

#' Simple difference-in-difference estimate for interaction effect for binary X
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects; not used here.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum; 
#'  not used here.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations; not used here.
#'
#'@return The usual point estimate for interaction effect.
tau.t <- function(Y, A, S, X, pi, q) {
  ## calculate Y_{11} - Y_{10} - Y_{01} + Y_{00}
  indicate11 = A * ifelse(X == 1, 1, 0)
  indicate10 = A * ifelse(X == 0, 1, 0)
  indicate01 = (1-A) * ifelse(X == 1, 1, 0)
  indicate00 = (1-A) * ifelse(X == 0, 1, 0)

  n11 = sum(indicate11)
  n10 = sum(indicate10)
  n01 = sum(indicate01)
  n00 = sum(indicate00)

  meanY11 = sum(indicate11 * Y) / n11
  meanY10 = sum(indicate10 * Y) / n10
  meanY01 = sum(indicate01 * Y) / n01
  meanY00 = sum(indicate00 * Y) / n00

  difference = meanY11 - meanY10 - meanY01 + meanY00
  return(difference)
}

#' Simple difference-in-difference estimate for interaction effect for categorical X
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects; not used here.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum; 
#'  not used here.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations; not used here.
#'
#'@return The usual point estimate for interaction effect.
tau.t.multi <- function(Y, A, S, X, pi, q) {
  ## calculate Y_{1x} - Y_{10} - Y_{0x} + Y_{00}
  ## X = 0, 1, 2, ...
  est.package = base.est.multi(Y, A, S, X, pi, q)
  n = est.package[["n"]]
  K = est.package[["K"]]
  meanY1x = est.package[["meanY1x"]]
  meanY0x = est.package[["meanY0x"]]
  difference = numeric(K)
  for (i in 1:K) {
    difference[i] = meanY1x[i+1] - meanY0x[i+1] - meanY1x[1] + meanY0x[1]
  }
  return(difference)
}

#' Stratified-adajusted estimate for interaction effect for binary X
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum; 
#'  not used here.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations; not used here.
#'
#'@return The stratified-adjusted point estimate for interaction effect.
tau.strata <- function(Y, A, S, X, pi, q) {
  ## strata combination estimator
  n = length(Y)

  indicate.x1 = ifelse(X == 1, 1, 0)
  indicate.x0 = ifelse(X == 0, 1, 0)

  indicate1x = matrix(0, nrow = n, ncol = 2)
  indicate0x = matrix(0, nrow = n, ncol = 2)

  indicate1x[,1]= A * indicate.x0
  indicate1x[,2]= A * indicate.x1
  indicate0x[,1]= (1-A) * indicate.x0
  indicate0x[,2]= (1-A) * indicate.x1

  n10 = sum(indicate1x[,1])
  n11 = sum(indicate1x[,2])
  n00 = sum(indicate0x[,1])
  n01 = sum(indicate0x[,2])

  nDot1 = (n11 + n01)
  nDot0 = (n10 + n00)

  meanY11 = sum(indicate1x[,2] * Y) / n11
  meanY10 = sum(indicate1x[,1] * Y) / n10
  meanY01 = sum(indicate0x[,2] * Y) / n01
  meanY00 = sum(indicate0x[,1] * Y) / n00

  mu1 <- function(x,s) {
    ## \hat \mu_{1x}(s)
    #indicate1xs = A * ifelse(S == s, 1, 0) * ifelse(X == x, 1, 0)
    indicate1xs = indicate.s[,s] * indicate1x[,x+1]
    n1xs = sum(indicate1xs)
    if (n1xs == 0) {
      return(0)
    }
    sum(indicate1xs*Y) / n1xs
  }
  mu0 <- function(x,s) {
    ## \hat \mu_{0x}(s)
    #indicate0xs = (1-A) * ifelse(S == s, 1, 0) * ifelse(X == x, 1, 0)
    indicate0xs = indicate.s[,s] * indicate0x[,x+1]
    n0xs = sum(indicate0xs)
    if (n0xs == 0) {
      return(0)
    }
    sum(indicate0xs*Y) / n0xs
  }

  ## n(s)
  ## nStrata = S %>% unique() %>% length()
  nStrata = max(S)
  ns = numeric(nStrata)
  indicate.s = matrix(0, nrow = n, ncol = nStrata)
  for (i in 1:nStrata) {
    indicate.s[,i] = ifelse(S == i, 1, 0)
    ns[i] = sum(indicate.s[,i])
  }

  ## p_{s|x}
  pDot1Cond <- function(s) {
    indicateDot1s = indicate.x1 * indicate.s[,s]
    sum(indicateDot1s) / nDot1
  }
  pDot0Cond <- function(s) {
    indicateDot0s = indicate.x0  * indicate.s[,s]
    sum(indicateDot0s) / nDot0
  }

  mu11 = numeric(nStrata)
  mu01 = numeric(nStrata)
  mu10 = numeric(nStrata)
  mu00 = numeric(nStrata)

  pDot1CondS = numeric(nStrata)
  pDot0CondS = numeric(nStrata)

  for (i in 1:nStrata) {
    mu11[i] = mu1(1, i)
    mu01[i] = mu0(1, i)
    mu10[i] = mu1(0, i)
    mu00[i] = mu0(0, i)
    pDot1CondS[i] = pDot1Cond(i)
    pDot0CondS[i] = pDot0Cond(i)
  }

  difference_strata_withX =  sum(pDot1CondS * (mu11 - mu01))
  difference_strata_withoutX =  sum(pDot0CondS * (mu10 - mu00))
  difference_strata = difference_strata_withX - difference_strata_withoutX
  return(difference_strata)
}

#' Stratified-adajusted estimate for interaction effect for categorical X
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum; 
#'  not used here.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations; not used here.
#'
#'@return The stratified-adjusted point estimate for interaction effect.
tau.strata.multi <- function(Y, A, S, X, pi, q) {
  ## strata combination estimator
  est.package = base.est.multi(Y, A, S, X, pi, q)
  n = est.package[["n"]]
  K = est.package[["K"]]
  meanY1x = est.package[["meanY1x"]]
  meanY0x = est.package[["meanY0x"]]
  n1x = est.package[["n1x"]]
  n0x = est.package[["n0x"]]
  indicate.x = est.package[["indicate.x"]]
  indicate1x = est.package[["indicate1x"]]
  indicate0x = est.package[["indicate0x"]]
  px = est.package[["px"]]

  nDotx = n1x + n0x

  mu1 <- function(x,s) {
    ## \hat \mu_{1x}(s)
    #indicate1xs = A * ifelse(S == s, 1, 0) * ifelse(X == x, 1, 0)
    indicate1xs = indicate.s[,s] * indicate1x[,x+1]
    n1xs = sum(indicate1xs)
    if (n1xs == 0) {
      return(0)
    }
    sum(indicate1xs*Y) / n1xs
  }
  mu0 <- function(x,s) {
    ## \hat \mu_{0x}(s)
    #indicate0xs = (1-A) * ifelse(S == s, 1, 0) * ifelse(X == x, 1, 0)
    indicate0xs = indicate.s[,s] * indicate0x[,x+1]
    n0xs = sum(indicate0xs)
    if (n0xs == 0) {
      return(0)
    }
    sum(indicate0xs*Y) / n0xs
  }
  ## n(s)
  ## nStrata = S %>% unique() %>% length()
  nStrata = max(S)
  ns = numeric(nStrata)
  indicate.s = matrix(0, nrow = n, ncol = nStrata)
  for (i in 1:nStrata) {
    indicate.s[,i] = ifelse(S == i, 1, 0)
    ns[i] = sum(indicate.s[,i])
  }

  ## p_{s|x}
  p.s.cond.x <- function(x,s) {
    indicateDotxs = indicate.x[,x+1] * indicate.s[,s]
    sum(indicateDotxs) / nDotx[x+1]
  }

  ## p_{s|x} with x by row and s by col
  p.sx.mat = matrix(0, nrow = K+1, ncol = nStrata)
  mu1.mat = matrix(0, nrow = K+1, ncol = nStrata)
  mu0.mat = matrix(0, nrow = K+1, ncol = nStrata)
  for (k in 0:K) {
    mu1.mat[k+1,] =  sapply(1:nStrata, function(i)mu1(k,i))
    mu0.mat[k+1,] =  sapply(1:nStrata, function(i)mu0(k,i))
    p.sx.mat[k+1,] = sapply(1:nStrata, function(i)p.s.cond.x(k, i))
  }

  tau = numeric(K + 1)
  tau = apply(p.sx.mat * (mu1.mat - mu0.mat), 1, sum)
  R = generate.R.matrix(K)
  difference_strata = R %*% tau
  return(difference_strata)
}

#' Huber White variance estimator for interaction effect for binary X
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects; not used here.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum; 
#'  not used here.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations; not used here.
#'
#'@return The Huber White variance estimator for interaction effect.
var.heter <- function(Y, A, S, X, pi, q) {
  ## calculate the Huber White variance estimators
  n = length(Y)
  indicate11 = A * ifelse(X == 1, 1, 0)
  indicate10 = A * ifelse(X == 0, 1, 0)
  indicate01 = (1-A) * ifelse(X == 1, 1, 0)
  indicate00 = (1-A) * ifelse(X == 0, 1, 0)

  n11 = sum(indicate11)
  n10 = sum(indicate10)
  n01 = sum(indicate01)
  n00 = sum(indicate00)

  meanY11 = sum(indicate11 * Y) / n11
  meanY10 = sum(indicate10 * Y) / n10
  meanY01 = sum(indicate01 * Y) / n01
  meanY00 = sum(indicate00 * Y) / n00

  varY11 = sum((Y - meanY11)^2 * indicate11) / n11
  varY10 = sum((Y - meanY10)^2 * indicate10) / n10
  varY01 = sum((Y - meanY01)^2 * indicate01) / n01
  varY00 = sum((Y - meanY00)^2 * indicate00) / n00

  var.est = n* (varY11 / n11 + varY10 / n10 + varY01 / n01 + varY00 / n00)
  #var.est =  varY11 * (n/n11) + varY10 * (n/n10) + varY01 * (n/n01) + varY00 * (n/n00)
  return(var.est)
}

#' The inverse of Huber White variance estimator for interaction effect for categorical X
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects; not used here.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum; 
#'  not used here.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations; not used here.
#'
#'@return The inverse of Huber White variance estimator for R\hat{tau}.
var.heter.multi.inv <- function(Y, A, S, X, pi, q) {
  ## calculate the Huber White variance estimators
  ## X = 0, 1, 2, ...
  est.package = base.est.multi(Y, A, S, X, pi, q)
  n = est.package[["n"]]
  K = est.package[["K"]]
  varY1x = est.package[["varY1x"]]
  varY0x = est.package[["varY0x"]]
  n1x = est.package[["n1x"]]
  n0x = est.package[["n0x"]]
  sigma.n0 = n * (varY1x[1] / n1x[1] +  varY0x[1] / n0x[1])
  sigma.nk = numeric(K)
  for (i in 1:K) {
    sigma.nk[i] = n * (varY1x[i+1] / n1x[i+1] +  varY0x[i+1] / n0x[i+1])
  }

  #return(list(sigma.n0 = sigma.n0, sigma.nk = sigma.nk))
  sigma.n0.inv = 1 / sigma.n0
  sigma.nk.inv = 1 / sigma.nk
  diag(sigma.nk.inv, nrow = K) - sigma.nk.inv %*% t(sigma.nk.inv) / (sigma.n0.inv + sum(sigma.nk.inv))
}

#' The modified variance estimator for interaction effect for binary X
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects;
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return The modified variance estimator for interaction effect.
var.modified <- function(Y, A, S, X, pi, q) {
  n = length(Y)

  indicate.x1 = ifelse(X == 1, 1, 0)
  indicate.x0 = ifelse(X == 0, 1, 0)


  indicate1x = matrix(0, nrow = n, ncol = 2)
  indicate0x = matrix(0, nrow = n, ncol = 2)

  indicate1x[,1]= A * indicate.x0
  indicate1x[,2]= A * indicate.x1
  indicate0x[,1]= (1-A) * indicate.x0
  indicate0x[,2]= (1-A) * indicate.x1

  n10 = sum(indicate1x[,1])
  n11 = sum(indicate1x[,2])
  n00 = sum(indicate0x[,1])
  n01 = sum(indicate0x[,2])

  nDot1 = (n11 + n01)
  nDot0 = (n10 + n00)

  meanY11 = sum(indicate1x[,2] * Y) / n11
  meanY10 = sum(indicate1x[,1] * Y) / n10
  meanY01 = sum(indicate0x[,2] * Y) / n01
  meanY00 = sum(indicate0x[,1] * Y) / n00

  ## p_x
  p1 = nDot1 / n
  p0 = nDot0 / n

  mu1 <- function(x,s) {
    ## \hat \mu_{1x}(s)
    #indicate1xs = A * ifelse(S == s, 1, 0) * ifelse(X == x, 1, 0)
    indicate1xs = indicate.s[,s] * indicate1x[,x+1]
    n1xs = sum(indicate1xs)
    if (n1xs == 0) {
      return(0)
    }
    sum(indicate1xs*Y) / n1xs
  }
  mu0 <- function(x,s) {
    ## \hat \mu_{0x}(s)
    #indicate0xs = (1-A) * ifelse(S == s, 1, 0) * ifelse(X == x, 1, 0)
    indicate0xs = indicate.s[,s] * indicate0x[,x+1]
    n0xs = sum(indicate0xs)
    if (n0xs == 0) {
      return(0)
    }
    sum(indicate0xs*Y) / n0xs
  }

  ## n(s)
  ## nStrata = S %>% unique() %>% length()
  nStrata = max(S)
  ns = numeric(nStrata)
  indicate.s = matrix(0, nrow = n, ncol = nStrata)
  for (i in 1:nStrata) {
    indicate.s[,i] = ifelse(S == i, 1, 0)
    ns[i] = sum(indicate.s[,i])
  }

  ## p_{x|s}
  p1Cond <- function(s) {
    indicateDot1s = indicate.x1 * indicate.s[,s]
    sum(indicateDot1s) / ns[s]
  }
  p0Cond <- function(s) {
    indicateDot0s = indicate.x0 * indicate.s[,s]
    sum(indicateDot0s) / ns[s]
  }

  varY11 = sum((Y - meanY11)^2 * indicate1x[,2]) / n11
  varY10 = sum((Y - meanY10)^2 * indicate1x[,1]) / n10
  varY01 = sum((Y - meanY01)^2 * indicate0x[,2]) / n01
  varY00 = sum((Y - meanY00)^2 * indicate0x[,1]) / n00

  nsN = ns / n
  mu11 = numeric(nStrata)
  mu01 = numeric(nStrata)
  mu10 = numeric(nStrata)
  mu00 = numeric(nStrata)
  p1CondS = numeric(nStrata)
  p0CondS = numeric(nStrata)

  for (i in 1:nStrata) {
    mu11[i] = mu1(1, i)
    mu01[i] = mu0(1, i)
    mu10[i] = mu1(0, i)
    mu00[i] = mu0(0, i)
    p1CondS[i] = p1Cond(i)
    p0CondS[i] = p0Cond(i)
  }

  zetaE1 = 1/pi * (p1*varY11 - crossprod(nsN*p1CondS^2, (mu11 - meanY11)^2 )) +
    1/(1-pi) * (p1*varY01 - crossprod(nsN*p1CondS^2, (mu01 - meanY01)^2 ))

  zetaU1 = crossprod(nsN*q*p1CondS^2, ((mu11 - meanY11)/pi + (mu01 - meanY01)/(1-pi))^2)

  zetaS1 = crossprod(nsN*p1CondS^2, (mu11 - meanY11 - mu01 + meanY01)^2)

  zetaE0 = 1/pi * (p0*varY10 - crossprod(nsN*p0CondS^2, (mu10 - meanY10)^2 )) +
    1/(1-pi) * (p0*varY00 - crossprod(nsN*p0CondS^2, (mu00 - meanY00)^2 ))

  zetaU0 = crossprod(nsN*q*p0CondS^2, ((mu10 - meanY10)/pi + (mu00 - meanY00)/(1-pi))^2)

  zetaS0 = crossprod(nsN*p0CondS^2, (mu10 - meanY10 - mu00 + meanY00)^2)

  zetaE10 =  - crossprod(nsN*p1CondS*p0CondS,  1/pi*(mu11 - meanY11)*(mu10 - meanY10) +
                           1/(1-pi)*(mu01 - meanY01)*(mu00 - meanY00))

  zetaU10 = crossprod(nsN*q*p1CondS*p0CondS,
                      ((mu11 - meanY11)/pi + (mu01 - meanY01)/(1-pi))* ((mu10 - meanY10)/pi + (mu00 - meanY00)/(1-pi)))

  zetaS10 = crossprod(nsN*p1CondS*p0CondS,
                      (mu11 - meanY11 - mu01 + meanY01)* (mu10 - meanY10 - mu00 + meanY00))

  ## strong balance without the third term
  var.est = 1/p1^2 * (zetaE1 + zetaU1 + zetaS1) - 2/(p1*p0) * (zetaE10 + zetaU10 + zetaS10) +
    1/p0^2 * (zetaE0 + zetaU0 + zetaS0)
  return(var.est)
}

#' The inverse of modified variance estimator for interaction effect for categorical X
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return The inverse of modified variance estimator for R\hat{tau}.
var.modified.multi.inv <- function(Y, A, S, X, pi, q) {
  ## calculate the modified variance estimators
  ## X = 0, 1, 2, ...
  est.package = base.est.multi(Y, A, S, X, pi, q)
  n = est.package[["n"]]
  K = est.package[["K"]]
  meanY1x = est.package[["meanY1x"]]
  meanY0x = est.package[["meanY0x"]]
  varY1x = est.package[["varY1x"]]
  varY0x = est.package[["varY0x"]]
  n1x = est.package[["n1x"]]
  n0x = est.package[["n0x"]]
  indicate.x = est.package[["indicate.x"]]
  indicate1x = est.package[["indicate1x"]]
  indicate0x = est.package[["indicate0x"]]
  px = est.package[["px"]]

  mu1 <- function(x,s) {
    ## \hat \mu_{1x}(s)
    #indicate1xs = A * ifelse(S == s, 1, 0) * ifelse(X == x, 1, 0)
    indicate1xs = indicate.s[,s] * indicate1x[,x+1]
    n1xs = sum(indicate1xs)
    if (n1xs == 0) {
      return(0)
    }
    sum(indicate1xs*Y) / n1xs
  }
  mu0 <- function(x,s) {
    ## \hat \mu_{0x}(s)
    #indicate0xs = (1-A) * ifelse(S == s, 1, 0) * ifelse(X == x, 1, 0)
    indicate0xs = indicate.s[,s] * indicate0x[,x+1]
    n0xs = sum(indicate0xs)
    if (n0xs == 0) {
      return(0)
    }
    sum(indicate0xs*Y) / n0xs
  }

  ## n(s)
  ## nStrata = S %>% unique() %>% length()
  nStrata = max(S)
  ns = numeric(nStrata)
  indicate.s = matrix(0, nrow = n, ncol = nStrata)
  for (i in 1:nStrata) {
    indicate.s[,i] = ifelse(S == i, 1, 0)
    ns[i] = sum(indicate.s[,i])
  }

  ## p_{x|s}
  pCond <- function(x,s) {
    indicateDotxs = indicate.x[,x+1] * indicate.s[,s]
    sum(indicateDotxs) / ns[s]
  }
  nsN = ns / n
  ##\mu_{ax}(s) with x by row and s by col
  mu1.mat = matrix(0, nrow = K+1, ncol = nStrata)
  mu0.mat = matrix(0, nrow = K+1, ncol = nStrata)
  ## p_{x|s} with x by row and s by col
  p.mat = matrix(0, nrow = K+1, ncol = nStrata)
  for (k in 0:K) {
    mu1.mat[k+1,] =  sapply(1:nStrata, function(i)mu1(k,i))
    mu0.mat[k+1,] =  sapply(1:nStrata, function(i)mu0(k,i))
    p.mat[k+1,] = sapply(1:nStrata, function(i)pCond(k, i))
  }
  main.part.mat = matrix(0, nrow = K+1, ncol = nStrata)
  main.part.mat = p.mat * (1/pi * (sweep(mu1.mat,1, meanY1x)) + 1/(1-pi) * (sweep(mu0.mat,1, meanY0x)) )
  main.part.mat = sweep(main.part.mat, 2, FUN = "*", sqrt((pi*(1-pi)-q) * nsN))

  sigma.n = numeric(K+1)
  for (i in 0:K) {
    sigma.n[i+1] = 1 / pi * varY1x[i+1] + 1 / (1-pi) * varY0x[i+1]
  }

  Sigma.mod = diag(sigma.n * px, nrow = K+1) - main.part.mat %*% t(main.part.mat)
  Sigma.mod = Sigma.mod * (1/px) %*% t(1/px)
  R = generate.R.matrix(K)
  solve(R %*%Sigma.mod %*% t(R))
}

#' The stratified-adjusted variance estimator for interaction effect for binary X
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects;
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return The stratified-adajusted variance estimator for interaction effect.
var.strata <- function(Y, A, S, X, pi, q) {
  ## strata combination estimator
  n = length(Y)

  indicate.x1 = ifelse(X == 1, 1, 0)
  indicate.x0 = ifelse(X == 0, 1, 0)

  indicate1x = matrix(0, nrow = n, ncol = 2)
  indicate0x = matrix(0, nrow = n, ncol = 2)

  indicate1x[,1]= A * indicate.x0
  indicate1x[,2]= A * indicate.x1
  indicate0x[,1]= (1-A) * indicate.x0
  indicate0x[,2]= (1-A) * indicate.x1

  n10 = sum(indicate1x[,1])
  n11 = sum(indicate1x[,2])
  n00 = sum(indicate0x[,1])
  n01 = sum(indicate0x[,2])

  nDot1 = (n11 + n01)
  nDot0 = (n10 + n00)

  meanY11 = sum(indicate1x[,2] * Y) / n11
  meanY10 = sum(indicate1x[,1] * Y) / n10
  meanY01 = sum(indicate0x[,2] * Y) / n01
  meanY00 = sum(indicate0x[,1] * Y) / n00

  nDot1 = (n11 + n01)
  nDot0 = (n10 + n00)
  p1 = nDot1 / n
  p0 = nDot0 / n

  mu1 <- function(x,s) {
    ## \hat \mu_{1x}(s)
    #indicate1xs = A * ifelse(S == s, 1, 0) * ifelse(X == x, 1, 0)
    indicate1xs = indicate.s[,s] * indicate1x[,x+1]
    n1xs = sum(indicate1xs)
    if (n1xs == 0) {
      return(0)
    }
    sum(indicate1xs*Y) / n1xs
  }
  mu0 <- function(x,s) {
    ## \hat \mu_{0x}(s)
    #indicate0xs = (1-A) * ifelse(S == s, 1, 0) * ifelse(X == x, 1, 0)
    indicate0xs = indicate.s[,s] * indicate0x[,x+1]
    n0xs = sum(indicate0xs)
    if (n0xs == 0) {
      return(0)
    }
    sum(indicate0xs*Y) / n0xs
  }

  ## n(s)
  ## nStrata = S %>% unique() %>% length()
  nStrata = max(S)
  ns = numeric(nStrata)
  indicate.s = matrix(0, nrow = n, ncol = nStrata)
  for (i in 1:nStrata) {
    indicate.s[,i] = ifelse(S == i, 1, 0)
    ns[i] = sum(indicate.s[,i])
  }
  nsN = ns / n

  ## p_{s|x}
  pDot1Cond <- function(s) {
    indicateDot1s = indicate.x1 * indicate.s[,s]
    sum(indicateDot1s) / nDot1
  }
  pDot0Cond <- function(s) {
    indicateDot0s = indicate.x0  * indicate.s[,s]
    sum(indicateDot0s) / nDot0
  }

  varY11 = sum((Y - meanY11)^2 * indicate1x[,2]) / n11
  varY10 = sum((Y - meanY10)^2 * indicate1x[,1]) / n10
  varY01 = sum((Y - meanY01)^2 * indicate0x[,2]) / n01
  varY00 = sum((Y - meanY00)^2 * indicate0x[,1]) / n00

  mu11 = numeric(nStrata)
  mu01 = numeric(nStrata)
  mu10 = numeric(nStrata)
  mu00 = numeric(nStrata)

  pDot1CondS = numeric(nStrata)
  pDot0CondS = numeric(nStrata)
  p1CondS = numeric(nStrata)
  p0CondS = numeric(nStrata)

  for (i in 1:nStrata) {
    mu11[i] = mu1(1, i)
    mu01[i] = mu0(1, i)
    mu10[i] = mu1(0, i)
    mu00[i] = mu0(0, i)
    pDot1CondS[i] = pDot1Cond(i)
    pDot0CondS[i] = pDot0Cond(i)
  }
  p1CondS = pDot1CondS * nDot1 / ns
  p0CondS = pDot0CondS * nDot0 / ns

  zetaE1 = 1/pi * (p1*varY11 - crossprod(nsN*p1CondS, (mu11 - meanY11)^2 )) +
    1/(1-pi) * (p1*varY01 - crossprod(nsN*p1CondS, (mu01 - meanY01)^2 ))

  zetaS1 = crossprod(nsN*p1CondS, (mu11 - meanY11 - mu01 + meanY01)^2)

  zetaE0 = 1/pi * (p0*varY10 - crossprod(nsN*p0CondS, (mu10 - meanY10)^2 )) +
    1/(1-pi) * (p0*varY00 - crossprod(nsN*p0CondS, (mu00 - meanY00)^2 ))


  zetaS0 = crossprod(nsN*p0CondS, (mu10 - meanY10 - mu00 + meanY00)^2)

  var.est = 1/p1^2 * (zetaE1 + zetaS1) + 1/p0^2 * (zetaE0 + zetaS0)
  return(var.est)
}

#' The stratified-adjusted variance estimator for interaction effect for categorical X
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param S a categorical vector of stratum labels. Its length should be the same as
#'  the number of subjects;
#'@param X a categorical vector of covariate levels, whose
#'  treatment-covariate interaction is of interest.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations.
#'
#'@return The inverser of the stratified-adajusted variance estimator for R\hat{tau}.
var.strata.multi.inv <- function(Y, A, S, X, pi, q) {
  ## calculate the modified variance estimators
  ## X = 0, 1, 2, ...
  est.package = base.est.multi(Y, A, S, X, pi, q)
  n = est.package[["n"]]
  K = est.package[["K"]]
  meanY1x = est.package[["meanY1x"]]
  meanY0x = est.package[["meanY0x"]]
  varY1x = est.package[["varY1x"]]
  varY0x = est.package[["varY0x"]]
  n1x = est.package[["n1x"]]
  n0x = est.package[["n0x"]]
  indicate.x = est.package[["indicate.x"]]
  indicate1x = est.package[["indicate1x"]]
  indicate0x = est.package[["indicate0x"]]
  px = est.package[["px"]]

  mu1 <- function(x,s) {
    ## \hat \mu_{1x}(s)
    #indicate1xs = A * ifelse(S == s, 1, 0) * ifelse(X == x, 1, 0)
    indicate1xs = indicate.s[,s] * indicate1x[,x+1]
    n1xs = sum(indicate1xs)
    if (n1xs == 0) {
      return(0)
    }
    sum(indicate1xs*Y) / n1xs
  }
  mu0 <- function(x,s) {
    ## \hat \mu_{0x}(s)
    #indicate0xs = (1-A) * ifelse(S == s, 1, 0) * ifelse(X == x, 1, 0)
    indicate0xs = indicate.s[,s] * indicate0x[,x+1]
    n0xs = sum(indicate0xs)
    if (n0xs == 0) {
      return(0)
    }
    sum(indicate0xs*Y) / n0xs
  }

  ## n(s)
  ## nStrata = S %>% unique() %>% length()
  nStrata = max(S)
  ns = numeric(nStrata)
  indicate.s = matrix(0, nrow = n, ncol = nStrata)
  for (i in 1:nStrata) {
    indicate.s[,i] = ifelse(S == i, 1, 0)
    ns[i] = sum(indicate.s[,i])
  }

  ## p_{x|s}
  pCond <- function(x,s) {
    indicateDotxs = indicate.x[,x+1] * indicate.s[,s]
    sum(indicateDotxs) / ns[s]
  }
  nsN = ns / n
  ##\mu_{ax}(s) with x by row and s by col
  mu1.mat = matrix(0, nrow = K+1, ncol = nStrata)
  mu0.mat = matrix(0, nrow = K+1, ncol = nStrata)
  ## p_{x|s} with x by row and s by col
  p.mat = matrix(0, nrow = K+1, ncol = nStrata)
  for (k in 0:K) {
    mu1.mat[k+1,] =  sapply(1:nStrata, function(i)mu1(k,i))
    mu0.mat[k+1,] =  sapply(1:nStrata, function(i)mu0(k,i))
    p.mat[k+1,] = sapply(1:nStrata, function(i)pCond(k, i))
  }
  main.part.mat = matrix(0, nrow = K+1, ncol = nStrata)
  main.part.mat = sqrt(p.mat) * (1/pi * (sweep(mu1.mat,1, meanY1x)) + 1/(1-pi) * (sweep(mu0.mat,1, meanY0x)) )
  main.part.mat = sweep(main.part.mat, 2, FUN = "*", sqrt( nsN))
  long.sum = apply(main.part.mat^2,1,sum)

  sigma.n = numeric(K+1)
  for (i in 0:K) {
    sigma.n[i+1] = 1 / pi * varY1x[i+1] + 1 / (1-pi) * varY0x[i+1]
  }

  Sigma.heter.diag = (sigma.n * px - pi * (1-pi) * long.sum) / px^2

  sigma.n0.inv = 1 / Sigma.heter.diag[1]
  sigma.nk.inv = 1 / Sigma.heter.diag[-1]
  diag(sigma.nk.inv, nrow = K) - sigma.nk.inv %*% t(sigma.nk.inv) / (sigma.n0.inv + sum(sigma.nk.inv))
}
