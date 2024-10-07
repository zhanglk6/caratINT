stratify <- function(...) {
  return(as.numeric(interaction(...)))
}


#'The usual interaction test
#'
#'Testing the interaction effect based on difference in means interaction effect
#'estimator and heteroscedasticity-robust variance estimator (Huber--White)
#'
#'Testing the interaction effect based on difference in means. It
#'implements the methods as described in Sections 3.1 and 4.1, Zhang and Ma (2024).
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
#'  randomizations. Detailed information can be found in Section 2, Ma et
#'  al.(2020) and Zhang and Ma (2024).
#'
#'@return The p-value for the test.
#'
#'@references Zhang, L. & Ma, W. (2024). \emph{Interaction tests with
#'  covariate-adaptive randomization}.
#'  arXiv preprint arXiv:2311.17445.
#'
#'#'@examples
#'#The code replicates the simulation setting of Model 2 in Section 5, Zhang and Ma (2024).
#'n <- 800
#'pi <- 0.5
#'q <- pi*(1-pi)
#'X_star = runif(N, -1, 1)
#'X_d = ifelse(X_star > 0, 1, 0)
#'W = rnorm(N,0,2)
#'W_d = ifelse(W > 0, 1, 0)
#'error1 = exp(0.5*X_star)*rnorm(N)
#'error0 = 0.5*exp(0.5*X_star)*rnorm(N)
#'S <- stratify(X_d, W_d)
#'X <- X_d
#'A <- sample(c(0,1),n,replace=TRUE,prob=c(1-pi,pi))
#'alphavec <- c(5, 4, 0.5, 1.2, 2, 6)
#'Y0 <- alphavec[2] +  exp(alphavec[3]*X_star) + alphavec[5]*W + error0
#'Y1 <- alphavec[1] + exp((alphavec[3] + alphavec[4])*X_star) + alphavec[5]*W+ alphavec[6] * W *X_star + error1
#'Y <- Y0*(1-A)+Y1*A
#'usual.test(Y, A, S, X, pi, q)
#'@export
usual.test <- function(Y, A, S, X, pi, q) {
  A = cat.to.int(A)
  S = cat.to.int(S) + 1
  X = cat.to.int(X)

  if (max(X) == 1) {
    return(old.test(Y, A, S, X, pi, q))
  } else {
      return(old.test.multi(Y, A, S, X, pi, q))
    }
}

modified.test <- function(Y, A, S, X, pi, q) {
  A = cat.to.int(A)
  S = cat.to.int(S) + 1
  X = cat.to.int(X)

  if (max(X) == 1) {
    return(new.test(Y, A, S, X, pi, q))
  } else {
    return(new.test.multi(Y, A, S, X, pi, q))
  }
}

stratified.adjusted.test <- function(Y, A, S, X, pi, q) {
  A = cat.to.int(A)
  S = cat.to.int(S) + 1
  X = cat.to.int(X)

  if (max(X) == 1) {
    return(strata.test(Y, A, S, X, pi, q))
  } else {
    return(strata.test.multi(Y, A, S, X, pi, q))
  }
}


