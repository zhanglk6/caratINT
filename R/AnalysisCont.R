# new for cont 25.10.29
ols.test.cont <- function(Y, A, X) {
  ## OLS test for no interaction
  lm.result = lm(Y~X + A + X:A)
  summary(lm.result)$coefficients[4,4]
}

usual.test.cont <- function(Y, A, X) {
  ## product term t-test using Huber White variance estimators
  lm.result = lm(Y~X + A + X:A)
  vcov = vcovHC(lm.result, type = "HC")
  delta.est = summary(lm.result)$coefficients[4,1]
  var.est = vcov[4,4]
  
  t.value =  delta.est / sqrt(var.est)
  return((2 * (1- pnorm(abs(t.value))) ))
}


mod.test.cont <- function(Y, A, S, X, pi, q) {
  ## product term t-test using modified variance estimators
  n = length(Y)
  lm.result = lm(Y~X + A + X:A)
  delta.est = summary(lm.result)$coefficients[4,1]
  var.est = var.modified(Y, A, S, X, pi, q)
  
  t.value = sqrt(n) * delta.est / sqrt(var.est)
  return((2 * (1- pnorm(abs(t.value))) ))
}

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
